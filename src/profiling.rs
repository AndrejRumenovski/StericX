//! Optional, allocation-free scope probes and requested-heap accounting.
//!
//! Enable `profiling`, install [`TrackingAllocator`] in the executable, and set
//! `STERICX_PROFILE_PATH` to a JSON path (or `-` for stderr). Capture a [`Session`]
//! before application work and finish it after all scopes and workers complete.
//! The uninstrumented executable remains the authority for benchmark wall time.
//! Scope times include probe overhead; summed worker elapsed time is not process
//! wall time or CPU time. Hardware counters and RSS require external tools.

use serde::Serialize;
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::io::{self, Write};
use std::marker::PhantomData;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

const MAX_THREADS: usize = 256;
const MAX_FUNCTIONS: usize = 256;
const MAX_DEPTH: usize = 64;
const INVALID_SITE: usize = usize::MAX;

static ENABLED: AtomicBool = AtomicBool::new(false);
static NEXT_THREAD: AtomicUsize = AtomicUsize::new(0);
static DROPPED_SCOPES: AtomicU64 = AtomicU64::new(0);
static FUNCTIONS: Mutex<[Option<Function>; MAX_FUNCTIONS]> = Mutex::new([None; MAX_FUNCTIONS]);
static TIMINGS: [[Timing; MAX_FUNCTIONS]; MAX_THREADS] =
    [const { [const { Timing::new() }; MAX_FUNCTIONS] }; MAX_THREADS];

#[derive(Clone, Copy)]
struct Function {
    stage: &'static str,
    name: &'static str,
}

struct Timing {
    calls: AtomicU64,
    inclusive_ns: AtomicU64,
    exclusive_ns: AtomicU64,
}

impl Timing {
    const fn new() -> Self {
        Self {
            calls: AtomicU64::new(0),
            inclusive_ns: AtomicU64::new(0),
            exclusive_ns: AtomicU64::new(0),
        }
    }
}

#[derive(Clone, Copy)]
struct Frame {
    function: usize,
    started: Instant,
    children_ns: u64,
}

struct Stack {
    thread: Option<usize>,
    frames: [Option<Frame>; MAX_DEPTH],
    depth: usize,
}

impl Stack {
    const fn new() -> Self {
        Self {
            thread: None,
            frames: [None; MAX_DEPTH],
            depth: 0,
        }
    }

    fn push(&mut self, function: usize, started: Instant) -> bool {
        if self.depth == MAX_DEPTH {
            return false;
        }
        self.frames[self.depth] = Some(Frame {
            function,
            started,
            children_ns: 0,
        });
        self.depth += 1;
        true
    }

    fn pop(&mut self, finished: Instant) -> (usize, u64, u64) {
        self.depth -= 1;
        let frame = self.frames[self.depth]
            .take()
            .expect("balanced scope stack");
        let inclusive_ns = finished.duration_since(frame.started).as_nanos() as u64;
        let exclusive_ns = inclusive_ns.saturating_sub(frame.children_ns);
        if self.depth > 0 {
            self.frames[self.depth - 1]
                .as_mut()
                .expect("parent scope exists")
                .children_ns += inclusive_ns;
        }
        (frame.function, inclusive_ns, exclusive_ns)
    }
}

thread_local! {
    static STACK: RefCell<Stack> = const { RefCell::new(Stack::new()) };
}

/// A lexical timing scope. It must be dropped on the thread that created it.
pub struct Span {
    active: bool,
    _not_send: PhantomData<*mut ()>,
}

impl Span {
    /// Called by `profile_scope!`; each call site registers only once.
    pub fn enter(
        site: &'static OnceLock<usize>,
        stage: &'static str,
        function: &'static str,
    ) -> Self {
        let mut span = Self {
            active: false,
            _not_send: PhantomData,
        };
        if !ENABLED.load(Ordering::Relaxed) {
            return span;
        }
        let function = *site.get_or_init(|| {
            let mut functions = FUNCTIONS.lock().unwrap_or_else(|error| error.into_inner());
            if let Some(index) = functions.iter().position(Option::is_none) {
                functions[index] = Some(Function {
                    stage,
                    name: function,
                });
                index
            } else {
                INVALID_SITE
            }
        });
        if function == INVALID_SITE {
            DROPPED_SCOPES.fetch_add(1, Ordering::Relaxed);
            return span;
        }
        span.active = STACK.with(|stack| {
            let mut stack = stack.borrow_mut();
            let thread = *stack
                .thread
                .get_or_insert_with(|| NEXT_THREAD.fetch_add(1, Ordering::Relaxed));
            if thread >= MAX_THREADS || !stack.push(function, Instant::now()) {
                DROPPED_SCOPES.fetch_add(1, Ordering::Relaxed);
                false
            } else {
                true
            }
        });
        span
    }
}

impl Drop for Span {
    fn drop(&mut self) {
        if !self.active {
            return;
        }
        let finished = Instant::now();
        STACK.with(|stack| {
            let mut stack = stack.borrow_mut();
            let (function, inclusive_ns, exclusive_ns) = stack.pop(finished);
            let timing = &TIMINGS[stack.thread.expect("scope has a thread")][function];
            timing.calls.fetch_add(1, Ordering::Relaxed);
            timing
                .inclusive_ns
                .fetch_add(inclusive_ns, Ordering::Relaxed);
            timing
                .exclusive_ns
                .fetch_add(exclusive_ns, Ordering::Relaxed);
        });
    }
}

static ALLOCATION_CALLS: AtomicU64 = AtomicU64::new(0);
static REALLOCATION_CALLS: AtomicU64 = AtomicU64::new(0);
static DEALLOCATION_CALLS: AtomicU64 = AtomicU64::new(0);
static FAILED_ALLOCATION_CALLS: AtomicU64 = AtomicU64::new(0);
static REQUESTED_BYTES: AtomicU64 = AtomicU64::new(0);
static LIVE_BYTES: AtomicU64 = AtomicU64::new(0);
static PEAK_BYTES: AtomicU64 = AtomicU64::new(0);

/// Counts successful requested heap operations while forwarding to `System`.
///
/// This does not measure allocator metadata, stack memory, file mappings or RSS.
/// Accounting starts at process startup, including allocations before `Session`.
/// Timing probes allocate no heap memory; session-path metadata is included, and
/// report creation happens after snapshot.
pub struct TrackingAllocator;

fn record_growth(bytes: u64) {
    let live = LIVE_BYTES.fetch_add(bytes, Ordering::Relaxed) + bytes;
    PEAK_BYTES.fetch_max(live, Ordering::Relaxed);
}

// SAFETY: each operation forwards the unchanged allocation contract to System.
// Bookkeeping uses only atomics, so it cannot recurse through this allocator.
unsafe impl GlobalAlloc for TrackingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc(layout) };
        if pointer.is_null() {
            FAILED_ALLOCATION_CALLS.fetch_add(1, Ordering::Relaxed);
        } else {
            ALLOCATION_CALLS.fetch_add(1, Ordering::Relaxed);
            REQUESTED_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
            record_growth(layout.size() as u64);
        }
        pointer
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc_zeroed(layout) };
        if pointer.is_null() {
            FAILED_ALLOCATION_CALLS.fetch_add(1, Ordering::Relaxed);
        } else {
            ALLOCATION_CALLS.fetch_add(1, Ordering::Relaxed);
            REQUESTED_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
            record_growth(layout.size() as u64);
        }
        pointer
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        unsafe { System.dealloc(pointer, layout) };
        DEALLOCATION_CALLS.fetch_add(1, Ordering::Relaxed);
        LIVE_BYTES.fetch_sub(layout.size() as u64, Ordering::Relaxed);
    }

    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let replacement = unsafe { System.realloc(pointer, layout, new_size) };
        if replacement.is_null() {
            FAILED_ALLOCATION_CALLS.fetch_add(1, Ordering::Relaxed);
        } else {
            REALLOCATION_CALLS.fetch_add(1, Ordering::Relaxed);
            REQUESTED_BYTES.fetch_add(new_size as u64, Ordering::Relaxed);
            if new_size >= layout.size() {
                record_growth((new_size - layout.size()) as u64);
            } else {
                LIVE_BYTES.fetch_sub((layout.size() - new_size) as u64, Ordering::Relaxed);
            }
        }
        replacement
    }
}

#[derive(Debug, Serialize)]
pub struct AllocationSnapshot {
    pub allocation_calls: u64,
    pub reallocation_calls: u64,
    pub deallocation_calls: u64,
    pub failed_allocation_calls: u64,
    pub requested_bytes: u64,
    pub live_requested_bytes: u64,
    pub peak_live_requested_bytes: u64,
}

/// Capture process-wide counters without allocating.
pub fn allocation_snapshot() -> AllocationSnapshot {
    AllocationSnapshot {
        allocation_calls: ALLOCATION_CALLS.load(Ordering::Relaxed),
        reallocation_calls: REALLOCATION_CALLS.load(Ordering::Relaxed),
        deallocation_calls: DEALLOCATION_CALLS.load(Ordering::Relaxed),
        failed_allocation_calls: FAILED_ALLOCATION_CALLS.load(Ordering::Relaxed),
        requested_bytes: REQUESTED_BYTES.load(Ordering::Relaxed),
        live_requested_bytes: LIVE_BYTES.load(Ordering::Relaxed),
        peak_live_requested_bytes: PEAK_BYTES.load(Ordering::Relaxed),
    }
}

#[derive(Serialize)]
struct FunctionReport {
    thread_index: usize,
    stage: &'static str,
    function: &'static str,
    calls: u64,
    inclusive_ns: u64,
    exclusive_ns: u64,
}

#[derive(Serialize)]
struct Report {
    schema_version: u32,
    timing_semantics: &'static str,
    allocation_semantics: &'static str,
    session_wall_ns: u64,
    observed_threads: usize,
    dropped_scopes: u64,
    allocations: AllocationSnapshot,
    stage_exclusive_ns: BTreeMap<&'static str, u64>,
    functions: Vec<FunctionReport>,
}

/// One diagnostic session per CLI process.
pub struct Session {
    path: PathBuf,
    started: Instant,
}

impl Session {
    /// Enable timers when `STERICX_PROFILE_PATH` is present.
    pub fn from_env() -> Option<Self> {
        let path = PathBuf::from(std::env::var_os("STERICX_PROFILE_PATH")?);
        ENABLED.store(true, Ordering::Relaxed);
        Some(Self {
            path,
            started: Instant::now(),
        })
    }

    /// Write JSON after application scopes have been dropped and workers joined.
    /// A `-` destination writes to stderr; stdout is never changed.
    pub fn finish(self) -> io::Result<()> {
        ENABLED.store(false, Ordering::Relaxed);
        let session_wall_ns = self.started.elapsed().as_nanos() as u64;
        let allocations = allocation_snapshot();
        let mut report = Report {
            schema_version: 1,
            timing_semantics: "Monotonic elapsed nanoseconds, including probe overhead. Exclusive times subtract nested scopes on the same thread. Worker totals are not process wall time or CPU time.",
            allocation_semantics: "Successful requested System heap allocations from process startup to pre-report snapshot, including realloc new sizes in requested_bytes and session-path metadata. Timing probes allocate no heap memory. Excludes report construction; not RSS or allocator usable-size accounting.",
            session_wall_ns,
            observed_threads: NEXT_THREAD.load(Ordering::Relaxed),
            dropped_scopes: DROPPED_SCOPES.load(Ordering::Relaxed),
            allocations,
            stage_exclusive_ns: [
                "file_parsing",
                "donor_bond_detection",
                "sterimol",
                "buried_volume",
                "pyramidalization",
                "conformer_processing",
                "model_inference",
                "output",
                "orchestration",
            ]
            .into_iter()
            .map(|stage| (stage, 0))
            .collect(),
            functions: Vec::new(),
        };
        let functions = FUNCTIONS.lock().unwrap_or_else(|error| error.into_inner());
        for (thread_index, timings) in TIMINGS.iter().enumerate().take(report.observed_threads) {
            for (index, timing) in timings.iter().enumerate() {
                let calls = timing.calls.load(Ordering::Relaxed);
                if calls == 0 {
                    continue;
                }
                let function = functions[index].expect("registered function");
                let exclusive_ns = timing.exclusive_ns.load(Ordering::Relaxed);
                *report.stage_exclusive_ns.entry(function.stage).or_default() += exclusive_ns;
                report.functions.push(FunctionReport {
                    thread_index,
                    stage: function.stage,
                    function: function.name,
                    calls,
                    inclusive_ns: timing.inclusive_ns.load(Ordering::Relaxed),
                    exclusive_ns,
                });
            }
        }
        drop(functions);
        let mut output: Box<dyn Write> = if self.path.as_os_str() == "-" {
            Box::new(io::stderr().lock())
        } else {
            Box::new(std::fs::File::create(&self.path)?)
        };
        serde_json::to_writer_pretty(&mut output, &report)?;
        writeln!(output)?;
        output.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn nested_scopes_partition_elapsed_time_without_double_counting() {
        let started = Instant::now();
        let mut stack = Stack::new();
        assert!(stack.push(0, started));
        assert!(stack.push(1, started + Duration::from_nanos(10)));
        assert!(stack.push(2, started + Duration::from_nanos(20)));
        assert_eq!(stack.pop(started + Duration::from_nanos(40)), (2, 20, 20));
        assert_eq!(stack.pop(started + Duration::from_nanos(60)), (1, 50, 30));
        assert_eq!(stack.pop(started + Duration::from_nanos(100)), (0, 100, 50));
        assert_eq!(stack.depth, 0);
    }

    #[test]
    fn scope_capacity_refuses_extra_frames_without_corrupting_stack() {
        let started = Instant::now();
        let mut stack = Stack::new();
        for function in 0..MAX_DEPTH {
            assert!(stack.push(function, started));
        }
        assert!(!stack.push(99, started));
        for function in (0..MAX_DEPTH).rev() {
            assert_eq!(stack.pop(started).0, function);
        }
        assert_eq!(stack.depth, 0);
    }

    #[test]
    fn allocator_preserves_data_and_counts_grow_shrink_and_free() {
        let before = allocation_snapshot();
        let allocator = TrackingAllocator;
        let layout = Layout::from_size_align(32, 8).unwrap();
        unsafe {
            let pointer = allocator.alloc_zeroed(layout);
            assert!(!pointer.is_null());
            assert!((0..32).all(|index| *pointer.add(index) == 0));
            pointer.write(73);
            let pointer = allocator.realloc(pointer, layout, 64);
            assert!(!pointer.is_null());
            assert_eq!(pointer.read(), 73);
            let pointer = allocator.realloc(pointer, Layout::from_size_align(64, 8).unwrap(), 16);
            assert!(!pointer.is_null());
            assert_eq!(pointer.read(), 73);
            allocator.dealloc(pointer, Layout::from_size_align(16, 8).unwrap());
        }
        let after = allocation_snapshot();
        assert_eq!(after.allocation_calls - before.allocation_calls, 1);
        assert_eq!(after.reallocation_calls - before.reallocation_calls, 2);
        assert_eq!(after.deallocation_calls - before.deallocation_calls, 1);
        assert_eq!(after.requested_bytes - before.requested_bytes, 112);
        assert_eq!(after.live_requested_bytes, before.live_requested_bytes);
    }
}
