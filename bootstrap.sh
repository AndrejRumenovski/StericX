#!/usr/bin/env bash
#
# StericX bootstrap: prepare a fresh checkout for building and reproduction.
#
# Idempotent. Safe to re-run. Steps:
#   1. Ensure uv (installs to ~/.local/bin if missing).
#   2. Ensure a Rust toolchain (installs rustup if missing).
#   3. Sync the pinned Python environment with the science extras.
#   4. Build the release binary for THIS machine's CPU baseline.
#   5. Smoke-test the binary so an ISA mismatch fails loudly and early.
#   6. Optionally (--with-quantum) install the checksum-pinned CREST/xTB
#      toolchain needed only for Study 003.
#
# Usage:
#   ./bootstrap.sh [--with-quantum] [--native]
#
#   --with-quantum  Also download and verify CREST 2.12 / xTB 6.4.0 (Study 003).
#   --native        Build with -C target-cpu=native (faster, non-portable).
#
set -euo pipefail

PROJECT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$PROJECT_ROOT"

WITH_QUANTUM=0
NATIVE=0
for arg in "$@"; do
  case "$arg" in
    --with-quantum) WITH_QUANTUM=1 ;;
    --native) NATIVE=1 ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "error: unknown argument: $arg" >&2; exit 2 ;;
  esac
done

log() { printf '\n\033[1m==> %s\033[0m\n' "$1"; }

# 1. uv --------------------------------------------------------------------
if ! command -v uv >/dev/null 2>&1; then
  log "Installing uv"
  curl -LsSf https://astral.sh/uv/install.sh | sh
  export PATH="$HOME/.local/bin:$PATH"
fi
command -v uv >/dev/null 2>&1 || { echo "error: uv still not on PATH; add ~/.local/bin" >&2; exit 2; }
log "uv $(uv --version)"

# 2. Rust ------------------------------------------------------------------
if ! command -v cargo >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  [ -f "$HOME/.cargo/env" ] && . "$HOME/.cargo/env"
fi
if ! command -v cargo >/dev/null 2>&1; then
  log "Installing Rust toolchain via rustup"
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
  # shellcheck disable=SC1091
  . "$HOME/.cargo/env"
fi
command -v cargo >/dev/null 2>&1 || { echo "error: cargo not found after install" >&2; exit 2; }
log "cargo $(cargo --version)"

# 3. Python environment ----------------------------------------------------
log "Syncing Python environment (uv sync --extra science)"
uv sync --extra science

# 4. Build -----------------------------------------------------------------
if [ "$NATIVE" -eq 1 ]; then
  log "Building release binary (target-cpu=native)"
  RUSTFLAGS="-C target-cpu=native" cargo build --release
else
  log "Building release binary (portable baseline)"
  cargo build --release
fi

# 5. Smoke test ------------------------------------------------------------
log "Smoke-testing the binary"
if ! ./target/release/stericx simulate --ddg 1.82 --temp 298.15 >/dev/null; then
  echo "error: the stericx binary failed to run on this machine." >&2
  echo "       If it crashed with 'Illegal instruction', rebuild without" >&2
  echo "       --native (a binary built for another CPU will not run here)." >&2
  exit 1
fi
echo "binary_smoke_test=ok"

# 6. Quantum toolchain (optional) -----------------------------------------
if [ "$WITH_QUANTUM" -eq 1 ]; then
  log "Installing quantum toolchain (CREST 2.12 / xTB 6.4.0)"
  ./install_quantum_tools.sh
fi

log "Bootstrap complete"
cat <<'EOF'
Next steps (see REPRODUCE.md for expected numbers):
  uv run --extra science python scripts/validate_stericx.py     # Sterimol vs morfeus
  uv run --extra science python studies/study_001_ni_hda.py --offline
  uv run --extra science python studies/study_002_buried_volume.py
  # Study 003 (needs ./bootstrap.sh --with-quantum first):
  uv run --extra science python scripts/prepare_quantum_data.py --mode lmo
  uv run --extra science python studies/study_003_quantum_geometry.py --no-build
EOF
