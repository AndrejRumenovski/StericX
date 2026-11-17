#!/usr/bin/env python3
"""Explicitly authorized temporary UI suspension, with independent automatic resume."""
import argparse, json, os, signal, subprocess, sys, time
from pathlib import Path

def identity(pid):
    raw = Path(f'/proc/{pid}/stat').read_text()
    fields = raw[raw.rindex(')') + 2:].split()
    return {'name': raw[raw.index('(') + 1:raw.rindex(')')],
            'state': fields[0], 'parent': int(fields[1]), 'start': fields[19]}

def same(pid, start):
    try:
        return identity(pid)['start'] == start
    except (OSError, ValueError):
        return False

def resume(pid, start):
    if same(pid, start):
        try:
            os.kill(pid, signal.SIGCONT)
            return True
        except ProcessLookupError:
            pass
    return False

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pid', type=int, required=True)
    p.add_argument('--start', required=True)
    p.add_argument('--seconds', type=int, default=1800)
    p.add_argument('--watchdog', action='store_true')
    p.add_argument('--receipt', type=Path)
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    assert 0 < a.seconds <= 3600
    if a.watchdog:
        deadline = time.monotonic() + a.seconds
        while time.monotonic() < deadline and same(a.pid, a.start):
            time.sleep(min(5, max(0, deadline - time.monotonic())))
        resume(a.pid, a.start)
        return
    assert a.receipt and not a.receipt.exists()
    command = a.command[1:] if a.command[:1] == ['--'] else a.command
    assert command
    target = identity(a.pid)
    assert target['name'] in ('warp', 'warp-terminal') and target['start'] == a.start
    assert target['state'] not in ('T', 't'), 'Do not take over an existing suspension'
    ancestor = os.getpid()
    while ancestor > 1:
        assert ancestor != a.pid, 'Do not suspend an ancestor of the measurement session'
        ancestor = identity(ancestor)['parent']
    receipt = {'scope': 'User-authorized suspension of one Warp UI process only; no descendants signalled.',
               'pid': a.pid, 'identity': target, 'maximum_seconds': a.seconds,
               'command': command, 'resumed': False}
    def save():
        a.receipt.write_text(json.dumps(receipt, indent=2) + '\n')
    save()
    watchdog = subprocess.Popen([sys.executable, str(Path(__file__).resolve()),
        '--watchdog', '--pid', str(a.pid), '--start', a.start, '--seconds', str(a.seconds)],
        start_new_session=True, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL)
    receipt['watchdog_pid'] = watchdog.pid
    def interrupted(signum, _frame):
        raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    try:
        assert same(a.pid, a.start), 'Process identity changed'
        os.kill(a.pid, signal.SIGSTOP)
        receipt['suspended_at'] = time.time(); save()
        result = subprocess.run(command)
        receipt['returncode'] = result.returncode
    finally:
        receipt['resumed'] = resume(a.pid, a.start)
        receipt['resume_attempted_at'] = time.time()
        receipt['resumed_at'] = receipt['resume_attempted_at'] if receipt['resumed'] else None
        receipt['identity_after'] = identity(a.pid) if same(a.pid, a.start) else None
        save()
        watchdog.terminate(); watchdog.wait(timeout=5)
    raise SystemExit(result.returncode)

if __name__ == '__main__':
    main()
