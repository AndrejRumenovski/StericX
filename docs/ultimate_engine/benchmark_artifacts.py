"""Compact journal references; full commands and outputs remain in each raw trial."""
import hashlib
import json
import os
import stat
import time
from pathlib import Path


def file_record(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'path': str(path), 'sha256': digest}


def compact_metrics(metrics, directory):
    # harness.measure has already retained the entire dictionary in this file.
    result = {key: value for key, value in metrics.items() if key != 'argv'}
    result['raw_metrics'] = file_record(Path(directory) / 'metrics.json')
    return result


def retain_fingerprint(value, path):
    path = Path(path)
    encoded = (json.dumps(value, sort_keys=True) + '\n').encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        with path.open('xb') as stream:
            stream.write(encoded)
    record = file_record(path)
    assert record['sha256'] == hashlib.sha256(encoded).hexdigest()
    return record


def drain_artifacts(directory):
    """Flush only this measurement's regular files/directories, outside all timers."""
    root = Path(directory).resolve()
    assert root.is_dir()
    started = time.monotonic()
    files = directories = bytes_observed = skipped = 0
    for parent, children, names in os.walk(root, topdown=False, followlinks=False):
        for name in names:
            path = Path(parent) / name
            info = path.lstat()
            if not stat.S_ISREG(info.st_mode):
                skipped += 1
                continue
            fd = os.open(path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
            files += 1
            bytes_observed += info.st_size
        fd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
        directories += 1
    return {'scope': 'Benchmark artifacts only; fsync outside every measured interval; no cache eviction.',
            'directory': str(root), 'files': files, 'directories': directories,
            'file_bytes_observed': bytes_observed, 'nonregular_files_skipped': skipped,
            'elapsed_seconds': time.monotonic() - started}
