"""Compact journal references; full commands and outputs remain in each raw trial."""
import hashlib
import json
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
