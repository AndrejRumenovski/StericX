"""Regenerate only the prospectively admitted report; all gates are checked."""
from pathlib import Path
import runpy
if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parent/"revision_policy/publish.py"),run_name="__main__")
