import sys
from pathlib import Path

# Make the harness modules importable as top-level modules (run.py, metrics.py, ...).
HARNESS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HARNESS_DIR))
