import sys
from pathlib import Path

# make benchmarks/simulate importable as plain modules (trivial, strec, ...)
SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))
