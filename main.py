"""Source-checkout launcher for the shared command-line interface."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from green_vs_docking.cli.main import main

if __name__ == "__main__":
    main()
