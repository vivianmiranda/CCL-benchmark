"""Portable archive and scratch paths for the frozen QAG diagnostic."""
import argparse
import os
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
INPUTS = PACKAGE / "inputs"
VENDOR = PACKAGE / "vendor"
RESULTS = PACKAGE / "results"
DEFAULT_OUTPUT = PACKAGE.parents[1] / "work" / "limber_false_success_rerun"


def output_dir():
    """Keep all newly calculated results outside the committed archive."""
    path = Path(os.environ.get("QAG_DIAGNOSTIC_OUTPUT", DEFAULT_OUTPUT)).resolve()
    if path == PACKAGE or PACKAGE in path.parents:
        raise ValueError("QAG_DIAGNOSTIC_OUTPUT must be outside the saved diagnostic")
    path.mkdir(parents=True, exist_ok=True)
    return path


def plot_paths():
    """Plot saved results by default, without rewriting their archive."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT / "figures")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    return args.results.resolve(), args.output.resolve()
