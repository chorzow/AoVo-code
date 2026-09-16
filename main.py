"""Runs the analysis pipeline."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = ROOT_DIR / "scripts"

# interpreter=None means "an R script" -> resolved to r_runner in main()
PIPELINE = [
    ("01-preload.R", None),
    ("02-rarecurve.R", None),
    ("03-alpha-diversity.R", None),
    ("04-beta-diversity.R", None),
    ("05-taxonomic-analysis.R", None),
    ("06-venn-diagrams.py", [sys.executable]),
    ("07-indicator-species.R", None),
]

R_RUNNERS = {
    # `uvr run` uses this project's uvr-managed R version/package library.
    "uvr": ["uvr", "run"],
    # Plain Rscript, for environments without uvr (relies on whatever
    # packages are already installed for the R on PATH).
    "rscript": ["Rscript"],
}


def resolve_r_runner():
    """Pick how to invoke R scripts: $R_RUNNER if set, else uvr if this
    project uses it and uvr is installed, else plain Rscript."""
    choice = os.environ.get("R_RUNNER", "").strip().lower()
    if choice:
        if choice not in R_RUNNERS:
            raise SystemExit(
                f"R_RUNNER={choice!r} is not one of {sorted(R_RUNNERS)}"
            )
        return R_RUNNERS[choice]

    if (ROOT_DIR / "uvr.toml").exists() and shutil.which("uvr"):
        return R_RUNNERS["uvr"]

    if shutil.which("Rscript"):
        return R_RUNNERS["rscript"]

    raise SystemExit(
        "No R runner found: install uvr, or make sure Rscript is on PATH."
    )


def main():
    r_runner = resolve_r_runner()
    print(f"Using R runner: {' '.join(r_runner)}")

    for script_name, interpreter in PIPELINE:
        script_path = SCRIPTS_DIR / script_name
        command = (interpreter or r_runner) + [str(script_path)]
        print(f"\n=== Running {script_name} ===")
        subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
