#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Script to calculate download hindcasts, calculate terciles and plot verification measures.

SLURM job submission parameters (use with sbatch):
    --qos=normal
    --mem=128G
    --ntasks=8
    --output=master.txt
    --error=master.err
    --time=00-06:00:00
    --export=NONE
    --mail-type=ALL
"""

import os
from pathlib import Path
import subprocess
import sys

import yaml


def setup_directories(scratch_dir):
    """Create necessary directories for the workflow."""
    base_dir = Path(scratch_dir) / "seafoam/data/master/hindcast"

    dirs = {
        "downloaddir": base_dir / "downloads",
        "productsdir": base_dir / "products",
        "scoresdir": base_dir / "scores",
        "plotdir": base_dir / "plots",
        "logdir": base_dir / "logfiles",
    }

    for dir_path in dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    return dirs


def setup_pythonpath():
    """Set PYTHONPATH relative to script location."""
    script_dir = Path(__file__).resolve().parent
    lib_path = (script_dir / "../lib").resolve()

    pythonpath = os.environ.get("PYTHONPATH", "")
    if pythonpath:
        os.environ["PYTHONPATH"] = f"{pythonpath}:{lib_path}"
    else:
        os.environ["PYTHONPATH"] = str(lib_path)


def create_services_yaml(parseyml_path):
    """Create YML file to pass dictionary parameters."""
    services = {
        "Services": {
            "ecmwf": 51,
            "meteo_france": 9,
            "dwd": 22,
            "cmcc": 35,
            "ncep": 2,
            "jma": 3,
            "eccc_can": 4,
            "eccc_gem5": 5,
            "ukmo": 604,
            "mme": 1,
        }
    }

    with open(parseyml_path, "w") as f:
        yaml.dump(services, f, default_flow_style=False)

    print(f"YML file created: {parseyml_path}")


def run_python_script(script_name, args, description):
    """Run a Python script with given arguments.

    Args:
        script_name: Name of the Python script to run
        args: List of command-line arguments
        description: Description of what the script does

    Returns
    -------
        bool: True if successful, False otherwise
    """
    cmd = [sys.executable, script_name] + args

    try:
        result = subprocess.run(cmd, check=True, capture_output=False)
        print(f"{description} successful")
        return True
    except subprocess.CalledProcessError as e:
        print(f"{description} failed")
        return False


def get_era5_data(month, leads, area, downloaddir, logdir, variable):
    """Download ERA5 data."""
    args = [
        "--month",
        str(month),
        "--leads",
        leads,
        "--area",
        area,
        "--downloaddir",
        str(downloaddir),
        "--logdir",
        str(logdir),
        "--variable",
        variable,
    ]

    return run_python_script("get_era5.py", args, "era5 download")


def process_centre(
    centre,
    month,
    leads,
    area,
    variable,
    downloaddir,
    productsdir,
    scoresdir,
    plotdir,
    logdir,
    location,
    method,
):
    """Process a single centre through all stages of the workflow."""
    # Download hindcast data (skip for mme)
    if centre != "mme":
        args = [
            "--centre",
            centre,
            "--month",
            str(month),
            "--leads",
            leads,
            "--area",
            area,
            "--variable",
            variable,
            "--downloaddir",
            str(downloaddir),
            "--logdir",
            str(logdir),
        ]
        success = run_python_script("get_any_hindcast.py", args, f"{centre} : download")
        if not success:
            print(f"{centre} : download failed")

    # Compute terciles and anomalies
    args = [
        "--centre",
        centre,
        "--month",
        str(month),
        "--leads",
        leads,
        "--area",
        area,
        "--variable",
        variable,
        "--downloaddir",
        str(downloaddir),
        "--productsdir",
        str(productsdir),
        "--logdir",
        str(logdir),
    ]
    success = run_python_script(
        "compute_products.py", args, f"{centre} : products generation"
    )

    # Calculate verification scores
    args = [
        "--centre",
        centre,
        "--month",
        str(month),
        "--leads",
        leads,
        "--area",
        area,
        "--downloaddir",
        str(downloaddir),
        "--scoresdir",
        str(scoresdir),
        "--productsdir",
        str(productsdir),
        "--variable",
        variable,
        "--logdir",
        str(logdir),
    ]
    success = run_python_script(
        "compute_scores.py", args, f"{centre} : scores generation"
    )

    # Plot scores
    args = [
        "--location",
        location,
        "--centre",
        centre,
        "--month",
        str(month),
        "--leads",
        leads,
        "--area",
        area,
        "--downloaddir",
        str(downloaddir),
        "--scoresdir",
        str(scoresdir),
        "--plotdir",
        str(plotdir),
        "--variable",
        variable,
        "--method",
        method,
        "--logdir",
        str(logdir),
    ]
    success = run_python_script(
        "plot_verification.py", args, f"{centre} : plots generation"
    )


def main():
    """Run the workflow."""
    # Get SCRATCH directory from environment
    scratch_dir = os.environ.get("SCRATCH")
    if not scratch_dir:
        print("Error: SCRATCH environment variable not set")
        sys.exit(1)

    # Setup directories
    dirs = setup_directories(scratch_dir)

    # Setup PYTHONPATH
    setup_pythonpath()

    # Create YAML configuration file
    parseyml = dirs["downloaddir"] / "parseyml.yml"
    create_services_yaml(parseyml)

    # Set parameters
    month = 5  # initialisation month
    leads = "2,3,4"  # e.g. if month=5 and leads="2,3,4", valid months are JJA (6,7,8)
    area = "45,-30,-2.5,60"  # sub-area in degrees (N,W,S,E)
    variable = "2m_temperature"  # typically "2m_temperature" or "total_precipitation"
    location = "Morocco"  # Options: 'None', 'UK', 'Morocco', 'SAU' (Saudi Arabia)
    method = "pmesh"  # Remove for smooth plotting on correlation plots

    # Get ERA5 data
    get_era5_data(month, leads, area, dirs["downloaddir"], dirs["logdir"], variable)

    # List of centres to process
    centres = [
        "meteo_france",
        "dwd",
        "cmcc",
        "ncep",
        "ukmo",
        "ecmwf",
        "jma",
        "eccc",
        "mme",
    ]

    # Process each centre
    for centre in centres:
        print(f"\nProcessing {centre}...")
        process_centre(
            centre=centre,
            month=month,
            leads=leads,
            area=area,
            variable=variable,
            downloaddir=dirs["downloaddir"],
            productsdir=dirs["productsdir"],
            scoresdir=dirs["scoresdir"],
            plotdir=dirs["plotdir"],
            logdir=dirs["logdir"],
            location=location,
            method=method,
        )

    print("\nWorkflow complete!")


if __name__ == "__main__":
    main()
