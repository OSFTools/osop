#!/usr/bin/env python3
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.
"""Cross-platform top-level runner for OSOP workflows.

This script replaces shell orchestration with a single YAML configuration file.
It preserves the existing step order and subprocess interfaces used by the
legacy shell runner while avoiding shell-specific behavior.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

import yaml

ALLOWED_VARIABLES = {"2m_temperature", "total_precipitation"}


class ConfigError(Exception):
    """Raised when the orchestrator configuration is invalid."""


def _parse_bool(value: str) -> bool:
    lowered = value.strip().lower()
    if lowered in {"1", "true", "yes", "y", "on"}:
        return True
    if lowered in {"0", "false", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")


def parse_leads(leads: str) -> list[int]:
    """Parse comma-separated lead months into positive integers.

    Parameters
    ----------
    leads : str
        Comma-separated lead month values.

    Returns
    -------
    list of int
        Parsed lead month values as positive integers.
    """
    try:
        parsed = [int(item.strip()) for item in leads.split(",") if item.strip()]
    except ValueError as exc:
        raise ConfigError(f"Invalid leads string: {leads}") from exc

    if not parsed or any(val <= 0 for val in parsed):
        raise ConfigError(f"Leads must be positive integers, got: {leads}")
    return parsed


def parse_bbox(bbox: str, name: str) -> list[float]:
    """Parse a bounding-box string into four numeric values.

    Parameters
    ----------
    bbox : str
        Bounding-box string in ``N,W,S,E`` order.
    name : str
        Configuration field name used in validation error messages.

    Returns
    -------
    list of float
        Bounding-box values as ``[north, west, south, east]``.
    """
    parts = [item.strip() for item in bbox.split(",") if item.strip()]
    if len(parts) != 4:
        raise ConfigError(f"{name} must contain 4 comma-separated numbers, got: {bbox}")
    try:
        return [float(item) for item in parts]
    except ValueError as exc:
        raise ConfigError(f"{name} must contain only numbers, got: {bbox}") from exc


def load_config(config_path: Path) -> dict[str, Any]:
    """Load a YAML configuration mapping from disk.

    Parameters
    ----------
    config_path : pathlib.Path
        Path to the YAML configuration file.

    Returns
    -------
    dict of str to Any
        Parsed top-level configuration mapping.
    """
    if not config_path.exists():
        raise ConfigError(f"Config file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)

    if not isinstance(loaded, dict):
        raise ConfigError("Top-level config must be a YAML mapping")
    return loaded


def apply_cli_overrides(
    config: dict[str, Any], args: argparse.Namespace
) -> dict[str, Any]:
    """Apply command-line overrides to the loaded configuration.

    Parameters
    ----------
    config : dict of str to Any
        Base configuration loaded from YAML.
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    dict of str to Any
        Configuration mapping with CLI overrides applied.
    """
    result = dict(config)
    result["workflow"] = dict(config.get("workflow", {}))
    result["parameters"] = dict(config.get("parameters", {}))

    if args.hindcast_only and args.forecast_only:
        raise ConfigError("Choose only one of --hindcast-only or --forecast-only")

    if args.hindcast_only:
        result["workflow"]["hindcast"] = True
        result["workflow"]["forecast"] = False
    if args.forecast_only:
        result["workflow"]["hindcast"] = False
        result["workflow"]["forecast"] = True
    if args.test_mode:
        result["workflow"]["test_mode"] = True

    if args.month is not None:
        result["parameters"]["month"] = args.month
    if args.leads is not None:
        result["parameters"]["leads"] = args.leads
    if args.area is not None:
        result["parameters"]["area"] = args.area
    if args.variable is not None:
        result["parameters"]["variable"] = args.variable
    if args.location is not None:
        result["parameters"]["location"] = args.location
    if args.method is not None:
        result["parameters"]["method"] = args.method
    if args.pycpt is not None:
        result["parameters"]["pycpt"] = args.pycpt
    if args.predictor_area is not None:
        result["parameters"]["predictor_area"] = args.predictor_area
    if args.forecast_year is not None:
        result["parameters"]["forecast_year"] = args.forecast_year

    return result


def _resolve_paths(paths_cfg: dict[str, Any]) -> dict[str, str]:
    if "base" not in paths_cfg:
        raise ConfigError("paths.base is required")

    base = Path(str(paths_cfg["base"]))
    formatted: dict[str, str] = {"base": str(base)}

    for key, raw in paths_cfg.items():
        if key == "base":
            continue
        if isinstance(raw, dict):
            formatted[key] = {}
            for inner_key, inner_value in raw.items():
                formatted[key][inner_key] = str(inner_value).format(base=str(base))
        else:
            formatted[key] = str(raw).format(base=str(base))

    return formatted


def validate_config(config: dict[str, Any]) -> dict[str, Any]:
    """Validate and normalize orchestrator configuration values.

    Raises ``ConfigError`` if any required fields are missing or invalid.

    Parameters
    ----------
    config : dict of str to Any
        Raw configuration mapping.

    Returns
    -------
    dict of str to Any
        Normalized configuration mapping ready for execution.
    """
    required_sections = ["workflow", "parameters", "paths", "centres", "services"]
    for section in required_sections:
        if section not in config:
            raise ConfigError(f"Missing required section: {section}")

    workflow = config["workflow"]
    parameters = config["parameters"]

    month = int(parameters["month"])
    if month < 1 or month > 12:
        raise ConfigError(f"Month must be 1-12, got: {month}")

    leads = str(parameters["leads"])
    parse_leads(leads)

    area = str(parameters["area"])
    parse_bbox(area, "parameters.area")

    variable = str(parameters["variable"])
    if variable not in ALLOWED_VARIABLES:
        raise ConfigError(
            f"Unsupported variable: {variable}. Allowed: {sorted(ALLOWED_VARIABLES)}"
        )

    pycpt_value = parameters["pycpt"]
    if isinstance(pycpt_value, str):
        pycpt_enabled = _parse_bool(pycpt_value)
    else:
        pycpt_enabled = bool(pycpt_value)

    predictor_area = str(parameters.get("predictor_area", ""))
    if pycpt_enabled:
        parse_bbox(predictor_area, "parameters.predictor_area")

    forecast_enabled = bool(workflow.get("forecast", False))
    if forecast_enabled and "forecast_year" not in parameters:
        raise ConfigError(
            "parameters.forecast_year is required when forecast is enabled"
        )

    centres_cfg = config["centres"]
    if "test" not in centres_cfg or "full" not in centres_cfg:
        raise ConfigError("centres.test and centres.full are required")

    if not isinstance(config["services"], dict) or not config["services"]:
        raise ConfigError("services must be a non-empty mapping")

    resolved_paths = _resolve_paths(config["paths"])

    return {
        "workflow": {
            "hindcast": bool(workflow.get("hindcast", True)),
            "forecast": forecast_enabled,
            "test_mode": bool(workflow.get("test_mode", False)),
            "continue_on_error": bool(workflow.get("continue_on_error", True)),
        },
        "parameters": {
            "month": month,
            "leads": leads,
            "area": area,
            "variable": variable,
            "location": str(parameters.get("location", "None")),
            "method": (
                None
                if parameters.get("method") in (None, "None", "")
                else str(parameters.get("method"))
            ),
            "pycpt": pycpt_enabled,
            "predictor_area": predictor_area,
            "forecast_year": int(parameters.get("forecast_year", 0)),
        },
        "paths": resolved_paths,
        "centres": config["centres"],
        "services": config["services"],
    }


def write_services_yaml(services: dict[str, Any], destination: Path) -> None:
    """Write the ``parseyml.yml`` file consumed by downstream scripts.

    Parameters
    ----------
    services : dict of str to Any
        Service configuration payload.
    destination : pathlib.Path
        Output location for the YAML file.

    Returns
    -------
    None
        This function writes a file as a side effect.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = {"Services": services}
    with destination.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)


def ensure_directories(paths: dict[str, Any]) -> None:
    """Create all configured output directories.

    Parameters
    ----------
    paths : dict of str to Any
        Resolved paths mapping containing workflow output directories.

    Returns
    -------
    None
        This function creates directories as a side effect.
    """
    directory_targets = [
        paths["logdir"],
        paths["hindcast"]["downloads"],
        paths["hindcast"]["products"],
        paths["hindcast"]["scores"],
        paths["hindcast"]["plots"],
        paths["hindcast"]["pycpt"],
        paths["forecast"]["downloads"],
        paths["forecast"]["products"],
        paths["forecast"]["scores"],
        paths["forecast"]["plots"],
        paths["forecast"]["pycpt"],
    ]
    for target in directory_targets:
        Path(target).mkdir(parents=True, exist_ok=True)


def _run_step(command: list[str], dry_run: bool) -> int:
    cmd_display = " ".join(command)
    print(f"[RUN] {cmd_display}")
    if dry_run:
        return 0

    completed = subprocess.run(command, check=False)
    return completed.returncode


def _build_common_args(config: dict[str, Any], downloaddir: str) -> list[str]:
    params = config["parameters"]
    return [
        "--month",
        str(params["month"]),
        "--leads",
        str(params["leads"]),
        "--area",
        str(params["area"]),
        "--variable",
        str(params["variable"]),
        "--downloaddir",
        downloaddir,
        "--logdir",
        config["paths"]["logdir"],
    ]


def _select_centres(config: dict[str, Any]) -> list[str]:
    selected = (
        config["centres"]["test"]
        if config["workflow"]["test_mode"]
        else config["centres"]["full"]
    )
    if not isinstance(selected, list) or not selected:
        raise ConfigError("Selected centres list is empty")
    return [str(item) for item in selected]


def _bool_str(value: bool) -> str:
    return "True" if value else "False"


def run_pipeline(config: dict[str, Any], script_dir: Path, dry_run: bool) -> int:
    """Execute configured hindcast and forecast workflows.

    Parameters
    ----------
    config : dict of str to Any
        Validated configuration mapping.
    script_dir : pathlib.Path
        Directory containing the workflow scripts to execute.
    dry_run : bool
        If ``True``, print commands without running subprocesses.

    Returns
    -------
    int
        Process-style status code: ``0`` for success, ``1`` if any step fails.
    """
    paths = config["paths"]
    params = config["parameters"]
    centres = _select_centres(config)
    pycpt_arg = _bool_str(params["pycpt"])

    ensure_directories(paths)

    parseyml_path = Path(paths["hindcast"]["downloads"]) / "parseyml.yml"
    write_services_yaml(config["services"], parseyml_path)
    fc_parseyml = Path(paths["forecast"]["downloads"]) / "parseyml.yml"
    shutil.copy2(parseyml_path, fc_parseyml)

    failures: list[str] = []

    if config["workflow"]["hindcast"]:
        era5_cmd = [
            sys.executable,
            str(script_dir / "get_era5.py"),
            *_build_common_args(config, paths["hindcast"]["downloads"]),
            "--pycpt",
            pycpt_arg,
            "--pycptdir",
            paths["hindcast"]["pycpt"],
        ]
        if _run_step(era5_cmd, dry_run) != 0:
            failures.append("era5")

        for centre in centres:
            if centre != "mme":
                download_cmd = [
                    sys.executable,
                    str(script_dir / "get_any_hindcast.py"),
                    "--centre",
                    centre,
                    *_build_common_args(config, paths["hindcast"]["downloads"]),
                    "--predictor_area",
                    params["predictor_area"],
                    "--pycpt",
                    pycpt_arg,
                    "--pycptdir",
                    paths["hindcast"]["pycpt"],
                ]
                if _run_step(download_cmd, dry_run) != 0:
                    failures.append(f"hindcast-download:{centre}")
                    continue

            products_cmd = [
                sys.executable,
                str(script_dir / "compute_products.py"),
                "--centre",
                centre,
                *_build_common_args(config, paths["hindcast"]["downloads"]),
                "--productsdir",
                paths["hindcast"]["products"],
                "--predictor_area",
                params["predictor_area"],
                "--pycpt",
                pycpt_arg,
                "--pycptdir",
                paths["hindcast"]["pycpt"],
            ]
            if _run_step(products_cmd, dry_run) != 0:
                failures.append(f"hindcast-products:{centre}")
                continue

            scores_cmd = [
                sys.executable,
                str(script_dir / "compute_scores.py"),
                "--centre",
                centre,
                *_build_common_args(config, paths["hindcast"]["downloads"]),
                "--scoresdir",
                paths["hindcast"]["scores"],
                "--productsdir",
                paths["hindcast"]["products"],
            ]
            if _run_step(scores_cmd, dry_run) != 0:
                failures.append(f"hindcast-scores:{centre}")
                continue

            plots_cmd = [
                sys.executable,
                str(script_dir / "plot_verification.py"),
                "--location",
                params["location"],
                "--centre",
                centre,
                *_build_common_args(config, paths["hindcast"]["downloads"]),
                "--scoresdir",
                paths["hindcast"]["scores"],
                "--plotdir",
                paths["hindcast"]["plots"],
            ]
            if params["method"] is not None:
                plots_cmd.extend(["--method", params["method"]])

            if _run_step(plots_cmd, dry_run) != 0:
                failures.append(f"hindcast-plots:{centre}")

    if config["workflow"]["forecast"]:
        for centre in centres:
            if centre != "mme":
                forecast_download_cmd = [
                    sys.executable,
                    str(script_dir / "get_any_hindcast.py"),
                    "--centre",
                    centre,
                    *_build_common_args(config, paths["forecast"]["downloads"]),
                    "--years",
                    str(params["forecast_year"]),
                    "--predictor_area",
                    params["predictor_area"],
                    "--pycpt",
                    pycpt_arg,
                    "--pycptdir",
                    paths["forecast"]["pycpt"],
                ]
                if _run_step(forecast_download_cmd, dry_run) != 0:
                    failures.append(f"forecast-download:{centre}")
                    continue

            forecast_products_cmd = [
                sys.executable,
                str(script_dir / "forecast_products.py"),
                "--centre",
                centre,
                *_build_common_args(config, paths["forecast"]["downloads"]),
                "--downloadhcdir",
                paths["hindcast"]["downloads"],
                "--productshcdir",
                paths["hindcast"]["products"],
                "--productsfcdir",
                paths["forecast"]["products"],
                "--yearsfc",
                str(params["forecast_year"]),
                "--predictor_area",
                params["predictor_area"],
                "--pycpt",
                pycpt_arg,
                "--pycptdir",
                paths["forecast"]["pycpt"],
                "--hindcast_pycptdir",
                paths["hindcast"]["pycpt"],
            ]
            if _run_step(forecast_products_cmd, dry_run) != 0:
                failures.append(f"forecast-products:{centre}")
                continue

            forecast_plots_cmd = [
                sys.executable,
                str(script_dir / "forecast_plots.py"),
                "--location",
                params["location"],
                "--centre",
                centre,
                *_build_common_args(config, paths["forecast"]["downloads"]),
                "--productsfcdir",
                paths["forecast"]["products"],
                "--plotsdir",
                paths["forecast"]["plots"],
                "--yearsfc",
                str(params["forecast_year"]),
            ]
            if _run_step(forecast_plots_cmd, dry_run) != 0:
                failures.append(f"forecast-plots:{centre}")

    if failures:
        print("Run completed with failures:")
        for failure in failures:
            print(f" - {failure}")
        return 1

    print("Run completed successfully")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser configured for the orchestration CLI.
    """
    parser = argparse.ArgumentParser(description="Run OSOP workflow from YAML config")
    parser.add_argument(
        "--config",
        default=str(Path(__file__).resolve().parents[1] / "osop_config.yml"),
        help="Path to YAML configuration file",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print commands only")
    parser.add_argument("--test-mode", action="store_true", help="Use test centres")
    parser.add_argument("--hindcast-only", action="store_true")
    parser.add_argument("--forecast-only", action="store_true")

    parser.add_argument("--month", type=int)
    parser.add_argument("--leads")
    parser.add_argument("--area")
    parser.add_argument("--variable", choices=sorted(ALLOWED_VARIABLES))
    parser.add_argument("--location")
    parser.add_argument("--method")
    parser.add_argument("--pycpt", type=_parse_bool)
    parser.add_argument("--predictor-area")
    parser.add_argument("--forecast-year", type=int)

    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the orchestration CLI entrypoint.

    Parameters
    ----------
    argv : list of str or None, optional
        Optional argument vector passed to ``argparse``.

    Returns
    -------
    int
        Exit status code returned by validation and pipeline execution.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        config_path = Path(args.config).expanduser().resolve()
        loaded = load_config(config_path)
        overridden = apply_cli_overrides(loaded, args)
        validated = validate_config(overridden)

        script_dir = Path(__file__).resolve().parent
        return run_pipeline(validated, script_dir, dry_run=args.dry_run)
    except ConfigError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
