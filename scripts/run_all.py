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
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

import cptcore
import yaml

logger = logging.getLogger(__name__)

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


def _resolve_paths(paths_cfg: dict[str, Any]) -> dict[str, Any]:
    if "base" not in paths_cfg:
        raise ConfigError("paths.base is required")

    base = Path(str(paths_cfg["base"]))
    formatted: dict[str, Any] = {"base": str(base)}

    for key, raw in paths_cfg.items():
        if key == "base":
            continue
        try:
            if isinstance(raw, dict):
                formatted[key] = {
                    inner_key: str(inner_value).format(base=str(base))
                    for inner_key, inner_value in raw.items()
                }
            else:
                formatted[key] = str(raw).format(base=str(base))
        except (KeyError, ValueError) as exc:
            raise ConfigError(f"Invalid placeholder in paths.{key}: {exc}") from exc

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

    try:
        month = int(parameters["month"])
        leads = str(parameters["leads"])
        area = str(parameters["area"])
        variable = str(parameters["variable"])
        pycpt_value = parameters["pycpt"]
    except KeyError as exc:
        raise ConfigError(f"Missing required parameters field: {exc.args[0]}") from exc
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"Invalid parameters value: {exc}") from exc

    if month < 1 or month > 12:
        raise ConfigError(f"Month must be 1-12, got: {month}")

    parse_leads(leads)

    parse_bbox(area, "parameters.area")

    if variable not in ALLOWED_VARIABLES:
        raise ConfigError(
            f"Unsupported variable: {variable}. Allowed: {sorted(ALLOWED_VARIABLES)}"
        )

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
    if not isinstance(centres_cfg, list) or not centres_cfg:
        raise ConfigError("centres must be a non-empty list")

    if not isinstance(config["services"], dict) or not config["services"]:
        raise ConfigError("services must be a non-empty mapping")

    resolved_paths = _resolve_paths(config["paths"])

    required_top = {"logdir", "hindcast", "forecast"}
    missing_top = required_top - set(resolved_paths)
    if missing_top:
        raise ConfigError(f"Missing required paths entries: {sorted(missing_top)}")

    required_nested = {"downloads", "products", "scores", "plots", "pycpt"}
    for section in ("hindcast", "forecast"):
        if not isinstance(resolved_paths.get(section), dict):
            raise ConfigError(f"paths.{section} must be a mapping")
        missing = required_nested - set(resolved_paths[section])
        if missing:
            raise ConfigError(
                f"Missing required paths.{section} entries: {sorted(missing)}"
            )
    return {
        "workflow": {
            "hindcast": bool(workflow.get("hindcast", True)),
            "forecast": forecast_enabled,
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


def _run_step(command: list[str], env: dict[str, str] | None = None) -> int:
    cmd_display = " ".join(command)
    logger.info("[RUN] %s", cmd_display)
    completed = subprocess.run(command, check=False, env=env)
    return completed.returncode


def _build_subprocess_env(script_dir: Path) -> dict[str, str]:
    """Build environment for subprocesses with repo lib path on PYTHONPATH.

    Parameters
    ----------
    script_dir : pathlib.Path
        Directory containing the orchestration scripts.

    Returns
    -------
    dict of str to str
        Environment mapping for subprocess execution.
    """
    repo_root = script_dir.parent
    lib_path = str((repo_root / "lib").resolve())

    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    if existing:
        env["PYTHONPATH"] = os.pathsep.join([lib_path, existing])
    else:
        env["PYTHONPATH"] = lib_path

    conda_prefix = env.get("CONDA_PREFIX")
    if conda_prefix:
        if sys.platform.startswith("win") or os.name == "nt":
            cpt_bin_dir = Path(conda_prefix) / "Library" / "cpt"
            env["CPT_BIN_DIR"] = str(cpt_bin_dir)
        else:
            cpt_bin_dir = Path(conda_prefix) / "bin"

        existing_path = env.get("PATH", "")
        path_entries = [str(cpt_bin_dir)]
        if existing_path:
            path_entries.append(existing_path)
        env["PATH"] = os.pathsep.join(path_entries)

    return env


def configure_pycpt(subprocess_env: dict[str, str]) -> None:
    """Apply Windows-specific pyCPT configuration.

    Parameters
    ----------
    subprocess_env : dict of str to str
        Environment mapping used by workflow subprocesses.
    """
    if not (sys.platform.startswith("win") or os.name == "nt"):
        return

    try:
        cpt_bin_dir = subprocess_env.get("CPT_BIN_DIR")

        if not cpt_bin_dir:
            return

        cpt_exe = Path(cpt_bin_dir) / "CPT.exe"

        if not cpt_exe.exists():
            logger.info(f"WARNING: CPT executable not found: {cpt_exe}")
            return

        orig_init = cptcore.CPT.__init__

        def patched_init(self, *args, **kwargs):
            kwargs.setdefault("cpt_executable", str(cpt_exe))
            return orig_init(self, *args, **kwargs)

        cptcore.CPT.__init__ = patched_init
        logger.info("pyCPT Windows only monkeypatch applied in memory")
    except Exception as exc:
        print("WARNING: Windows pyCPT monkeypatch failed:", exc)


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
    selected = config["centres"]
    if not isinstance(selected, list) or not selected:
        raise ConfigError("Selected centres list is empty")
    return [str(item) for item in selected]


def _bool_str(value: bool) -> str:
    return "True" if value else "False"


def run_pipeline(config: dict[str, Any], script_dir: Path) -> int:
    """Execute configured hindcast and forecast workflows.

    Parameters
    ----------
    config : dict of str to Any
        Validated configuration mapping.
    script_dir : pathlib.Path
        Directory containing the workflow scripts to execute.

    Returns
    -------
    int
        Process-style status code: ``0`` for success, ``1`` if any step fails.
    """
    paths = config["paths"]
    params = config["parameters"]
    centres = _select_centres(config)
    pycpt_arg = _bool_str(params["pycpt"])
    subprocess_env = _build_subprocess_env(script_dir)
    configure_pycpt(subprocess_env)

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
        if _run_step(era5_cmd, env=subprocess_env) != 0:
            failures.append("era5")
        else:
            print("ERA5 download completed successfully")

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
                if _run_step(download_cmd, env=subprocess_env) != 0:
                    failures.append(f"hindcast-download:{centre}")
                    continue
                else:
                    print(f"Hindcast download for {centre} completed successfully")

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
            if _run_step(products_cmd, env=subprocess_env) != 0:
                failures.append(f"hindcast-products:{centre}")
                continue
            else:
                print(
                    f"Hindcast products computation for {centre} completed successfully"
                )

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
            if _run_step(scores_cmd, env=subprocess_env) != 0:
                failures.append(f"hindcast-scores:{centre}")
                continue
            else:
                print(
                    f"Hindcast scores computation for {centre} completed successfully"
                )

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
            if _run_step(plots_cmd, env=subprocess_env) != 0:
                failures.append(f"hindcast-plots:{centre}")
                continue
            else:
                print(f"Hindcast plots generation for {centre} completed successfully")

            if _run_step(plots_cmd, env=subprocess_env) != 0:
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
                if _run_step(forecast_download_cmd, env=subprocess_env) != 0:
                    failures.append(f"forecast-download:{centre}")
                    continue
                else:
                    print(f"Forecast download for {centre} completed successfully")

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
            if _run_step(forecast_products_cmd, env=subprocess_env) != 0:
                failures.append(f"forecast-products:{centre}")
                continue
            else:
                print(
                    f"Forecast products computation for {centre} completed successfully"
                )

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
            if _run_step(forecast_plots_cmd, env=subprocess_env) != 0:
                failures.append(f"forecast-plots:{centre}")
            else:
                print(f"Forecast plots generation for {centre} completed successfully")

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
        validated = validate_config(loaded)

        script_dir = Path(__file__).resolve().parent
        return run_pipeline(validated, script_dir)
    except ConfigError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
