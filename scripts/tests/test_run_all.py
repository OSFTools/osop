"""Tests for the Python top-level orchestrator."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import yaml


def _load_run_all_module():
    module_path = Path(__file__).resolve().parents[1] / "run_all.py"
    spec = importlib.util.spec_from_file_location("run_all", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _minimal_config(tmp_path: Path) -> dict:
    base = tmp_path / "out"
    return {
        "workflow": {
            "hindcast": False,
            "forecast": False,
            "test_mode": False,
            "continue_on_error": True,
        },
        "parameters": {
            "month": 5,
            "leads": "2,3,4",
            "area": "39,60,-11,141",
            "variable": "total_precipitation",
            "location": "None",
            "method": "pmesh",
            "pycpt": True,
            "predictor_area": "40,0,-40,359",
            "forecast_year": 2025,
        },
        "paths": {
            "base": str(base),
            "logdir": "{base}/logfiles",
            "hindcast": {
                "downloads": "{base}/hindcast/downloads",
                "products": "{base}/hindcast/products",
                "scores": "{base}/hindcast/scores",
                "plots": "{base}/hindcast/plots",
                "pycpt": "{base}/hindcast/pycpt",
            },
            "forecast": {
                "downloads": "{base}/forecast/downloads",
                "products": "{base}/forecast/products",
                "scores": "{base}/forecast/scores",
                "plots": "{base}/forecast/plots",
                "pycpt": "{base}/forecast/pycpt",
            },
        },
        "centres": {
            "test": ["meteo_france", "ukmo", "mme"],
            "full": ["meteo_france", "ukmo", "mme"],
        },
        "services": {
            "meteo_france": [9, 1],
            "ukmo": [604, 1],
            "mme": [1, 0],
            "jma": [3, 0],
        },
    }


def test_validate_config_accepts_minimal(tmp_path: Path):
    run_all = _load_run_all_module()
    config = _minimal_config(tmp_path)

    validated = run_all.validate_config(config)

    assert validated["workflow"]["hindcast"] is False
    assert validated["workflow"]["forecast"] is False
    assert validated["parameters"]["month"] == 5


def test_validate_config_rejects_bad_month(tmp_path: Path):
    run_all = _load_run_all_module()
    config = _minimal_config(tmp_path)
    config["parameters"]["month"] = 13

    try:
        run_all.validate_config(config)
        assert False, "Expected ConfigError"
    except run_all.ConfigError as exc:
        assert "Month must be 1-12" in str(exc)


def test_dry_run_writes_services_and_returns_success(tmp_path: Path):
    run_all = _load_run_all_module()
    config = _minimal_config(tmp_path)
    validated = run_all.validate_config(config)

    rc = run_all.run_pipeline(
        validated,
        script_dir=Path(__file__).resolve().parents[1],
        dry_run=True,
    )

    assert rc == 0
    parseyml_hc = Path(validated["paths"]["hindcast"]["downloads"]) / "parseyml.yml"
    parseyml_fc = Path(validated["paths"]["forecast"]["downloads"]) / "parseyml.yml"
    assert parseyml_hc.exists()
    assert parseyml_fc.exists()

    with parseyml_hc.open("r", encoding="utf-8") as handle:
        parsed = yaml.safe_load(handle)
    assert "Services" in parsed
    assert parsed["Services"]["ukmo"] == [604, 1]
