"""Tests for the Python top-level orchestrator."""

from __future__ import annotations

import argparse
import functools
import importlib.util
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import yaml


def _load_run_all_module():
    module_path = Path(__file__).resolve().parents[1] / "run_all.py"
    spec = importlib.util.spec_from_file_location("run_all", module_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@functools.lru_cache(maxsize=1)
def _get_run_all():
    return _load_run_all_module()


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


# ---------------------------------------------------------------------------
# _parse_bool
# ---------------------------------------------------------------------------


def test_parse_bool_truthy_values():
    """All recognised truthy string representations return True."""
    run_all = _get_run_all()
    for val in ("1", "true", "True", "TRUE", "yes", "YES", "y", "Y", "on", "ON"):
        assert run_all._parse_bool(val) is True


def test_parse_bool_falsy_values():
    """All recognised falsy string representations return False."""
    run_all = _get_run_all()
    for val in ("0", "false", "False", "no", "NO", "n", "N", "off", "OFF"):
        assert run_all._parse_bool(val) is False


def test_parse_bool_invalid_raises():
    """An unrecognised value raises ArgumentTypeError."""
    run_all = _get_run_all()
    with pytest.raises(argparse.ArgumentTypeError):
        run_all._parse_bool("maybe")


# ---------------------------------------------------------------------------
# parse_leads
# ---------------------------------------------------------------------------


def test_parse_leads_valid():
    """Comma-separated positive integers are parsed into a list."""
    run_all = _get_run_all()
    assert run_all.parse_leads("2,3,4") == [2, 3, 4]
    assert run_all.parse_leads("1") == [1]


def test_parse_leads_invalid_string_raises():
    """Non-integer tokens raise ConfigError."""
    run_all = _get_run_all()
    with pytest.raises(run_all.ConfigError, match="Invalid leads string"):
        run_all.parse_leads("a,b")


def test_parse_leads_zero_raises():
    """A zero value raises ConfigError."""
    run_all = _get_run_all()
    with pytest.raises(run_all.ConfigError, match="Leads must be positive integers"):
        run_all.parse_leads("0,1")


def test_parse_leads_negative_raises():
    """A negative value raises ConfigError."""
    run_all = _get_run_all()
    with pytest.raises(run_all.ConfigError, match="Leads must be positive integers"):
        run_all.parse_leads("-1,2")


def test_parse_leads_empty_raises():
    """Whitespace-only input raises ConfigError."""
    run_all = _get_run_all()
    with pytest.raises(run_all.ConfigError, match="Leads must be positive integers"):
        run_all.parse_leads("  ")


# ---------------------------------------------------------------------------
# parse_bbox
# ---------------------------------------------------------------------------


def test_parse_bbox_valid():
    """A valid N,W,S,E string is parsed into four floats."""
    run_all = _get_run_all()
    assert run_all.parse_bbox("39,60,-11,141", "area") == [39.0, 60.0, -11.0, 141.0]


def test_parse_bbox_wrong_count_raises():
    """Fewer than 4 values raises ConfigError."""
    run_all = _get_run_all()
    with pytest.raises(run_all.ConfigError, match="must contain 4 comma-separated"):
        run_all.parse_bbox("1,2,3", "area")


def test_parse_bbox_non_numeric_raises():
    """Non-numeric tokens raise ConfigError."""
    run_all = _get_run_all()
    with pytest.raises(run_all.ConfigError, match="must contain only numbers"):
        run_all.parse_bbox("N,S,E,W", "area")


# ---------------------------------------------------------------------------
# load_config
# ---------------------------------------------------------------------------


def test_load_config_loads_valid_file(tmp_path):
    """A well-formed YAML file is loaded as a dict."""
    run_all = _get_run_all()
    cfg_file = tmp_path / "config.yml"
    cfg_file.write_text("key: value\n", encoding="utf-8")
    assert run_all.load_config(cfg_file) == {"key": "value"}


def test_load_config_missing_file_raises(tmp_path):
    """A missing path raises ConfigError."""
    run_all = _get_run_all()
    with pytest.raises(run_all.ConfigError, match="Config file not found"):
        run_all.load_config(tmp_path / "nonexistent.yml")


def test_load_config_non_mapping_raises(tmp_path):
    """A YAML list at the top level raises ConfigError."""
    run_all = _get_run_all()
    cfg_file = tmp_path / "config.yml"
    cfg_file.write_text("- item1\n- item2\n", encoding="utf-8")
    with pytest.raises(
        run_all.ConfigError, match="Top-level config must be a YAML mapping"
    ):
        run_all.load_config(cfg_file)


# ---------------------------------------------------------------------------
# validate_config
# ---------------------------------------------------------------------------


def test_validate_config_accepts_minimal(tmp_path):
    """A minimal valid config passes validation without errors."""
    run_all = _get_run_all()
    validated = run_all.validate_config(_minimal_config(tmp_path))
    assert validated["workflow"]["hindcast"] is False
    assert validated["workflow"]["forecast"] is False
    assert validated["parameters"]["month"] == 5


def test_validate_config_rejects_bad_month(tmp_path):
    """A month value > 12 raises ConfigError."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["parameters"]["month"] = 13
    with pytest.raises(run_all.ConfigError, match="Month must be 1-12"):
        run_all.validate_config(config)


def test_validate_config_rejects_month_zero(tmp_path):
    """Month 0 raises ConfigError."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["parameters"]["month"] = 0
    with pytest.raises(run_all.ConfigError, match="Month must be 1-12"):
        run_all.validate_config(config)


def test_validate_config_rejects_missing_section(tmp_path):
    """A missing top-level section raises ConfigError."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    del config["services"]
    with pytest.raises(run_all.ConfigError, match="Missing required section: services"):
        run_all.validate_config(config)


def test_validate_config_rejects_unsupported_variable(tmp_path):
    """An unsupported variable name raises ConfigError."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["parameters"]["variable"] = "wind_speed"
    with pytest.raises(run_all.ConfigError, match="Unsupported variable"):
        run_all.validate_config(config)


def test_validate_config_pycpt_as_string(tmp_path):
    """A truthy string value for pycpt is coerced to True."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["parameters"]["pycpt"] = "true"
    assert run_all.validate_config(config)["parameters"]["pycpt"] is True


def test_validate_config_pycpt_false_skips_predictor_area(tmp_path):
    """When pycpt is False, predictor_area is not validated."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["parameters"]["pycpt"] = False
    config["parameters"]["predictor_area"] = ""
    assert run_all.validate_config(config)["parameters"]["pycpt"] is False


def test_validate_config_forecast_requires_year(tmp_path):
    """Enabling forecast without forecast_year raises ConfigError."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["forecast"] = True
    del config["parameters"]["forecast_year"]
    with pytest.raises(run_all.ConfigError, match="forecast_year is required"):
        run_all.validate_config(config)


def test_validate_config_rejects_missing_centres_keys(tmp_path):
    """Ensure centres dict missing 'test' or 'full' raises ConfigError."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    del config["centres"]["test"]
    with pytest.raises(
        run_all.ConfigError, match="centres.test and centres.full are required"
    ):
        run_all.validate_config(config)


def test_validate_config_rejects_empty_services(tmp_path):
    """An empty services dict raises ConfigError."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["services"] = {}
    with pytest.raises(
        run_all.ConfigError, match="services must be a non-empty mapping"
    ):
        run_all.validate_config(config)


def test_validate_config_rejects_non_dict_services(tmp_path):
    """A non-dict services value raises ConfigError."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["services"] = ["list"]
    with pytest.raises(
        run_all.ConfigError, match="services must be a non-empty mapping"
    ):
        run_all.validate_config(config)


def test_validate_config_method_none_string(tmp_path):
    """The string 'None' for method is normalised to Python None."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["parameters"]["method"] = "None"
    assert run_all.validate_config(config)["parameters"]["method"] is None


def test_validate_config_method_empty_string(tmp_path):
    """An empty string for method is normalised to Python None."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["parameters"]["method"] = ""
    assert run_all.validate_config(config)["parameters"]["method"] is None


def test_validate_config_paths_base_required(tmp_path):
    """Missing paths.base raises ConfigError."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    del config["paths"]["base"]
    with pytest.raises(run_all.ConfigError, match="paths.base is required"):
        run_all.validate_config(config)


# ---------------------------------------------------------------------------
# write_services_yaml
# ---------------------------------------------------------------------------


def test_write_services_yaml(tmp_path):
    """Ensure services are written to YAML under a 'Services' key, creating parent dirs."""
    run_all = _get_run_all()
    dest = tmp_path / "sub" / "parseyml.yml"
    services = {"ukmo": [604, 1], "mme": [1, 0]}
    run_all.write_services_yaml(services, dest)
    assert dest.exists()
    with dest.open("r", encoding="utf-8") as f:
        assert yaml.safe_load(f) == {"Services": services}


# ---------------------------------------------------------------------------
# ensure_directories
# ---------------------------------------------------------------------------


def test_ensure_directories_creates_all_dirs(tmp_path):
    """All configured output directories are created."""
    run_all = _get_run_all()
    base = str(tmp_path / "out")
    paths = {
        "logdir": f"{base}/logs",
        "hindcast": {
            "downloads": f"{base}/hc/dl",
            "products": f"{base}/hc/prod",
            "scores": f"{base}/hc/scores",
            "plots": f"{base}/hc/plots",
            "pycpt": f"{base}/hc/pycpt",
        },
        "forecast": {
            "downloads": f"{base}/fc/dl",
            "products": f"{base}/fc/prod",
            "scores": f"{base}/fc/scores",
            "plots": f"{base}/fc/plots",
            "pycpt": f"{base}/fc/pycpt",
        },
    }
    run_all.ensure_directories(paths)
    for d in [
        paths["logdir"],
        paths["hindcast"]["downloads"],
        paths["forecast"]["plots"],
    ]:
        assert Path(d).is_dir()


# ---------------------------------------------------------------------------
# _run_step
# ---------------------------------------------------------------------------


def test_run_step_dry_run_returns_zero():
    """In dry-run mode no subprocess is spawned and 0 is returned."""
    run_all = _get_run_all()
    assert run_all._run_step(["echo", "hello"], dry_run=True) == 0


def test_run_step_real_run_success():
    """A subprocess returning 0 causes _run_step to return 0."""
    run_all = _get_run_all()
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        assert run_all._run_step(["echo", "hello"], dry_run=False) == 0


def test_run_step_real_run_failure():
    """A subprocess returning 1 causes _run_step to return 1."""
    run_all = _get_run_all()
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=1)
        assert run_all._run_step(["false"], dry_run=False) == 1


# ---------------------------------------------------------------------------
# _build_subprocess_env
# ---------------------------------------------------------------------------


def test_build_subprocess_env_windows_uses_library_cpt(tmp_path):
    """Prepend CONDA_PREFIX/Library/cpt to PATH on Windows platforms."""
    run_all = _get_run_all()
    script_dir = tmp_path / "scripts"
    script_dir.mkdir()
    conda_prefix = str(tmp_path / "conda")

    with patch.dict(
        run_all.os.environ,
        {"CONDA_PREFIX": conda_prefix, "PATH": "orig-path", "PYTHONPATH": "base-py"},
        clear=True,
    ):
        with patch.object(run_all.sys, "platform", "win32"):
            env = run_all._build_subprocess_env(script_dir)

    assert env["PATH"] == os.pathsep.join(
        [str(Path(conda_prefix) / "Library" / "cpt"), "orig-path"]
    )
    assert env["PYTHONPATH"].startswith(str((script_dir.parent / "lib").resolve()))


def test_build_subprocess_env_non_windows_uses_bin(tmp_path):
    """Prepend CONDA_PREFIX/bin to PATH on non-Windows platforms."""
    run_all = _get_run_all()
    script_dir = tmp_path / "scripts"
    script_dir.mkdir()
    conda_prefix = str(tmp_path / "conda")

    with patch.dict(
        run_all.os.environ,
        {"CONDA_PREFIX": conda_prefix, "PATH": "orig-path", "PYTHONPATH": "base-py"},
        clear=True,
    ):
        with patch.object(run_all.sys, "platform", "linux"):
            env = run_all._build_subprocess_env(script_dir)

    assert env["PATH"] == os.pathsep.join(
        [str(Path(conda_prefix) / "bin"), "orig-path"]
    )
    assert env["PYTHONPATH"].startswith(str((script_dir.parent / "lib").resolve()))


def test_build_subprocess_env_without_conda_prefix_keeps_path(tmp_path):
    """If CONDA_PREFIX is missing, PATH is unchanged while PYTHONPATH is set."""
    run_all = _get_run_all()
    script_dir = tmp_path / "scripts"
    script_dir.mkdir()

    with patch.dict(
        run_all.os.environ,
        {"PATH": "orig-path", "PYTHONPATH": "base-py"},
        clear=True,
    ):
        env = run_all._build_subprocess_env(script_dir)

    assert env["PATH"] == "orig-path"
    assert env["PYTHONPATH"] == os.pathsep.join(
        [str((script_dir.parent / "lib").resolve()), "base-py"]
    )


# ---------------------------------------------------------------------------
# _select_centres
# ---------------------------------------------------------------------------


def test_select_centres_test_mode(tmp_path):
    """When test_mode is True, the 'test' centres list is used."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["test_mode"] = True
    validated = run_all.validate_config(config)
    assert run_all._select_centres(validated) == ["meteo_france", "ukmo", "mme"]


def test_select_centres_full_mode(tmp_path):
    """When test_mode is False, the 'full' centres list is used."""
    run_all = _get_run_all()
    validated = run_all.validate_config(_minimal_config(tmp_path))
    assert run_all._select_centres(validated) == ["meteo_france", "ukmo", "mme"]


def test_select_centres_empty_raises(tmp_path):
    """An empty centres list raises ConfigError."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["centres"]["full"] = []
    with pytest.raises(run_all.ConfigError, match="Selected centres list is empty"):
        run_all._select_centres(run_all.validate_config(config))


# ---------------------------------------------------------------------------
# run_pipeline – hindcast
# ---------------------------------------------------------------------------


def test_hindcast_dry_run(tmp_path):
    """Hindcast pipeline in dry-run mode returns 0 without spawning subprocesses."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = True
    config["workflow"]["forecast"] = False
    rc = run_all.run_pipeline(
        run_all.validate_config(config),
        script_dir=Path(__file__).resolve().parents[1],
        dry_run=True,
    )
    assert rc == 0


def test_hindcast_dry_run_test_mode(tmp_path):
    """Hindcast pipeline with test_mode uses the test centres list."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = True
    config["workflow"]["test_mode"] = True
    rc = run_all.run_pipeline(
        run_all.validate_config(config),
        script_dir=Path(__file__).resolve().parents[1],
        dry_run=True,
    )
    assert rc == 0


def test_hindcast_no_method_dry_run(tmp_path):
    """Hindcast pipeline without a plot method omits --method from the plots command."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = True
    config["parameters"]["method"] = None
    rc = run_all.run_pipeline(
        run_all.validate_config(config),
        script_dir=Path(__file__).resolve().parents[1],
        dry_run=True,
    )
    assert rc == 0


def test_hindcast_step_failure_returns_nonzero(tmp_path):
    """Any failing step causes the pipeline to return 1."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = True
    config["centres"]["full"] = ["ukmo"]
    validated = run_all.validate_config(config)
    with patch.object(run_all, "_run_step", return_value=1):
        assert (
            run_all.run_pipeline(
                validated, Path(__file__).resolve().parents[1], dry_run=False
            )
            == 1
        )


def test_hindcast_mme_skips_download(tmp_path):
    """The mme centre skips the download step but still runs products, scores, and plots."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = True
    config["centres"]["full"] = ["mme"]
    validated = run_all.validate_config(config)

    calls = []
    with patch.object(
        run_all,
        "_run_step",
        side_effect=lambda cmd, dry_run, env=None: calls.append(cmd) or 0,
    ):
        run_all.run_pipeline(
            validated, Path(__file__).resolve().parents[1], dry_run=False
        )

    script_names = [Path(cmd[1]).name for cmd in calls if len(cmd) > 1]
    assert "get_any_hindcast.py" not in script_names
    assert "compute_products.py" in script_names


def test_hindcast_products_failure_skips_scores_plots(tmp_path):
    """A products failure skips scores and plots for that centre."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = True
    config["centres"]["full"] = ["ukmo"]
    validated = run_all.validate_config(config)

    def selective_fail(command, dry_run, env=None):
        return 1 if Path(command[1]).name == "compute_products.py" else 0

    with patch.object(run_all, "_run_step", side_effect=selective_fail):
        assert (
            run_all.run_pipeline(
                validated, Path(__file__).resolve().parents[1], dry_run=False
            )
            == 1
        )


def test_hindcast_scores_failure_skips_plots(tmp_path):
    """A scores failure skips plots for that centre."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = True
    config["centres"]["full"] = ["ukmo"]
    validated = run_all.validate_config(config)

    def selective_fail(command, dry_run, env=None):
        return 1 if Path(command[1]).name == "compute_scores.py" else 0

    with patch.object(run_all, "_run_step", side_effect=selective_fail):
        assert (
            run_all.run_pipeline(
                validated, Path(__file__).resolve().parents[1], dry_run=False
            )
            == 1
        )


def test_hindcast_plots_failure_recorded(tmp_path):
    """A plots failure is recorded in the failures list."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = True
    config["centres"]["full"] = ["ukmo"]
    validated = run_all.validate_config(config)

    def selective_fail(command, dry_run, env=None):
        return 1 if Path(command[1]).name == "plot_verification.py" else 0

    with patch.object(run_all, "_run_step", side_effect=selective_fail):
        assert (
            run_all.run_pipeline(
                validated, Path(__file__).resolve().parents[1], dry_run=False
            )
            == 1
        )


# ---------------------------------------------------------------------------
# run_pipeline – forecast
# ---------------------------------------------------------------------------


def test_forecast_dry_run(tmp_path):
    """Forecast pipeline in dry-run mode returns 0 without spawning subprocesses."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = False
    config["workflow"]["forecast"] = True
    rc = run_all.run_pipeline(
        run_all.validate_config(config),
        script_dir=Path(__file__).resolve().parents[1],
        dry_run=True,
    )
    assert rc == 0


def test_forecast_step_failure_returns_nonzero(tmp_path):
    """Any failing step causes the forecast pipeline to return 1."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = False
    config["workflow"]["forecast"] = True
    config["centres"]["full"] = ["ukmo"]
    validated = run_all.validate_config(config)
    with patch.object(run_all, "_run_step", return_value=1):
        assert (
            run_all.run_pipeline(
                validated, Path(__file__).resolve().parents[1], dry_run=False
            )
            == 1
        )


def test_forecast_mme_skips_download(tmp_path):
    """The mme centre skips forecast download but still runs products and plots."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = False
    config["workflow"]["forecast"] = True
    config["centres"]["full"] = ["mme"]
    validated = run_all.validate_config(config)

    calls = []
    with patch.object(
        run_all,
        "_run_step",
        side_effect=lambda cmd, dry_run, env=None: calls.append(cmd) or 0,
    ):
        run_all.run_pipeline(
            validated, Path(__file__).resolve().parents[1], dry_run=False
        )

    script_names = [Path(cmd[1]).name for cmd in calls if len(cmd) > 1]
    assert "get_any_hindcast.py" not in script_names
    assert "forecast_products.py" in script_names
    assert "forecast_plots.py" in script_names


def test_forecast_products_failure_skips_plots(tmp_path):
    """A forecast products failure skips plots for that centre."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = False
    config["workflow"]["forecast"] = True
    config["centres"]["full"] = ["ukmo"]
    validated = run_all.validate_config(config)

    def selective_fail(command, dry_run, env=None):
        return 1 if Path(command[1]).name == "forecast_products.py" else 0

    with patch.object(run_all, "_run_step", side_effect=selective_fail):
        assert (
            run_all.run_pipeline(
                validated, Path(__file__).resolve().parents[1], dry_run=False
            )
            == 1
        )


def test_forecast_plots_failure_recorded(tmp_path):
    """A forecast plots failure is recorded in the failures list."""
    run_all = _get_run_all()
    config = _minimal_config(tmp_path)
    config["workflow"]["hindcast"] = False
    config["workflow"]["forecast"] = True
    config["centres"]["full"] = ["ukmo"]
    validated = run_all.validate_config(config)

    def selective_fail(command, dry_run, env=None):
        return 1 if Path(command[1]).name == "forecast_plots.py" else 0

    with patch.object(run_all, "_run_step", side_effect=selective_fail):
        assert (
            run_all.run_pipeline(
                validated, Path(__file__).resolve().parents[1], dry_run=False
            )
            == 1
        )


# ---------------------------------------------------------------------------
# run_pipeline – writes parseyml files
# ---------------------------------------------------------------------------


def test_dry_run_writes_services_and_returns_success(tmp_path):
    """Ensure run_pipeline writes parseyml.yml into both hindcast and forecast download dirs."""
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


# ---------------------------------------------------------------------------
# build_parser
# ---------------------------------------------------------------------------


def test_build_parser_defaults():
    """Parser defaults expose only YAML path and dry-run switch."""
    run_all = _get_run_all()
    args = run_all.build_parser().parse_args([])
    assert args.dry_run is False
    assert args.config


def test_build_parser_with_dry_run_flag():
    """--dry-run is parsed correctly."""
    run_all = _get_run_all()
    args = run_all.build_parser().parse_args(["--dry-run"])
    assert args.dry_run is True


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def test_main_dry_run(tmp_path):
    """Ensure main() loads config, validates it, and runs the pipeline in dry-run mode."""
    run_all = _get_run_all()
    cfg_path = tmp_path / "config.yml"
    with cfg_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(_minimal_config(tmp_path), f)
    assert run_all.main(["--config", str(cfg_path), "--dry-run"]) == 0


def test_main_config_error(tmp_path):
    """An invalid config causes main() to print an error and return 2."""
    run_all = _get_run_all()
    cfg_path = tmp_path / "bad.yml"
    cfg_path.write_text("workflow:\n  hindcast: true\n", encoding="utf-8")
    assert run_all.main(["--config", str(cfg_path)]) == 2


def test_main_missing_config(tmp_path):
    """A missing config file causes main() to return 2."""
    run_all = _get_run_all()
    assert run_all.main(["--config", str(tmp_path / "nonexistent.yml")]) == 2
