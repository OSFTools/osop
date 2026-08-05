# Running the toolkit

The recommended top-level entrypoint is now a Python script controlled by one
YAML file.

## Run command

From the repository root:

```bash
python scripts/run_all.py --config osop_config.yml
```

For a safe command preview that does not download or process data:

```bash
python scripts/run_all.py --config osop_config.yml --dry-run
```

## Configuration file

All user options are in `osop_config.yml`:

- workflow switches (`hindcast`, `forecast`)
- core parameters (`month`, `leads`, `area`, `variable`, `location`)
- pycpt options (`pycpt`, `predictor_area`)
- forecast year (`forecast_year`)
- output directories (`paths`)
- service IDs and weights (`services`)

This replaces editing shell variables directly.

## Viewing outputs

Outputs are written to the path configured in `paths.base`.

- Hindcast plots: `{base}/hindcast/plots`
- Forecast plots: `{base}/forecast/plots`
- Logs: `{base}/logfiles`

If the default config is unchanged, output is under `./output/single_script`.
