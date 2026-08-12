# Running the toolkit

Inside the Gitbash or Linux terminal, navigating to the osop-main file
and typing `ls` should bring up a list of directories contained within
the toolkit. 

## Run command

The top-level programme is a Python script controlled by a YAML file.

From the repository root run:

```bash
python scripts/run_all.py --config osop_config.yml
```

This runs the python code using the options selected in the
`osop_config.yml` file. If the `---config` option is not provided the
code defaults to osop_config.yml. You can create additional config
files as required.

Summary of run options:

- `--config`: Path to the YAML configuration file (example: `osop_config.yml`).
- The script uses the `workflow`, `parameters`, `paths`, `centres`, and `services` sections from the config.

## Configuration file

The script as downloaded is set up with default options for the
variable, season and area to provide a forecast for. However, the user
will want to set the time frame, location and variable of their
choice.

To do this we will need to edit the `osop_config.yml` file. The first
step here is to open this file which can be found at the top level of
the osop toolkit. It does not matter how this is opened. One option
for people who have followed this text up to this point is to open it
with nano from a terminal as was done earlier.

This file is written in [YAML](https://en.wikipedia.org/wiki/YAML) which
is commonly used for configuration files. It aims to be readable for
human beings while also being capable of being loaded by a computer.

Below we explain each of the contents of the configuration file.
Make the changes you need to run your experiment then run:

```bash
python scripts/run_all.py --config osop_config.yml
```

Keep an eye on the output. Once complete the plots should be available
in the dedicated directory. This will give you the hindcast
verification for your set up and the forecast if requested. 

The configuration file has a number of sections. Example file:

```yaml
workflow:
  hindcast: true
  forecast: false

parameters:
  month: 5
  leads: "2,3,4"
  area: "39,60,-11,141"
  variable: "total_precipitation"
  location: "None"
  method: "pmesh"
  pycpt: true
  predictor_area: "40,0,-40,359"
  forecast_year: 2025

paths:
  base: "./output/single_script"
  logdir: "{base}/logfiles"
  hindcast:
    downloads: "{base}/hindcast/downloads"
    products: "{base}/hindcast/products"
    scores: "{base}/hindcast/scores"
    plots: "{base}/hindcast/plots"
    pycpt: "{base}/hindcast/pycpt"
  forecast:
    downloads: "{base}/forecast/downloads"
    products: "{base}/forecast/products"
    scores: "{base}/forecast/scores"
    plots: "{base}/forecast/plots"
    pycpt: "{base}/forecast/pycpt"

centres:
  - meteo_france
  - dwd
  - cmcc
  - ncep
  - ukmo
  - ecmwf
  - jma
  - eccc
  - bom
  - mme

services:
  ecmwf: [51, 1]
  meteo_france: [9, 1]
  dwd: [22, 1]
  cmcc: [35, 1]
  ncep: [2, 1]
  jma: [3, 0]
  eccc_can: [4, 1]
  eccc_gem5: [5, 1]
  ukmo: [604, 1]
  bom: [2, 1]
  mme: [1, 0]

```

### Workflow

- **`hindcast`**: Run hindcast verification (true/false).
- **`forecast`**: Run generation of forecast outputs (true/false).

### Parameters

- **`month`**: Integer. Month used for initialisation time of the forecast (example: `5`).
- **`leads`**: Comma-separated lead times as a string (example: `"2,3,4"`). With an initialisation 
month of 5 leads 2,3,4 gives a forecast/hindcast averaged across the months of June, July and August. Note that the code follows the C3S conventions and so 0 is not a valid lead time. To use the data from the same month as the forecast initialisation time use a lead of 1.
- **`area`**: Bounding box for hindcast/forecast area as a string `latN,latW,lonS,lonE` (example: `"39,60,-11,141"`).
- **`variable`**: Forecast variable. Valid values are  `total_precipitation` and `2m_temperature`.  
- **`location`**: Optional named location (string) or `None`.
- **`method`**: Plotting method. `pmesh` gives colormesh plots, any other value gives contour plots.
- **`pycpt`**: Boolean toggle to enable pycpt processing.
- **`predictor_area`**: Predictor bbox string used by pycpt (example:
  `"40,0,-40,359"`).  string `latN,latW,lonS,lonE`
- **`forecast_year`**: Year used for forecast (integer).

### Paths

- **`base`**: Base output directory (example: `./output/single_script`). In the following entries `{base}` is replaced with the value of this setting.
- **`logdir`**: Log directory
- **`hindcast` / `forecast`**: Subfolders for `downloads`, `products`, `scores`, `plots`, and `pycpt`.

### Centres

- Ordered list of model centres to include (examples: `meteo_france`, `dwd`, `ecmwf`, `ukmo`, etc.).

### Services

- Mapping of centre names to a two-element array: `[service_id,
  enabled_flag_or_weight]` (example: `ecmwf: [51, 1]`).
  - Service ID is the a unique value which is assigned by CDS trying
     to map as close as possible the version numbering used by the
     forecasting centres.
  - Weight is the relvative weight given to each model when producing
    a multimodel ensemble (MME).


  For a given centre you will need to check two things. A) The service
  you would like to use is included and B) that the system version is
  correct for the time frame and use case. Generally, you want to be
  on the latest version. However, some systems do not run historic
  years until the month is reached for the new model. As such checking
  the data is available on Copernicus Climate Data Store is advised.
  See: [Description of the C3S seasonal
  multi-system](https://confluence.ecmwf.int/spaces/CKB/pages/77213502/Description+of+the+C3S+seasonal+multi-system)
  and/or [Seasonal forecast monthly statistics on single
  levels](https://cds.climate.copernicus.eu/datasets/seasonal-monthly-single-levels?tab=overview)

## Viewing outputs

After running the script you will probably want to review the output of your
work. Outputs are written to the path configured in the yml file as
explained above.

- Hindcast plots: `{base}/hindcast/plots`
- Forecast plots: `{base}/forecast/plots`
- Logs: `{base}/logfiles`

If the default config is unchanged, output is under `./output/single_script`.

The rest of these instructions assume the files are stored under the
output folder. These instructions show how to navigate via terminal
using Gitbash. Type `cd output/data/master/hindcast/plots` and then
type `ls`. This will bring up a list of appropriate plots. From here
copy a file name and type code `pastefile` and enter. This should
bring up an image output.

If you receive an error and want to debug it, you may wish to see the
logfiles. To do that you can go to `cd
./output/data/master/hindcast/logfiles`. This is where output of the
code is stored for checking each time it is run. If an error is
received it will be captured here. It can be accessed the same way as
the image plots.

Forecasts are stored in a slightly different location if nothing has
been changed by user inputs `cd ./output/data/master/forecasts/plots`.
