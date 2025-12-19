# extremesealevel-pointsoverthreshold
This module fits a General Pareto Distribution on preprocessed GESLA2 data. It also calculates the associated covariance matrix. Using input from the total module, it generates samples of local msl change and GPD parameters. From those samples it calculates historical and future return curves at user defined return periods. The return curves are used to calculate the amplification factor and allowance for a given station at user defined percentiles. The analysis is based on the MATLAB code of Thomas Frederikse used for SROCC and LocalizeSL (Buchanan et al. 2016). A Peak-Over-Threshold is used. Above the threshold a Pareto distribution is used to model the extremes. Below, a Gumbel distribution is assumed and cut off at MHHW. MHHW is calculated as the long-term mean of daily maxima.

Code by Tim Hermans

> [!CAUTION]
> This is a prototype. It is likely to change in breaking ways. It might delete all your data. Don't use it in production.

## Example

First, clone the repo, download required input data and prepare for the run, like:
```shell
git clone --single-branch --branch package git@github.com:e-marshall/extremesealevel-pointsoverthreshold`
# ^^ for now
# eventually:
git clone git@github.com:fact-sealevel/extremesealevel-pointsoverthreshold

mkdir -p ./data/input
curl -sL https://zenodo.org/record/7478192/files/extremesealevel_pointsoverthreshold_data.tgz | tar -zx -C ./data/input
#alt
# curl -sL https://zenodo.org/record/7478192/files/extremesealevel_pointsoverthreshold_fulldata.tgz | tar -zx -C ./data/input

mkdir -p ./data/output
```
Build docker container:

```shell
docker build -t extremesealevel-pointsoverthreshold .
```

Now run the container:
```shell
docker run --rm \
    -v /path/to/dir/with/totaled_local_output:/mnt/totaled_sealevel_data \
    -v /path/to/extremesealevel_input_data:/mnt/extremesealevel_data_in \
    extremesealevel-pointsoverthreshold \
    --total-localsl-file /mnt/totaled_sealevel_data/total_lslr_ns1000.nc \
    --gesla-dir /mnt/extremesealevel_data_in/gesla_data
```

## Features

```shell
Usage: extremesealevel-pointsoverthreshold [OPTIONS]

Options:
  --total-localsl-file PATH  Total localized sea-level projection file. Site
                             lats/lons are taken from this file and mapped to
                             the GESLA database [default='./total-
                             workflow_localsl.nc'  [required]
  --min-days INTEGER         Minimum number of days in a valid year  [default:
                             250]
  --min-years INTEGER        Minimum number of years available  [default: 20]
  --match-lim FLOAT          Radius around requested locations to find a
                             matching tide gauge in GESLA database  [default:
                             0.1]
  --center-year INTEGER      This year +/- 9 years for centering data must be
                             available  [default: 2000]
  --pct-pot INTEGER          Percentile for Point Over Threshold analysis
                             [default: 95]
  --gpd-pot-threshold FLOAT  Percentile for GPD analysis   [default: 99.7]
  --cluster-lim INTEGER      Maximum number of hours that define a cluster for
                             extreme events  [default: 72]
  --target-years TEXT        Comma-delimited list of projection years of
                             interest (e.g. 2050,2100)  [default: 2100]
  --gesla-dir PATH           Directory containing GESLA data base  [required]
  --pipeline-id TEXT         Unique identifier for this instance of the module
  --min-z FLOAT              Minimum heigh over which frequency is calculated
                             (m)  [default: -0.5]
  --max-z FLOAT              Maximum heigh over which frequency is calculated
                             (m)  [default: 8.0]
  --quantile-min FLOAT       Minimum quantile to assess  [default: 0.01]
  --quantile-max FLOAT       Maximum quantile to assess  [default: 0.99]
  --quantile-step FLOAT      Quantile step  [default: 0.01]
  --z-step FLOAT             Stepping from min-z to max-z  [default: 0.01]
  --nsamps INTEGER           Number of samples to draw  [default: 500]
  --allowance-freq FLOAT     Frequency at which allowances are calculated
                             [default: 0.01]
  --seed INTEGER             Seed for the random number generator  [default:
                             1234]
  --output-dir PATH          Directory to write output files to
  --debug / --no-debug
  --help                     Show this message and exit.
```

See the above by running:
```shell
docker run --rm extremesealevel-pointsoverthreshold --help
```

## Results
If this module runs successfully, it will write one NetCDF file for every element of the `--target-years` option. For example, if `--target-years = 2050`, one file will be written. If `--target-years = 2050,2100`, two will be written.

## Support 
Source code is available online at https://github.com/fact-sealevel/extremesealevel-pointsoverthreshold. This software is open source, available under the MIT license.

Please file issues in the issue tracker at https://github.com/extremesealevel-pointsoverthreshold/facts-total/issues.