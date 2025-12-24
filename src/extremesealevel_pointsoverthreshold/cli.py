import click
import logging

from extremesealevel_pointsoverthreshold.extremesealevel_pointsoverthreshold_preprocess import (
    esl_preprocess,
)

from extremesealevel_pointsoverthreshold.extremesealevel_pointsoverthreshold_fit import (
    extremesl_fit,
)
from extremesealevel_pointsoverthreshold.extremesealevel_pointsoverthreshold_project import (
    extremesl_project,
    project_prep,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


@click.command()
@click.option(
    "--total-localsl-file",
    type=click.Path(exists=True),
    required=True,
    help="Total localized sea-level projection file. Site lats/lons are taken from this file and mapped to the GESLA database [default='./total-workflow_localsl.nc'",
)
@click.option(
    "--min-days",
    type=int,
    help="Minimum number of days in a valid year",
    default=250,
    show_default=True,
)
@click.option(
    "--min-years",
    type=int,
    help="Minimum number of years available",
    default=20,
    show_default=True,
)
@click.option(
    "--match-lim",
    type=float,
    help="Radius around requested locations to find a matching tide gauge in GESLA database",
    default=0.1,
    show_default=True,
)
@click.option(
    "--center-year",
    help="This year +/- 9 years for centering data must be available",
    type=int,
    default=2000,
    show_default=True,
)
@click.option(
    "--pct-pot",
    help="Percentile for Point Over Threshold analysis",
    default=95,
    type=int,
    show_default=True,
)
@click.option(
    "--gpd-pot-threshold",
    help="Percentile for GPD analysis ",
    type=float,
    default=99.7,
    show_default=True,
)
@click.option(
    "--cluster-lim",
    help="Maximum number of hours that define a cluster for extreme events",
    default=72,
    type=int,
    show_default=True,
)
@click.option(
    "--target-years",
    help="Comma-delimited list of projection years of interest (e.g. 2050,2100)",
    default="2100",
    type=str,
    show_default=True,
)
@click.option(
    "--gesla-dir",
    help="Directory containing GESLA data base",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--pipeline-id",
    help="Unique identifier for this instance of the module",
    type=str,
    default="my_pipeline_id",
)
@click.option(
    "--min-z",
    help="Minimum heigh over which frequency is calculated (m)",
    default=-0.5,
    type=float,
    show_default=True,
)
@click.option(
    "--max-z",
    help="Maximum heigh over which frequency is calculated (m)",
    default=8.0,
    show_default=True,
    type=float,
)
@click.option(
    "--quantile-min",
    help="Minimum quantile to assess",
    default=0.01,
    show_default=True,
    type=float,
)
@click.option(
    "--quantile-max",
    help="Maximum quantile to assess",
    default=0.99,
    show_default=True,
    type=float,
)
@click.option(
    "--quantile-step",
    help="Quantile step",
    default=0.01,
    show_default=True,
    type=float,
)
@click.option(
    "--z-step",
    help="Stepping from min-z to max-z",
    default=0.01,
    type=float,
    show_default=True,
)
@click.option(
    "--nsamps",
    help="Number of samples to draw",
    default=500,
    type=int,
    show_default=True,
)
@click.option(
    "--allowance-freq",
    help="Frequency at which allowances are calculated",
    default=0.01,
    type=float,
    show_default=True,
)
@click.option(
    "--seed",
    help="Seed for the random number generator",
    default=1234,
    type=int,
    show_default=True,
)
@click.option(
    "--output-dir",
    type=click.Path(),
    help="Directory to write output files to",
    required=False,
)
@click.option(
    "--debug/--no-debug",
    default=False,
    envvar="EXTREMESEALEVEL_POINTSOVERTHRESHOLD_DEBUG",
)
def main(
    total_localsl_file,
    min_days,
    min_years,
    match_lim,
    center_year,
    pct_pot,
    gpd_pot_threshold,
    cluster_lim,
    target_years,
    gesla_dir,
    pipeline_id,
    min_z,
    max_z,
    quantile_min,
    quantile_max,
    quantile_step,
    z_step,
    nsamps,
    allowance_freq,
    seed,
    output_dir,
    debug,
):
    if debug:
        logging.root.setLevel(logging.DEBUG)
    else:
        logging.root.setLevel(logging.INFO)
    # make target years list
    target_years = [int(yr) for yr in target_years.split(",")]

    logger.info("Starting esl preprocessing...")
    slproj_output, config_output, station_data = esl_preprocess(
        total_localsl_file,
        min_days,
        min_years,
        match_lim,
        center_year,
        pct_pot,
        gpd_pot_threshold,
        cluster_lim,
        target_years,
        gesla_dir,
        pipeline_id,
    )
    logger.info("Finished esl preprocessing.")

    fitted_data = extremesl_fit(station_data=station_data, pipeline_id=pipeline_id)
    logger.info("Finished esl fitting.")
    proj_qnts = project_prep(
        quantile_min=quantile_min,
        quantile_max=quantile_max,
        quantile_step=quantile_step,
    )
    logger.info("Starting esl projection...")
    extremesl_project(
        esl_fit_file=fitted_data,
        slproj_data=slproj_output,
        min_z=min_z,
        max_z=max_z,
        z_step=z_step,
        allowance_freq=allowance_freq,
        nsamps=nsamps,
        seed=seed,
        proj_qnts=proj_qnts,
        pipeline_id=pipeline_id,
        # station_data=station_data,
        output_dir=output_dir,
    )
    logger.info("Finished esl projection.")
