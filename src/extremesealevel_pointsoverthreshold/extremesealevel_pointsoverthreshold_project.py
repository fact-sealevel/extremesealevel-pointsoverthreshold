#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""project_esl.py
This script runs the projecting task for extreme sea-level analysis. This tasks
generates samples of local msl change and GPD parameters. From those samples it
calculates historical and future return curves at user defined return periods.
The return curves are used to calculate the amplification factor and allowance
for a given station at user defined percentiles.

The analysis is based on the MATLAB code of Thomas Frederikse used for SROCC
and LocalizeSL (Buchanan et al. 2016). A Peak-Over-Threshold is used.
Above the threshold a Pareto distribution is used to model the extremes. Below,
a Gumbel distribution is assumed and cut off at MHHW. MHHW is calculated as the
long-term mean of daily maxima.

Input parameters:
station_data		= dictionary with relevant preprocessed GESLA 2 data
proj_slc	 		= gridded sea-level change samples
allowance_freq		= query n-yr event frequency
testz				= query heights
num_mc				= number of samples to draw

Output: Adds results of extreme sea-level analysis for all stations to the
station data dictionary.
lcl_slc_median		= median sea-level change at station
lcl_slc_qnts		= sea-level change quantiles at station
amp_factor_qnts		= amplification factor percentiles
allowance_qnts		= allowance percentiles
fut_freqs_qnts		= future return curves percentiles
hist_freqs			= historical return curve

Created on Tue Nov	5 09:37:57 2019
@author: Tim Hermans
tim(dot)hermans@nioz(dot)nl
"""

import argparse
import os
import time
import numpy.matlib
import numpy as np
import xarray as xr
import fnmatch
import tarfile
import logging

# from sample_from_quantiles import sample_from_quantiles
# import pandas as pd
# from datetime import datetime as dt
logger = logging.getLogger(__name__)


def getFreqFromZ_ESL(scale, shape, loc, avg_exceed, z, mhhw, mhhwFreq):
    # This function obtains frequencies for corresponding test heights.
    # It uses a General Pareto Distribution for heights that exceed the threshold.
    # A Gumbel distribution is assumed between the MHHW and the threshold.

    # input:
    # scale: scale factor GPD
    # shape: shape factor GPD
    # loc: threshold historical
    # avg_exceed: frequency that threshold is exceeded
    # z: test height minus (historical or future) location parameter
    # mhhw: Mean Higher High Water
    # mhhwFreq: frequency of MHHW

    # output:
    # Return frequencies

    z0 = z
    z0[z0 < (mhhw - loc)] = mhhw - loc  # cut off heights below loc at MHHW

    idx = np.ones(np.shape(z), dtype=bool)
    idx[z < 0] = (
        False  # indices where test heights are below loc (for large SLR or low testz)
    )
    idx[(shape * z / scale) < -1] = (
        False  # above upper bound GPD (for low/negative SLR or high testz)
    )

    sub = np.zeros(np.shape(z), dtype=bool)
    sub[z < 0] = True
    freq = np.nan * z  # initialize
    log_freq = np.nan * z

    if np.isscalar(scale) and np.isscalar(shape):
        freq[idx] = avg_exceed * np.power(
            (1 + (shape * z[idx] / scale)), (-1 / shape)
        )  # get frequencies from pareto distribution for testz>loc
    else:
        freq[idx] = avg_exceed * np.power(
            (1 + (shape[idx] * z[idx] / scale[idx])), (-1 / shape[idx])
        )  # get frequencies from pareto distribution for testz>loc

    log_freq[sub] = (
        np.log(avg_exceed)
        + (np.log(mhhwFreq) - np.log(avg_exceed)) * (z0[sub] / (mhhw - loc))
    )  # from Gumbel (see Buchanan et al. 2016 & LocalizeSL) for testz<loc (logN linear in z, assume MHHW once every 2 days)
    freq[sub] = np.exp(log_freq[sub])

    # freq[np.less(freq, 1e-6, where=~np.isnan(freq))] = np.nan #throw away frequencies lower than 1e-6 (following SROCC scripts)

    return freq


def create_extreme_sealevel_dataset_xr(
    testz,
    proj_qnts,
    nsamps,
    site_lat,
    site_lon,
    site_id,
    proj_years,
    lcl_msl_samples,
    loc,
    mhhw,
    mhhwFreq,
    fut_freqs_qnts,
    hist_freqs_qnts,
    ampfactors_qnts,
    allowances_qnts,
    shape_samples,
    scale_samples,
    seed,
    allowance_freq,
):
    """
    Create an xarray Dataset equivalent to the netCDF4 Dataset creation.

    Parameters match the variables used in the netCDF4 version.
    """

    # Create coordinate arrays
    heights_coord = testz
    quantiles_coord = proj_qnts
    samples_coord = np.arange(nsamps)
    scalar_coord = np.array([0])  # Single element for scalar dimension

    # Create data variables as dictionaries
    data_vars = {
        # Scalar variables (1D with scalar dimension)
        "lat": (["scalar"], np.array([site_lat])),
        "lon": (["scalar"], np.array([site_lon])),
        "id": (["scalar"], np.array([site_id], dtype=np.int32)),
        "year": (["scalar"], np.array([proj_years], dtype=np.int32)),
        "gpd_location": (["scalar"], np.array([loc])),
        "mhhw": (["scalar"], np.array([mhhw])),
        "mhhw_freq": (["scalar"], np.array([mhhwFreq])),
        # Data variables
        "localSL_quantiles": (["quantiles"], np.quantile(lcl_msl_samples, proj_qnts)),
        "projected_frequencies": (["quantiles", "heights"], fut_freqs_qnts),
        "historical_frequencies": (["quantiles", "heights"], hist_freqs_qnts),
        "amplification_factors_quantiles": (["quantiles"], ampfactors_qnts),
        "allowance_quantiles": (["quantiles"], allowances_qnts),
        "gpd_shape": (["samples"], shape_samples),
        "gpd_scale": (["samples"], scale_samples),
    }

    # Create coordinates dictionary
    coords = {
        "scalar": scalar_coord,
        "quantiles": quantiles_coord,
        "heights": heights_coord,
        "samples": samples_coord,
    }

    # Create the Dataset
    ds = xr.Dataset(
        data_vars=data_vars,
        coords=coords,
        attrs={
            "description": "Extreme Sea-Level",
            "history": f"Created {time.ctime(time.time())}, Seed {seed}",
            "source": f"FACTS - Extreme Sea-Level Module. Allowance Frequency = {allowance_freq}",
        },
    )

    # Add variable attributes (units)
    ds["lat"].attrs["units"] = "Degrees North"
    ds["lon"].attrs["units"] = "Degrees East"
    ds["localSL_quantiles"].attrs["units"] = "m"
    ds["heights"].attrs["units"] = "m"
    ds["gpd_location"].attrs["units"] = "m"
    ds["mhhw"].attrs["units"] = "m"
    ds["mhhw_freq"].attrs["units"] = "per annum"
    ds["projected_frequencies"].attrs["units"] = "per annum"
    ds["historical_frequencies"].attrs["units"] = "per annum"
    ds["allowance_quantiles"].attrs["units"] = "m"

    # Set encoding for compression (equivalent to zlib=True, least_significant_digit=6)
    encoding = {
        "localSL_quantiles": {
            "zlib": True,
            "least_significant_digit": 6,
        },
        "projected_frequencies": {
            "zlib": True,
            "least_significant_digit": 6,
        },
        "historical_frequencies": {
            "zlib": True,
            "least_significant_digit": 6,
        },
        "amplification_factors_quantiles": {
            "zlib": True,
            "least_significant_digit": 6,
        },
        "allowance_quantiles": {
            "zlib": True,
            "least_significant_digit": 6,
        },
        "gpd_shape": {
            "zlib": True,
            "least_significant_digit": 6,
        },
        "gpd_scale": {
            "zlib": True,
            "least_significant_digit": 6,
        },
    }

    # Save to NetCDF file with encoding
    # ds.to_netcdf(output_filename, encoding=encoding)

    return ds, encoding


def project_station(
    station_data,
    slproj_data,
    proj_qnts,
    testz,
    allowance_freq,
    nsamps,
    seed,
    output_filename,
):
    # Seed the RNG
    rng = np.random.default_rng(seed)

    # Extract the sea-level projection data
    lcl_msl_samples = slproj_data["proj_slc"] / 1000.0
    # should we clip the tails of the samples?

    proj_years = slproj_data["proj_years"]
    site_lat = slproj_data["site_lat"]
    site_lon = slproj_data["site_lon"]
    site_id = slproj_data["site_id"]

    # add SLC samples to location parameter
    logger.info(f"Station data keys: {list(station_data.keys())}")
    loc = station_data["loc"]
    loc_fut_samples = lcl_msl_samples + loc

    ###### GET GPD SAMPLES
    shape = station_data["gp_shape"]
    scale = station_data["gp_scale"]
    gp_cov = station_data["gp_cov"]  # covariance matrix
    avg_exceed = station_data[
        "avg_exceed"
    ]  # average number of exceedences of threshold per year (lambda)

    gp_samples = rng.multivariate_normal([shape, scale], gp_cov, size=nsamps)
    shape_samples = gp_samples[:, 0]
    scale_samples = gp_samples[:, 1]
    scale_samples[scale_samples < 0.001] = 0.001  # no negative scales

    # get MHHW for Gumbel distribution below Pareto
    mhhw = station_data["mhhw"]
    mhhwFreq = station_data["mhhwFreq"]
    # historical; use best estimate scale and shape factors
    # hist_freqs = getFreqFromZ_ESL(scale,shape,loc,avg_exceed,testz-loc,mhhw,mhhwFreq) #get frequencies for test heights
    # idx_alw_hist = np.nanargmin(abs(hist_freqs-allowance_freq)) #find the historical height with 'allowance_freq' frequency

    ###### QUERY FUTURE RETURN FREQUENCIES FOR TEST HEIGHTS
    # generate matrices for faster multiplication and division
    testz_mat = np.transpose(np.matlib.repmat(testz, nsamps, 1))  # test heights
    loc_fut_mat = np.matlib.repmat(
        loc_fut_samples, len(testz), 1
    )  # future location parameter (loc+slr)
    scale_mat = np.matlib.repmat(scale_samples, len(testz), 1)  # scale & shape samples
    shape_mat = np.matlib.repmat(shape_samples, len(testz), 1)

    hist_freqs_mat = getFreqFromZ_ESL(
        scale_mat, shape_mat, loc, avg_exceed, testz_mat - loc, mhhw, mhhwFreq
    )  # get frequencies for test heights
    fut_freqs_mat = getFreqFromZ_ESL(
        scale_mat, shape_mat, loc, avg_exceed, testz_mat - loc_fut_mat, mhhw, mhhwFreq
    )

    ###### CALCULATE ALLOWANCES AND AMPLIFICATION FACTORS
    idx_alwfreq_hist = np.nanargmin(
        abs(hist_freqs_mat - allowance_freq), axis=0
    )  # for each sample historical return curve, get index of event with user defined frequency (default 1/100yr)
    idx_alwfreq_fut = np.nanargmin(
        abs(fut_freqs_mat - allowance_freq), axis=0
    )  # similar for sample future return curves

    # get amplification factors (how more often will the historical allowance frequency event with testz x occur in the future)
    ampfactors = (
        fut_freqs_mat[idx_alwfreq_hist, np.arange(0, len(idx_alwfreq_hist))]
        / allowance_freq
    )

    # get allowances (how much should a structure be raised to maintain the same allowance frequency in the future)
    allowances = (
        testz_mat[idx_alwfreq_fut, np.arange(0, len(idx_alwfreq_fut))]
        - testz_mat[idx_alwfreq_hist, np.arange(0, len(idx_alwfreq_hist))]
    )

    # get quantiles from these
    hist_freqs_qnts = np.nanquantile(hist_freqs_mat, proj_qnts, axis=1)
    fut_freqs_qnts = np.nanquantile(fut_freqs_mat, proj_qnts, axis=1)
    allowances_qnts = np.nanquantile(allowances, proj_qnts)
    ampfactors_qnts = np.nanquantile(ampfactors, proj_qnts)

    ds, encoding = create_extreme_sealevel_dataset_xr(
        testz=testz,
        proj_qnts=proj_qnts,
        nsamps=nsamps,
        site_lat=site_lat,
        site_lon=site_lon,
        site_id=site_id,
        proj_years=proj_years,
        lcl_msl_samples=lcl_msl_samples,
        loc=loc,
        mhhw=mhhw,
        mhhwFreq=mhhwFreq,
        fut_freqs_qnts=fut_freqs_qnts,
        hist_freqs_qnts=hist_freqs_qnts,
        ampfactors_qnts=ampfactors_qnts,
        allowances_qnts=allowances_qnts,
        shape_samples=shape_samples,
        scale_samples=scale_samples,
        seed=seed,
        allowance_freq=allowance_freq,
    )
    ds.to_netcdf(output_filename, encoding=encoding)
    logger.info(f"Written extreme sea-level data to {output_filename}")
    # -------------------------------------------------------------------------------------
    # Write this station information out to a netCDF file --------------------------------
    # rootgrp = Dataset(output_filename, "w", format="NETCDF4")

    # Define Dimensions
    # nheights = len(testz)
    # nq = len(proj_qnts)
    # height_dim = rootgrp.createDimension("heights", nheights)
    # q_dim = rootgrp.createDimension("quantiles", nq)
    # sample_dim = rootgrp.createDimension("samples", nsamps)
    # scalar_dim = rootgrp.createDimension("scalar", 1)

    # Populate dimension variables
    # lat_var = rootgrp.createVariable("lat", "f4", ("scalar",))
    # lon_var = rootgrp.createVariable("lon", "f4", ("scalar",))
    # id_var = rootgrp.createVariable("id", "i4", ("scalar",))
    # year_var = rootgrp.createVariable("year", "i4", ("scalar",))
    # q_var = rootgrp.createVariable("quantiles", "f4", ("quantiles",))
    # heights_var = rootgrp.createVariable("heights", "f4", ("heights",))
    # loc_var = rootgrp.createVariable("gpd_location", "f4", ("scalar",))
    # mhhw_var = rootgrp.createVariable("mhhw", "f4", ("scalar",))
    # mhhwFreq_var = rootgrp.createVariable("mhhw_freq", "f4", ("scalar",))

    # Create a data variable
    # localslq = rootgrp.createVariable("localSL_quantiles", "f4", ("quantiles",), zlib=True, least_significant_digit=6)
    # futfreqsq = rootgrp.createVariable("projected_frequencies", "f4", ("quantiles","heights"), zlib=True, least_significant_digit=6)
    # histfreqsq = rootgrp.createVariable("historical_frequencies", "f4", ("quantiles","heights",), zlib=True, least_significant_digit=6)
    # ampfactorq = rootgrp.createVariable("amplification_factors_quantiles", "f4", ("quantiles",), zlib=True, least_significant_digit=6)
    # allowanceq = rootgrp.createVariable("allowance_quantiles", "f4", ("quantiles",), zlib=True, least_significant_digit=6)
    # shapesamps = rootgrp.createVariable("gpd_shape", "f4", ("samples",), zlib=True, least_significant_digit=6)
    # scalesamps = rootgrp.createVariable("gpd_scale", "f4", ("samples",), zlib=True, least_significant_digit=6)

    # Assign attributes
    # rootgrp.description = "Extreme Sea-Level"
    # rootgrp.history = "Created " + time.ctime(time.time()) + ", Seed {0}".format(seed)
    # rootgrp.source = "FACTS - Extreme Sea-Level Module. Allowance Frequency = {}".format(allowance_freq)

    # Variable units
    # lat_var.units = "Degrees North"
    # lon_var.units = "Degrees East"
    # localslq.units = "m"
    # heights_var.units = "m"
    # loc_var.units = "m"
    # mhhw_var.units = "m"
    # mhhwFreq_var.units = "per annum"
    # futfreqsq.units = "per annum"
    # histfreqsq.units = "per annum"
    # allowanceq.units = "m"

    # Put the data into the netcdf variables
    # lat_var[:] = site_lat
    # lon_var[:] = site_lon
    # id_var[:] = site_id
    # year_var[:] = proj_years
    # q_var[:] = proj_qnts
    # localslq[:] = np.quantile(lcl_msl_samples,proj_qnts)
    # heights_var[:] = testz
    # loc_var[:] = loc
    # mhhw_var[:] = mhhw
    # mhhwFreq_var[:] = mhhwFreq
    # futfreqsq[:,:] = fut_freqs_qnts
    # histfreqsq[:,:] = hist_freqs_qnts
    # ampfactorq[:] = ampfactors_qnts
    # allowanceq[:] = allowances_qnts
    # shapesamps[:] = shape_samples
    # scalesamps[:] = scale_samples

    # Close the netcdf
    # rootgrp.close()

    # Done
    return 0


def project_prep(quantile_min, quantile_max, quantile_step):
    proj_qnts = np.arange(quantile_min, quantile_max, quantile_step)
    return proj_qnts


def extremesl_project(
    esl_fit_file,  # esl_fit_data,
    slproj_data,
    min_z,
    max_z,
    z_step,
    allowance_freq,
    nsamps,
    seed,
    proj_qnts,
    pipeline_id,
    # station_data,
    output_dir,
):
    station_data = esl_fit_file
    logging.info("esl fit file keys: %s", list(station_data.keys()))

    # Make the test heights
    testz = np.arange(min_z, max_z, z_step)

    # Loop through the available stations
    seed_offset = 0
    for station_id in station_data:
        # Extract the necessary items for this station
        this_station = station_data[station_id]
        this_slproj = slproj_data[station_id]

        # Set the seed for this station
        this_seed = seed + (seed_offset * 10)
        seed_offset += 1

        # Loop through the target years
        for i in np.arange(len(this_slproj["proj_years"])):
            # proj_slc_qnts = np.quantile(this_slproj["proj_slc"], proj_qnts, axis=0)

            # Build dictionary for this sea-level projection year
            this_slproj_year = {
                "site_lat": this_slproj["site_lat"],
                "site_lon": this_slproj["site_lon"],
                "site_id": this_slproj["site_id"],
                "proj_years": this_slproj["proj_years"][i],
                "proj_qnts": proj_qnts,
                "proj_slc": this_slproj["proj_slc"][::, 0, i],
            }

            # Output file name for this station id and target year combination
            this_output_filename = os.path.join(
                output_dir,
                "{}_id{}_year{}_extremesl.nc".format(
                    pipeline_id, station_id, this_slproj["proj_years"][i]
                ),
            )

            # Run the projection for this station
            project_station(
                this_station,
                this_slproj_year,
                proj_qnts,
                testz,
                allowance_freq,
                nsamps,
                this_seed,
                this_output_filename,
            )

    # Collect all the extremesl.nc files into an archive that can be retrieved
    
    curdir_files = os.listdir(output_dir) #os.listdir(".")
    logger.info(f"Current directory files: {curdir_files}")
    archive_files = fnmatch.filter(curdir_files, "*_extremesl.nc")
    with tarfile.open("{}/{}_extremesl.tgz".format(output_dir, pipeline_id), "w:gz") as tar:
        for name in archive_files:
            tar.add(os.path.join(output_dir, name))

    # Done
    return 0


if __name__ == "__main__":
    # Initialize the command-line argument parser
    parser = argparse.ArgumentParser(
        description="Run the projection stage for the extreme sea-level workflow",
        epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)",
    )

    # Define the command line arguments to be expected
    parser.add_argument(
        "--esl_fit_file",
        help="Fitted station data produced from fitting stage",
        default=None,
    )
    parser.add_argument(
        "--slproj_file",
        help="Sea-level rise projections extracted from the preprocessing stage",
        default=None,
    )
    parser.add_argument(
        "--min_z",
        help="Minimum height over which frequency is calculated (m) [default = -0.5]",
        type=float,
        default=-0.5,
    )
    parser.add_argument(
        "--max_z",
        help="Minimum height over which frequency is calculated (m) [default = 8.0]",
        type=float,
        default=8.0,
    )
    parser.add_argument(
        "--quantile_min",
        help="Minimum quantile to assess [default = 0.01]",
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--quantile_max",
        help="Maximum quantile to assess [default = 0.99]",
        type=float,
        default=0.99,
    )
    parser.add_argument(
        "--quantile_step",
        help="Quantile step [default = 0.01]",
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--z_step",
        help="Stepping from min_z to max_z [default = 0.01]",
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--allowance_freq",
        help="Frequency at which allowances are calculated [defalut = 0.01; 1/100]",
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--nsamps",
        help="Number of samples to draw [default = 20000]",
        type=int,
        default=20000,
    )
    parser.add_argument(
        "--pipeline_id", help="Unique identifier for this instance of the module"
    )
    parser.add_argument(
        "--seed",
        help="Seed for the random number generator [default = 1234]",
        default=1234,
        type=int,
    )

    # Parse the arguments
    args = parser.parse_args()

    # Use default for station data file if necessary
    if args.esl_fit_file is None:
        esl_fit_file = os.path.join(
            os.path.dirname(__file__), "{}_fit.pkl".format(args.pipeline_id)
        )
    else:
        esl_fit_file = args.esl_fit_file

    # Use default for slr projection data if necessary
    if args.slproj_file is None:
        slproj_file = os.path.join(
            os.path.dirname(__file__), "{}_slproj_data.pkl".format(args.pipeline_id)
        )
    else:
        slproj_file = args.slproj_file

    proj_qnts = np.arange(args.quantile_min, args.quantile_max, args.quantile_step)

    # Run the fitting process on the
    extremesl_project(
        esl_fit_file,
        slproj_file,
        args.min_z,
        args.max_z,
        args.z_step,
        args.allowance_freq,
        args.nsamps,
        args.seed,
        proj_qnts,
        args.pipeline_id,
    )

    # Done
    exit()
