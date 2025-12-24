import os
import numpy as np
import pandas as pd
from datetime import datetime as dt
from netCDF4 import Dataset
from typing import Tuple, Any, Optional

# Do not warn about chained assignments
pd.options.mode.chained_assignment = None  # default='warn'


""" extremesealevel_preprocess_pointsoverthreshold.py

Runs the pre-processing stage for the default extreme sea-level workflow.

Parameters:
total_localsl_file = Total localized sea-level projection file. Site lats/lons are taken from this file and mapped to the GESLA database [default="./total-workflow_localsl.nc"]
minDays = Minimum number of days in a valid year [default=250]
minYears = Minimum number of years available [default=20]
match_lim = Radius around requested locations to find a matching tide gauge in GESLA database [default=0.1]
center_year = This year +/- 9 years for centering data must be available [default=2000]
pctPot = Percentile for Point Over Threshold analysis [default=95]
gpd_pot_threshold = Percentile for GPD analysis [default=99.7]
cluster_lim = Maximum number of hours that define a cluster for extreme events [default=72]
target_years = Space-delimited list of projection years of interest (e.g. 2050 2100) [default=2100]
gesla_dir = Directory containing GESLA database [default="./gesla_data"]
pipeline_id = Unique identifier for this instance of the module [default="1234test"]

"""


def readmeta(filename: str) -> Tuple[str, float, float]:
    with open(filename, encoding="raw_unicode_escape") as myfile:
        head = [next(myfile) for x in range(6)]
    station_name = head[1][12:-1].strip()
    station_lat = float(head[4][11:-1].strip())
    station_lon = float(head[5][12:-1].strip())

    return (station_name, station_lat, station_lon)


def extract_gesla_locations(gesladir: str) -> Tuple[list, list, list, list]:
    # Get a list of the gesla database files
    geslafiles = os.listdir(gesladir)

    # Initialize the station variables
    station_names = []
    station_lats = []
    station_lons = []
    station_filenames = []

    # Loop over the gesla files
    for this_file in geslafiles:
        # Extract the station header information
        this_name, this_lat, this_lon = readmeta(os.path.join(gesladir, this_file))

        # Append this information to appropriate lists
        station_names.append(this_name)
        station_lats.append(float(this_lat))
        station_lons.append(float(this_lon))
        station_filenames.append(this_file)

    return (station_names, station_lats, station_lons, station_filenames)


def angd(lat0: list[Any], lon0: list[Any], lat: float, lon: float) -> float:
    # Convert the input from degrees to radians
    (lat0, lon0) = np.radians((lat0, lon0))
    (lat, lon) = np.radians((lat, lon))

    # Calculate the angle between the vectors
    temp = np.arctan2(
        np.sqrt(
            (np.cos(lat) * np.sin(lon - lon0)) ** 2
            + (
                np.cos(lat0) * np.sin(lat)
                - np.sin(lat0) * np.cos(lat) * np.cos(lon - lon0)
            )
            ** 2
        ),
        (np.sin(lat0) * np.sin(lat) + np.cos(lat0) * np.cos(lat) * np.cos(lon - lon0)),
    )

    # Convert the results from radians to degrees and return
    return np.degrees(temp)


def mindist(
    qlat: float, qlon: float, lats: list, lons: list, limit=0.1
) -> Optional[list]:
    # Calculate the angular distance
    dist = angd(lats, lons, qlat, qlon)
    # Log distance calculation results

    # If the minimum distance is beyond the limit, print a warning and return None
    if np.amin(dist) > limit:
        return None

    else:
        # Perform an indirect sort of the distances
        sort_idx = np.argsort(dist)
        # Find which ones fall within the radius limit
        min_idx = sort_idx[np.flatnonzero(dist[sort_idx] <= limit)]
        return min_idx


def load_tg(
    geslafile: str,
    minDays: int,
    minYears: int,
    center_year: int,
    pctPot: float,
    gpd_pot_threshold: float,
    cluster_lim: int,
):
    # Initialize data constraint flags
    # include_file = False
    pass_nyears = False
    pass_centeryears = False

    # read meta data
    with open(geslafile, encoding="raw_unicode_escape") as myfile:
        head = [next(myfile) for x in range(6)]
    station_name = head[1][12:-1].strip()
    station_lat = float(head[4][11:-1].strip())
    station_lon = float(head[5][12:-1].strip())

    # read raw data
    rawdata = np.loadtxt(
        geslafile,
        skiprows=32,
        dtype={
            "names": ("date", "hour", "val", "qf", "ua"),
            "formats": ("S10", "S8", "f4", "i4", "i4"),
        },
        encoding="raw_unicode_escape",
    )

    # remove invalid or missing data (see also GESLA2 definitions)
    rawdata = rawdata[rawdata["val"] > -9.9]  # skip remaining missing data
    rawdata = rawdata[rawdata["qf"] == 1]  # data with quality-flag 1 is correct

    # store date and time in datetime object with hourly resolution
    time_hourly = [
        dt(
            year=int(mydate[0:4]),
            month=int(mydate[5:7]),
            day=int(mydate[8:10]),
            hour=int(myhour[0:2]),
        )
        for mydate, myhour in zip(rawdata["date"], rawdata["hour"])
    ]

    # compute hourly means as means of all time entries starting with same hour
    df = pd.DataFrame({"mytime": time_hourly, "height": rawdata["val"]})
    data_hourly = df.groupby("mytime", as_index=False).mean()

    # unique years in datetimes
    unique_years = np.array(list({(i.year) for i in data_hourly["mytime"]}))

    # calculate annual means
    annual_means = data_hourly.groupby(data_hourly.mytime.dt.year).mean()

    # count number of observations in each year
    obs_pyear = data_hourly.groupby(data_hourly.mytime.dt.year, as_index=False).count()

    # select only years with enough obs
    goodYears_ind = (
        obs_pyear > minDays * 24
    )  # nad hoc, no inferences made about how these hours are distributed in time
    good_years = unique_years[goodYears_ind.mytime.values]
    nYears = len(good_years)

    # Make sure there are 19 good year surrounding year 2000
    center_years = np.arange(center_year - 9, center_year + 10)
    if np.sum(np.isin(center_years, good_years)) >= 16:
        pass_centeryears = True

    # Does this file pass the minimum data check
    if nYears >= minYears:
        pass_nyears = True

    # Does everything pass?
    if not (pass_centeryears and pass_nyears):
        raise Exception("File does not pass data constraint")

    # loop over unique years with enough observations and compute anomalies wrt annual mean
    good_data = []
    for myyear in good_years:
        yearly_anom = data_hourly[data_hourly["mytime"].dt.year == myyear]
        yearly_anom.height -= annual_means.loc[myyear].height  # subtract annual means

        good_data.append(yearly_anom)

    good_data = pd.concat(good_data)
    good_data.reset_index(drop=True, inplace=True)

    # compute Peak-over-Threshold (POT) height and fetch observed heights exceeding POT
    pot_threshold = np.percentile(good_data.height.values, pctPot)

    extremes = good_data[good_data["height"] > pot_threshold]
    extremes.reset_index(drop=True, inplace=True)
    extremes = extremes.sort_values(
        by=["height"], ascending=False, kind="mergesort"
    )  # note mergesort required to get same output as matlab, this has stable ordered duplicate values

    # compute range of pots
    pot_vals = np.percentile(good_data.height.values, gpd_pot_threshold)

    ##### DECLUSTER EXTREME EVENTS
    # initialize declustered extremes with first extreme
    decl_extremes = extremes.iloc[
        [0]
    ]  # double [] makes pandas dataframe instead of series

    # check iteratively if new extreme is too close (in time) to declustered extremes
    for i in range(1, len(extremes)):
        next_extreme = extremes.iloc[[i]]  # pick new extreme
        next_extreme.reset_index(drop=True, inplace=True)  # reset indices

        timediffs = (
            decl_extremes["mytime"] - next_extreme.mytime[0]
        )  # calculate time difference between declustered extremes and next extreme

        # if none of the current extremes are less than 'cluster_lim' hours away from next extreme
        if all(abs(timediffs / np.timedelta64(3600, "s")) >= cluster_lim):
            decl_extremes = pd.concat([decl_extremes, next_extreme])  # add next extreme

    decl_extremes.reset_index(drop=True, inplace=True)

    # calculate return periods in years
    decl_extremes["return_period"] = [
        (1 + len(good_data) / 365.25 / 24) / i for i in range(1, len(decl_extremes) + 1)
    ]

    # Return the information for this station
    return (
        station_name,
        station_lat,
        station_lon,
        pot_vals,
        nYears,
        good_data,
        decl_extremes,
    )


def preprocess_gesla(
    gesladir: str,
    minDays: int,
    minYears: int,
    match_limit: float,
    center_year: int,
    pctPot: float,
    gpd_pot_threshold: float,
    cluster_lim: int,
    site_lats: list,
    site_lons: list,
    site_ids: list,
    pipeline_id: str,
):
    # Extract the gesla station information
    (station_names, station_lats, station_lons, station_filenames) = (
        extract_gesla_locations(gesladir)
    )

    # Match the target lats/lons with their nearest stations
    min_idx = [
        mindist(x, y, station_lats, station_lons, match_limit)
        for x, y in zip(site_lats, site_lons)
    ]

    # Generate a list of input files and matched IDSs for the matched tide gauge data
    matched_filenames = []
    matched_ids = []
    for i in np.arange(len(min_idx)):
        if min_idx[i] is not None:
            matched_filenames.append([station_filenames[x] for x in min_idx[i]])
            matched_ids.append(site_ids[i])

    # If no matches are found, quit with error
    if not matched_filenames:
        raise ValueError(
            "No matches found within {} degrees for provided lat/lon list".format(
                match_limit
            )
        )

    # Define the output directory
    # outdir = os.path.dirname(__file__)

    # Initialize station data dictionary to hold matched location data
    station_data = {}

    # Initialize variables to track the files that have been tested
    pass_files = {}
    fail_files = []

    # Loop over the matched files
    for i in np.arange(len(matched_filenames)):
        # This ID
        this_id = matched_ids[i]
        this_id_passed = False

        # Loop over the files within the match radius for this location
        for this_file in matched_filenames[i]:
            # This file was tested and passed
            if np.isin(this_file, list(pass_files.keys())):
                print(
                    "{0} already PASSED a previous check on data constraints. Mapping site ID {1} to site ID {2}.".format(
                        this_file, pass_files[this_file], this_id
                    )
                )
                station_data[this_id] = station_data[pass_files[this_file]]
                this_id_passed = True

            # If this ID already has data, skip to the next ID
            if this_id_passed:
                continue

            # Move on to the next file if this file was already tested and failed
            if np.isin(this_file, fail_files):
                print(
                    "{0} already FAILED previous check on data constraints. Continuing with the next file".format(
                        this_file
                    )
                )
                continue

            # This file has not been tested yet, go ahead and try to load it in.
            try:
                # Load this tide gauge data
                (
                    match_name,
                    match_lat,
                    match_lon,
                    pot_vals,
                    nyears,
                    obs,
                    decl_extremes,
                ) = load_tg(
                    os.path.join(gesladir, this_file),
                    minDays,
                    minYears,
                    center_year,
                    pctPot,
                    gpd_pot_threshold,
                    cluster_lim,
                )

                # Put this entry into the station data dictionary
                station_data[this_id] = {
                    "station_name": match_name,
                    "station_id": this_id,
                    "lon": match_lon,
                    "lat": match_lat,
                    "pot_vals": pot_vals,
                    "nyears": nyears,
                    "obs": obs,
                    "decl_extremes": decl_extremes,
                }

                # This file passed, add it to the pass file dictionary
                pass_files[this_file] = this_id

                # This ID has data, set the flag to move onto the next ID
                this_id_passed = True

            except ValueError as e:
                # This file failed, add it to the fail list and continue with the next file
                print(
                    "{} did not pass the data constraints: {}. Moving on to the next file.".format(
                        this_file, e
                    )
                )
                fail_files.append(this_file)
                continue

        # Let the user know we didn't find a file that passes the data constraints
        if not this_id_passed:
            print(
                "No locations within {0} degrees pass the data constraints for ID {1}.".format(
                    match_limit, this_id
                )
            )

    # Exit with error if we didn't find any gesla location suitable
    if not station_data:
        raise ValueError("No data found for any requested locations")

    return (station_data, list(station_data.keys()))


def esl_preprocess(
    total_localsl_file: str,
    min_days: int,
    min_years: int,
    match_lim: float,
    center_year: int,
    pct_pot: float,
    gpd_pot_threshold: float,
    cluster_lim: int,
    target_years: list,
    gesla_dir: str,
    pipeline_id: str,
):
    # total_sl_ds = xr.open_dataset(total_localsl_file)

    # downcast to float32 so all match (rn some modules have float64 some float32)
    # total_sl_ds['lat'] = total_sl_ds['lat'].astype('float32')
    # total_sl_ds['lon'] = total_sl_ds['lon'].astype('float32')
    # this is just a stand in for now when there's one loc
    # will need to change when there are multiple locations in the output
    # and do this within every location (there's a much cleaner way to do
    # this but may need to wait for a bigger refactor)
    # for loc in total_sl_ds['locations'].values:
    # for coord in ['lat','lon']:
    # assert len(np.unique(total_sl_ds[coord].sel(locations=loc).values)) == 1, "Something's wrong and the lat lon values do not match for multiple module output files for the same location."

    # site_lats = total_sl_ds['lat'].data

    # site_lons = total_sl_ds['lon'].data
    # site_ids = total_sl_ds['locations'].data
    # proj_yrs = total_sl_ds['years'].data
    # proj_slc = total_sl_ds['totaled_sea_level_change'].data
    # Load the sea-level projection file
    nc = Dataset(total_localsl_file, "r")
    # pdb.set_trace()
    # Extract data from NetCDF file
    site_lats = nc.variables["lat"][:].data
    site_lons = nc.variables["lon"][:].data
    site_ids = nc.variables["locations"][:].data
    proj_yrs = nc.variables["years"][:].data
    # proj_qnts = nc.variables['quantiles'][:].data
    # proj_slc_qnts = nc.variables['localSL_quantiles'][::,::,::].data
    proj_slc = nc.variables["sea_level_change"][::, ::, ::].data
    nc.close()

    # Make sure lats and lons are equal in length
    if not len(site_lats) == len(site_lons):
        raise ValueError(
            "Number of latitudes and longitudes not equal ({0} != {1})".format(
                len(site_lats), len(site_lons)
            )
        )

    # Test to make sure the list of site ids matches the length of the lats/lons
    if not len(site_lats) == len(site_ids):
        raise ValueError(
            "Number of site IDs not equal to number of locations provided in lat/lon list ({0} != {1})".format(
                len(site_ids), len(site_lats)
            )
        )

    # Find the target years that overlap the projection years
    slproj_targ_year_idx = np.flatnonzero(np.isin(proj_yrs, target_years))
    if len(slproj_targ_year_idx) == 0:
        raise ValueError("Cannot find target years in projection years")

        # Run the gesla preprocessing and return the matching side IDs
    station_data, matched_ids = preprocess_gesla(
        gesla_dir,
        min_days,
        min_years,
        match_lim,
        center_year,
        pct_pot,
        gpd_pot_threshold,
        cluster_lim,
        site_lats,
        site_lons,
        site_ids,
        pipeline_id,
    )
    slproj_matched_ids_idx = np.flatnonzero(np.isin(site_ids, matched_ids))

    slproj_output = {}
    for this_id_idx in slproj_matched_ids_idx:
        proj_slc_subset = proj_slc[::, [slproj_targ_year_idx], this_id_idx]

        slproj_output[site_ids[this_id_idx]] = {
            "site_lat": site_lats[this_id_idx],
            "site_lon": site_lons[this_id_idx],
            "site_id": site_ids[this_id_idx],
            "proj_years": proj_yrs[slproj_targ_year_idx],
            "proj_slc": proj_slc_subset,
        }

    config_output = {
        "minDays": min_days,
        "minYears": min_years,
        "match_lim": match_lim,
        "center_year": center_year,
        "pctPot": pct_pot,
        "cluster_lim": cluster_lim,
        "total_localsl_file": total_localsl_file,
        "target_years": target_years,
        "gesla_dir": gesla_dir,
        "pipeline_id": pipeline_id,
        "gpd_pot_threshold": gpd_pot_threshold,
    }

    return slproj_output, config_output, station_data
