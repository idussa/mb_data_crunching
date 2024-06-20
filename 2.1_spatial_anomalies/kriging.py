from __future__ import annotations
from typing import Callable

import numpy as np
from scipy.spatial.distance import pdist, squareform
import pandas as pd
import pyproj
from pykrige.ok import OrdinaryKriging

#################################
# ESTIMATED CORRELATION FUNCTIONS
#################################

def ba_anom_spatialcorr(d: np.ndarray):
    """
    Spatial correlation of annual anomaly in mass balance (this study = Dussaillant et al., 2024).

    :param d: Distance between two glaciers (meters).

    :return: Spatial correlation function (input = distance in meters, output = correlation between 0 and 1).
    """

    # Three ranges and partial sill for three exponential models
    r1 = 5.051675e+03
    r2 = 9.985127e+04
    r3 = 5.000000e+06
    ps1 = 0.224308
    ps2 = 0.140278
    ps3 = 0.635414

    exp1 = ps1 * (1 - np.exp(-3 * d / r1))
    exp2 = ps2 * (1 - np.exp(-3 * d / r2))
    exp3 = ps3 * (1 - np.exp(-3 * d / r3))

    # Spatial correlation
    return 1 - (exp1 + exp2 + exp3)

############################################################################
# COORDINATE TRANSFORMATIONS (FOR SPEED COMPUTING DISTS DURING ERROR PROPAG)
############################################################################

def latlon_to_utm(lat: float, lon: float) -> str:
    """
    Get UTM zone for a given latitude and longitude coordinates.

    :param lat: Latitude coordinate.
    :param lon: Longitude coordinate.

    :returns: UTM zone.
    """

    if not (
        isinstance(lat, (float, np.floating, int, np.integer))
        and isinstance(lon, (float, np.floating, int, np.integer))
    ):
        raise TypeError("Latitude and longitude must be floats or integers.")

    if not -180 <= lon < 180:
        raise ValueError("Longitude value is out of range [-180, 180[.")
    if not -90 <= lat < 90:
        raise ValueError("Latitude value is out of range [-90, 90[.")

    # Get UTM zone from name string of crs info
    utm_zone = pyproj.database.query_utm_crs_info(
        "WGS 84", area_of_interest=pyproj.aoi.AreaOfInterest(lon, lat, lon, lat)
    )[0].name.split(" ")[-1]

    return str(utm_zone)


def utm_to_epsg(utm: str) -> int:
    """
    Get EPSG code of UTM zone.

    :param utm: UTM zone.

    :return: EPSG of UTM zone.
    """

    if not isinstance(utm, str):
        raise TypeError("UTM zone must be a str.")

    # Whether UTM is passed as single or double digits, homogenize to single-digit
    utm = str(int(utm[:-1])) + utm[-1].upper()

    # Get corresponding EPSG
    epsg = pyproj.CRS(f"WGS 84 / UTM Zone {utm}").to_epsg()

    return int(epsg)


def reproject_points(
    points: list[list[float]] | list[float] | tuple[list[float], list[float]] | np.ndarray, in_crs: pyproj.CRS, out_crs: pyproj.CRS
) -> tuple[list[float], list[float]]:
    """
    Reproject a set of point from input_crs to output_crs.

    :param points: Input points to be reprojected. Must be of shape (2, N), i.e (x coords, y coords)
    :param in_crs: Input CRS
    :param out_crs: Output CRS

    :returns: Reprojected points, of same shape as points.
    """
    assert np.shape(points)[0] == 2, "points must be of shape (2, N)"

    x, y = points
    transformer = pyproj.Transformer.from_crs(in_crs, out_crs)
    xout, yout = transformer.transform(x, y)
    return (xout, yout)


def reproject_from_latlon(
    points: list[list[float]] | tuple[list[float], list[float]] | np.ndarray, out_crs: pyproj.CRS, round_: int = 2
) -> tuple[list[float], list[float]]:
    """
    Reproject a set of point from lat/lon to out_crs.

    :param points: Input points to be reprojected. Must be of shape (2, N), i.e (x coords, y coords)
    :param out_crs: Output CRS
    :param round_: Output rounding. Default of 2 ensures cm accuracy

    :returns: Reprojected points, of same shape as points.
    """
    crs_4326 = pyproj.CRS.from_epsg(4326)
    proj_points = reproject_points(points, crs_4326, out_crs)
    proj_points = np.round(proj_points, round_)
    return proj_points


##################################################
# SPATIAL INTERPOLATION WITH CORRELATION = KRIGING
##################################################

def krige_ba_anom(xobs: np.ndarray, yobs: np.ndarray, ba_anom_obs: np.ndarray, xpred: np.ndarray, ypred: np.ndarray):
    """
    Interpolate annual mass balance anomaly using kriging.

    :param xobs: X coordinates of observed glaciers.
    :param yobs: Y coordinates of observed glaciers.
    :param ba_anom_obs: Annual mass balance anomalies of observed glaciers.
    :param xpred: X coordinates of glaciers to predict.
    :param ypred: Y coordinates of glaciers to predict.

    :return: Annual mass balance anomalies of glacier to predict, Error (1-sigma) of predicted anomalies.
    """

    # TODO: Either pass standardized anomalies, or multiply model by variance here

    # Ordinary kriging = kriging with a mean function
    OK = OrdinaryKriging(
        xobs,
        yobs,
        ba_anom_obs,
        variogram_model="custom",
        variogram_parameters=[],
        variogram_function=ba_anom_spatialcorr,
        verbose=False,
        enable_plotting=False,
    )

    # Predict on grid, with uncertainty
    ba_anom_pred, sig_anom_pred = OK.execute("points", xpred, ypred)

    return ba_anom_pred, sig_anom_pred