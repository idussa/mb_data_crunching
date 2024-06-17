from typing import Iterable
import numpy as np


MEAN_EARTH_RADIUS: float = 6371.0
"""
Average radius of the Earth (km).

See https://en.wikipedia.org/wiki/Earth_radius.
"""


def great_circle_distance(
  lat0: float,
  lng0: float,
  lat: Iterable[float],
  lng: Iterable[float],
  method: str = 'vincenty'
) -> np.ndarray:
  """
  Calculate the great-circle distance from one to many points.

  Arguments
  ---------
  lat0
    Latitude of origin (decimal degrees).
  lng0
    Longitude of origin (decimal degrees).
  lat
    Latitude of destination points (decimal degrees).
  lng
    Longitude of destination points (decimal degrees).
  method
    Method to use for the computation:

    - 'cosines': Spherical law of cosines. Not precise for distances less than a
      few meters.
    - 'vincenty': Vincenty formula. Precise for small distances and antipodal
      points.

    See https://en.wikipedia.org/wiki/Great-circle_distance.

  Returns
  -------
  np.ndarray
    Distance (km) from origin `lat0, lng0` to each point in `lat, lng`.
    Assumes the average radius of the Earth, which actually varies by 0.3%.

  Examples
  --------
  >>> great_circle_distance(47.4, -122.3, [47.4], [-122.3])
  array([0.])
  >>> great_circle_distance(47.4, -122.3, [47.4, 47.41], [-122.3, -122.29])
  array([0., 1.34268696])
  """
  # Coerce to numpy arrays (as needed)
  lat, lng = np.asarray(lat), np.asarray(lng)
  # Convert from degrees to radians
  lat0, lng0 = np.radians(lat0), np.radians(lng0)
  lat, lng = np.radians(lat), np.radians(lng)
  # Compute central angle
  if method == 'cosines':
    # Use the spherical law of cosines
    # https://en.wikipedia.org/wiki/Great-circle_distance#Formulae
    theta = np.arccos(
      np.sin(lat0) * np.sin(lat) +
      np.cos(lat0) * np.cos(lat) * np.cos(np.abs(lng - lng0))
    )
  elif method == 'vincenty':
    # Use the Vincenty formula
    # https://en.wikipedia.org/wiki/Great-circle_distance#Computational_formulas
    sin_lat0, cos_lat0 = np.sin(lat0), np.cos(lat0)
    sin_lat, cos_lat = np.sin(lat), np.cos(lat)
    delta_lng = lng - lng0
    cos_delta_lng, sin_delta_lng = np.cos(delta_lng), np.sin(delta_lng)
    theta = np.arctan2(
      np.sqrt(
        (cos_lat * sin_delta_lng) ** 2 +
        (cos_lat0 * sin_lat - sin_lat0 * cos_lat * cos_delta_lng) ** 2
      ),
      sin_lat0 * sin_lat + cos_lat0 * cos_lat * cos_delta_lng
    )
  return MEAN_EARTH_RADIUS * theta
