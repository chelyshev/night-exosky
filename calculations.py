import numpy as np
from kelvin_to_rgb import convert_K_to_RGB
from PyAstronomy import pyasl
from trianglesolver import solve, degree

def apparent_to_absolute_magnitude(dist_pc, app_mag):
    # dist = 10.0 if dist_pc <= 0.0 else dist_pc
    v = 5.0 * np.log10(dist_pc)
    abs_mag = np.zeros(len(app_mag), dtype=float)
    for i in range(len(app_mag)):
        abs_mag[i] = app_mag[i] - v[i] + 5.0
    return abs_mag

def apparent_to_absolute_magnitude_t(dist_pc, app_mag):
    # dist = 10.0 if dist_pc <= 0.0 else dist_pc
    v = 5.0 * np.log10(dist_pc)
    return app_mag - v + 5.0

def absolute_magnitude_to_pseudosize(abs_mag):
    # Pseudo-luminosity. Usually L = L0 * 10^(-0.4*Mbol). We omit M0 and approximate Mbol = M
    pseudo_l = np.power(10, -0.4 * abs_mag)
    size_factor = 3.08567758149137e13 * 1e-13 * 0.15
    return np.power(pseudo_l, 0.5) * size_factor

def convert_teff_to_rgb(t_eff):
    tuple_shape = (len(t_eff),)
    dtype = [('r', 'f8'), ('g', 'f8'), ('b', 'f8')]
    stellar_color = np.zeros(tuple_shape, dtype=dtype)
    for i in range(len(t_eff)):
        stellar_color[i] = convert_K_to_RGB(t_eff[i].value)
    return stellar_color

def angular_distance_between_stars_and_exo_pl(ra1, dec1, ra2, dec2):
    # https://en.wikipedia.org/wiki/Angular_distance
    # https://pyastronomy.readthedocs.io/en/latest/pyaslDoc/aslDoc/angularDistance.html

    # theta = pyasl.getAngDist(ra1, dec1, ra2, dec2)

    theta = np.zeros(len(ra1))
    for i in range(len(ra1)):
        theta[i] = pyasl.getAngDist(ra1[i], dec1[i], ra2, dec2)

    return theta

def possition_angle_of_exo_pl_relative_to_stars(pl_ra, pl_dc, star_ra, star_dc):
    # https://pyastronomy.readthedocs.io/en/latest/pyaslDoc/aslDoc/posAngle.html
    # Position angle of stars from Exo-planet (pl_ra, pl_dc)
    # Calculate position angle
    r = np.zeros(len(star_ra))
    for i in range(len(star_ra)):
        r[i] = pyasl.positionAngle(pl_ra, pl_dc, star_ra[i], star_dc[i])
    return r

def distance_to_stars_from_exo_pl(theta, gaia_dist, pl_dist):
    # https://pypi.org/project/trianglesolver/

    tuple_shape = (len(gaia_dist),)
    # dtype = [('a', 'f8'), ('b', 'f8'), ('c', 'f8'), ('A', 'f8'), ('B', 'f8'), ('C', 'f8')]
    dtype = [('a', 'f8'), ('B', 'f8'), ('C', 'f8'), ('A_C', 'f8')]
    star_dist = np.zeros(tuple_shape, dtype=dtype)
    for i in range(len(gaia_dist)):
        a,b,c,A,B,C = solve(b=gaia_dist[i], c=pl_dist, A=theta[i]*degree)
        star_dist[i] = a, B / degree, C / degree, theta[i] + (C / degree)

    # print(star_dist)
    return star_dist