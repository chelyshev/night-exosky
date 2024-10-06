import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from kelvin_to_rgb import convert_K_to_RGB

from astropy import units as u
from astropy.coordinates import (SkyCoord, Distance, Galactic,
                                 EarthLocation, AltAz)
from astropy.table import QTable

import ipywidgets as widgets
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, clear_output

def apparent_to_absolute_magnitude(dist_pc, app_mag):
    # dist = 10.0 if dist_pc <= 0.0 else dist_pc
    v = 5.0 * np.log10(dist_pc)
    return app_mag - v + 5.0

def absolute_magnitude_to_pseudosize(abs_mag):
    # Pseudo-luminosity. Usually L = L0 * 10^(-0.4*Mbol). We omit M0 and approximate Mbol = M
    pseudo_l = np.power(10, -0.4 * abs_mag)
    size_factor = 3.08567758149137e13 * 1e-13 * 0.15
    return np.power(pseudo_l, 0.5) * size_factor

def build_data_visualization(data):
    tbl = QTable.read(data, format='fits')

    print(len(tbl))
    # only keep bright stars (G magnitude < 13.1)
    # tbl = tbl[tbl['phot_g_mean_mag'] < 13.1*u.mag]
    # print(len(tbl))
    gaia_dist = Distance(parallax=tbl['parallax'])

    t_eff = tbl['teff_gspphot']
    print(t_eff)
    tuple_shape = (len(t_eff),)
    dtype = [('r', 'f8'), ('g', 'f8'), ('b', 'f8')]
    stelar_color = np.zeros(tuple_shape, dtype=dtype)
    for i in range(len(t_eff)):
        # print(t_eff[i])
        stelar_color[i] = convert_K_to_RGB(t_eff[i].value)

    print(stelar_color)

    sky_coords_icrs = SkyCoord(
        ra=tbl['ra'],
        dec=tbl['dec'],
        distance=Distance(parallax=tbl['parallax']),
        frame='icrs')


    scg = sky_coords_icrs.galactic
    magnitude = tbl['phot_g_mean_mag']
    parralax = tbl['parallax']
    # print(parralax)
    # print(gaia_dist)
    # print(1000 / tbl['parallax'].value)
    print(magnitude.value)
    abs_mag = apparent_to_absolute_magnitude(scg.distance.value, magnitude.value)
    print(abs_mag)
    pseudosize = absolute_magnitude_to_pseudosize(abs_mag)

    print(pseudosize)

    # fig, ax = plt.subplots(figsize=(18, 9), constrained_layout=True)
    # fig.patch.set_facecolor('black')
    # border = plt.Rectangle((0, 0), 1, 1, color='navy', fill = True)
    # ax.add_patch(border)

    # ax.patch.set_facecolor('black')

    # marker_size = max_star_size * 10 ** (magnitude / -2.5)
    # im1 = ax.scatter(sky_coords_3d.ra.degree, sky_coords_3d.dec.degree, s=(max_star_size * 10 ** (magnitude.value / -2.5)), color='white')

    # im1 = ax.scatter(-scg.l.wrap_at(180 * u.deg).degree, scg.b.degree, s=pseudosize, c='white')

    # im = ax.imshow(scg, cmap='viridis',  norm=LogNorm(vmin=2, vmax=50))
    # im1 = ax.scatter(scg.l.degree, scg.b.degree, s=(max_star_size * 10 ** (magnitude.value / -2.5)), c=scg.distance.kpc, cmap='twilight')
    # ax.set_xlabel('RA [deg]')
    # ax.set_ylabel('Dec [deg]')

    # Add axes labels
    # ax.set_xlabel("Right Ascension", fontsize=16)
    # ax.set_ylabel("Declination", fontsize=16)
    fig, ax = plt.subplots(figsize=(18, 9), constrained_layout=True)
    ax.patch.set_facecolor('black')

    # sph = scg.spherical
    # im1 = ax.scatter(-sph.lon.wrap_at(180*u.deg).radian, sph.lat.radian, s=pseudosize, c='white')
    im1 = ax.scatter(-scg.l.wrap_at(180 * u.deg).degree, scg.b.degree, s=pseudosize, color=stelar_color, marker='.', linewidths=0)

    # def fmt_func(x, pos):
    #     val = -sph.lat.wrap_at(360*u.deg).degree
    #     return f'${val}'
    #
    # ticker = mpl.ticker.FuncFormatter(fmt_func)
    # ax.xaxis.set_major_formatter(ticker)

    # ax.grid()

    ax.grid(ls='solid', lw=.1)

    # Overlay set of Galactic Coordinate Axes
    # overlay = ax.get_coords_overlay('galactic')
    # overlay.grid(color='black', ls='dotted', lw=1)
    # overlay[0].set_axislabel('Galactic Longitude', fontsize=14)
    # overlay[1].set_axislabel('Galactic Latitude', fontsize=14)

    plt.xlabel('Galactic Longitude')
    plt.ylabel('Galactic Latitude')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    # plt.axis('off')
    plt.show()



if __name__ == '__main__':
    build_data_visualization("gaia_data-t6.fits")

