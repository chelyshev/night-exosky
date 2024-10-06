import sys
import matplotlib

from calculations import apparent_to_absolute_magnitude, absolute_magnitude_to_pseudosize, convert_teff_to_rgb, \
    angular_distance_between_stars_and_exo_pl, distance_to_stars_from_exo_pl, \
    possition_angle_of_exo_pl_relative_to_stars, apparent_to_absolute_magnitude_t

matplotlib.use('Qt5Agg')

from PySide2.QtWidgets import QMainWindow, QVBoxLayout, QWidget, QApplication, QLabel, QComboBox, QPushButton

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PySide2.QtCore import Qt, QSize
from astropy.table import QTable
from kelvin_to_rgb import convert_K_to_RGB

from astropy import units as u
from astropy.coordinates import (SkyCoord, Distance, Galactic,
                                 EarthLocation, Galactocentric)
import astropy.coordinates as coord

import matplotlib.pyplot as plt

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=18, height=9, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi, constrained_layout=True)
        fig.supxlabel('Galactic Longitude', fontsize=10)
        fig.supylabel('Galactic Latitude', fontsize=10)
        self.axes = fig.add_subplot(111, facecolor="black")
        self.axes.patch.set_facecolor('black')
        super(MplCanvas, self).__init__(fig)


class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.setWindowTitle("Night Exo-Sky app")

        self.gaia_tbl = QTable.read('gaia_data.fits', format='fits')
        self.nasa_ps_tbl = QTable.read('nasa_ps_data.vot', format='votable')

        sky_coords_icrs = SkyCoord(
            ra=self.gaia_tbl['ra'],
            dec=self.gaia_tbl['dec'],
            distance=Distance(parallax=self.gaia_tbl['parallax']),
            frame='icrs')
        self.sky_coords_galactic = sky_coords_icrs.galactic
        self.night_exosky_coords_galactic = self.sky_coords_galactic

        pl_coords_icrs = SkyCoord(
            ra=self.nasa_ps_tbl['ra'],
            dec=self.nasa_ps_tbl['dec'],
            distance=self.nasa_ps_tbl['sy_dist'],
            frame='icrs')
        self.pl_coords_galactic = pl_coords_icrs.galactic

        magnitude = self.gaia_tbl['phot_g_mean_mag']
        # print(magnitude.value)
        abs_mag = apparent_to_absolute_magnitude(self.sky_coords_galactic.distance.value, magnitude.value)
        # print(abs_mag)
        self.pseudosize = absolute_magnitude_to_pseudosize(abs_mag)
        # print(self.pseudosize)

        t_eff = self.gaia_tbl['teff_gspphot']
        # print(t_eff)
        self.stellar_color = convert_teff_to_rgb(t_eff)
        # print(self.stellar_color)

        # pl_magnitude = self.nasa_ps_tbl['sy_gaiamag']
        # print(pl_magnitude.value)
        # pl_abs_mag = apparent_to_absolute_magnitude(self.sky_coords_galactic.distance.value, pl_magnitude.value)
        # print(pl_abs_mag)
        # pl_pseudosize = absolute_magnitude_to_pseudosize(pl_abs_mag)
        # print(pl_pseudosize)
        #
        # pl_t_eff = self.nasa_ps_tbl['st_teff']
        # print(pl_t_eff)
        # pl_stellar_color = convert_teff_to_rgb(pl_t_eff)
        # print(pl_stellar_color)

        self.sc = MplCanvas(self, width=18, height=9, dpi=100)
        # self.axes.scatter(-scg.l.wrap_at(180 * u.deg).degree, scg.b.degree, s=pseudosize, color=stellar_color, marker='.',
        #                  linewidths=0)
        self.sc.axes.scatter(self.sky_coords_galactic.l.degree, self.sky_coords_galactic.b.degree, s=self.pseudosize, color=self.stellar_color,
                          marker='.', linewidths=0)
        self.sc.axes.grid(ls='solid', lw=.1)
        self.sc.axes.set_xlim(0, 360)
        self.sc.axes.set_ylim(-90, 90)
        self.sc.axes.set_title(label='Default: Galactic View from Earth')

        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar = NavigationToolbar(self.sc, self)

        dropdown_list = QComboBox()
        pl_name = self.nasa_ps_tbl['pl_name']
        dropdown_list.addItem("Earth")
        dropdown_list.addItems(pl_name)

        # The default signal from currentIndexChanged sends the index
        dropdown_list.currentIndexChanged.connect(self.index_changed)

        # The same signal can send a text string
        # dropdown_list.currentTextChanged.connect(self.text_changed)

        self.setCentralWidget(dropdown_list)

        button = QPushButton("Construct visualization")

        # self.setFixedSize(QSize(400, 300))

        # Set the central widget of the Window.
        # self.setCentralWidget(button)

        layout = QVBoxLayout()
        layout.addWidget(toolbar)
        layout.addWidget(self.sc)
        layout.addWidget(dropdown_list)
        layout.addWidget(button)


        # Create a placeholder widget to hold our toolbar and canvas.
        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)
        # self.updatesEnabled()
        self.show() #showFullScreen()


    def index_changed(self, i):
        i = i - 1

        if i == -1:
            print("Selected Earth to construct night sky visualization")
            label = 'Galactic View from Earth'
            fig, ax = plt.subplots(figsize=(18, 9), constrained_layout=True)
            ax.patch.set_facecolor('black')
            plt.title(label=label)
            ax.scatter(self.sky_coords_galactic.l.degree, self.sky_coords_galactic.b.degree,
                       s=self.pseudosize, color=self.stellar_color, marker='.',
                       linewidths=0, edgecolors='face')
            plt.xlabel('Galactic Longitude')
            plt.ylabel('Galactic Latitude')
            ax.set_xlim(0, 360)
            ax.set_ylim(-90, 90)
            ax.grid(ls='solid', lw=.1)
            plt.show()

        else:
            pl_name = self.nasa_ps_tbl['pl_name']
            print("Selected the following exo-planet to construct night sky visualization")
            print(self.nasa_ps_tbl[i])
            label = 'Galactic View from "' + str(pl_name[i] + '"')

            ra = self.gaia_tbl['ra'].value
            dec = self.gaia_tbl['dec'].value
            pl_ra = self.nasa_ps_tbl['ra'].value
            pl_dec = self.nasa_ps_tbl['dec'].value
            theta = angular_distance_between_stars_and_exo_pl(ra, dec, pl_ra[i], pl_dec[i])
            print("theta:", theta)

            l = self.sky_coords_galactic.l.degree
            b = self.sky_coords_galactic.b.degree
            # print(l, b)
            pl_l = self.pl_coords_galactic.l.degree
            pl_b = self.pl_coords_galactic.b.degree
            # print(pl_l[i], pl_b[i])
            theta = angular_distance_between_stars_and_exo_pl(l, b, pl_l[i], pl_b[i])
            print("theta:", theta)

            # pos_angle = possition_angle_of_exo_pl_relative_to_stars(pl_ra[i], pl_dec[i], ra, dec)
            # print("pos_angle:", pos_angle)

            gaia_dist = Distance(parallax=self.gaia_tbl['parallax']).value
            pl_dist = self.pl_coords_galactic.distance.value
            stars_dist = distance_to_stars_from_exo_pl(theta, gaia_dist, pl_dist[i])
            print(theta[0], stars_dist['B'][0], stars_dist['C'][0], stars_dist['A_C'][0], ra[0], dec[0], gaia_dist[0], pl_dist[i], stars_dist['a'][0])

            # self.night_exosky_coords_galactic = SkyCoord(
            #     l=stars_dist['l']*u.deg,
            #     b=b*u.deg,
            #     distance=stars_dist['a']*u.pc,
            #     frame='galactic'
            #     )
            # night_exosky_coords_icrs = SkyCoord(
            #     ra=ra * u.deg, # stars_dist['A_C'] * u.deg,
            #     dec=dec * u.deg,
            #     distance=Distance(parallax=self.gaia_tbl['parallax']), # stars_dist['a'] * u.pc,
            #     frame='icrs'
            # )
            night_exosky_coords_icrs = SkyCoord(
                ra=stars_dist['A_C'] * u.deg,
                dec=dec * u.deg,
                distance=stars_dist['a'] * u.pc,
                frame='icrs'
            )

            print(night_exosky_coords_icrs)

            magnitude = self.gaia_tbl['phot_g_mean_mag']
            # print(magnitude.value)
            abs_mag = apparent_to_absolute_magnitude(night_exosky_coords_icrs.distance.value, magnitude.value)
            # print('abs_mag:', abs_mag)
            pseudosize = absolute_magnitude_to_pseudosize(abs_mag)
            print(pseudosize)
            print(self.pseudosize)

            print(self.stellar_color)

            self.night_exosky_coords_galactic = night_exosky_coords_icrs.transform_to(Galactic())

            # self.sc.axes.set_title(label=label)
            # self.sc.axes.scatter(-self.night_exosky_coords_galactic.l.degree, self.night_exosky_coords_galactic.b.degree,
            #                      s=pseudosize, color=self.stellar_color,
            #                      marker='.', linewidths=0)
            # self.sc.axes.grid(ls='solid', lw=.2)
            # self.update()
            # self.show()
            # self.showFullScreen()


            fig, ax = plt.subplots(figsize=(18, 9), constrained_layout=True)
            ax.patch.set_facecolor('black')
            plt.title(label=label)


            # Galactic coords in Longitude(l), Latitude(b) degrees
            # print(self.night_exosky_coords_galactic)
            ax.scatter(self.night_exosky_coords_galactic.l.degree, self.night_exosky_coords_galactic.b.degree, s=pseudosize, color=self.stellar_color, marker='.',
                             linewidths=0, edgecolors='face')
            plt.xlabel('Galactic Longitude')
            plt.ylabel('Galactic Latitude')
            # ax.set_xlim(0, 360)
            # ax.set_ylim(-90, 90)


            # ICRS coords in RA, DEC degrees
            # im1 = ax.scatter(night_exosky_coords_icrs.ra.value,  night_exosky_coords_icrs.dec.value, s=self.pseudosize, color=self.stellar_color, marker='.',
            #                  linewidths=0)
            # plt.xlabel('RA')
            # plt.ylabel('DEC')
            # ax.set_xlim(0, 360)
            # ax.set_ylim(-90, 90)
            # plt.axis('off')


            # Plot x,y,x coords in pc from Galactic center
            # self.night_exosky_coords_galactic.transform_to(Galactocentric())
            # print(self.night_exosky_coords_galactic)
            # m1 = ax.scatter(self.night_exosky_coords_galactic.x, self.night_exosky_coords_galactic.y,
            #                 s=self.pseudosize, color=self.stellar_color, marker='.',
            #                 linewidths=0)
            # plt.xlabel('pc')
            # plt.ylabel('pc')

            ax.grid(ls='solid', lw=.1)
            plt.show()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MainWindow()
    app.exec_()