# Night-ExoSky

### Requirements

- Hard drive: ~500 MB of free disk space (including downloaded datasets)
- RAM: 2-4 GB
- CPU: Intel Core i5 3rd Gen. 4+ cores recommended
- Operating system: Linux / Windows 7+ / macOS

### Installation:

```shell script
cd night-exosky
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt 
```

### Preparation of datasets

```shell script
python download_data.py
```
- Note: download time ~5 minutes

### Run

```shell script
python main.py
```
- Note: loading time of the Gaia "Light" dataset can be about 1 minute 
- Warning: loading time of the Gaia "Big" dataset can be about 10 minutes 

### References

- https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html
- https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=PS
- https://www.cosmos.esa.int/web/gaia-users/archive/extract-data#bulk_download
- https://www.cosmos.esa.int/web/gaia-users/archive/datalink-products#datalink_jntb_get_above_lim
- https://gaia.ari.uni-heidelberg.de/gaiasky/docs/master/STIL-loader.html#magnitudes
- https://gaia.ari.uni-heidelberg.de/gaiasky/docs/master/STIL-loader.html#colors
- https://en.wikipedia.org/wiki/Color_index
- http://www.vendian.org/mncharity/dir3/blackbody/
- https://tannerhelland.com/2012/09/18/convert-temperature-rgb-algorithm-code.html
- https://gist.github.com/petrklus/b1f427accdf7438606a6
- https://gaia.ari.uni-heidelberg.de/gaiasky/docs/master/LOD-catalogs.html#catalogs
- https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
- https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html#PS
- https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html
- https://astroquery.readthedocs.io/en/latest/#astroquery
- https://pyvo.readthedocs.io/en/latest/#pyvo
- https://pyvo.readthedocs.io/en/latest/dal/index.html#pyvo-data-access
- https://en.wikipedia.org/wiki/Angular_distance
- https://pyastronomy.readthedocs.io/en/latest/pyaslDoc/aslDoc/angularDistance.html
- https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
- https://science.nasa.gov/learn/basics-of-space-flight/chapter2-2/
- https://stellarium.org/
- https://www.ivoa.net/documents/VOTable/20130920/REC-VOTable-1.3-20130920.html
- https://fits.gsfc.nasa.gov/fits_wcs.html
- http://tdc-www.harvard.edu/software/wcstools/wcstools.wcs.html
- https://learn.astropy.org/tutorials/1-Coordinates-Intro.html
- https://learn.astropy.org/tutorials/2-Coordinates-Transforms.html
- https://learn.astropy.org/tutorials/3-Coordinates-Velocities.html
- https://learn.astropy.org/tutorials/gaia-galactic-orbits.html
- https://learn.astropy.org/tutorials/FITS-cubes.html
- https://viyaleta.medium.com/how-to-make-a-sky-map-in-python-a362bf722bb2
- https://iopscience.iop.org/article/10.3847/1538-3881/aabc4f/pdf
- https://academic.oup.com/mnras/article/529/4/3816/7627451
- https://learn.astropy.org/tutorials/celestial_coords1.html
- https://blog.finxter.com/5-best-ways-to-create-numpy-arrays-of-tuples-in-python/
- https://nsweb.tn.tudelft.nl/~gsteele/TN2513_Tips/Creating%20a%20numpy%20array%20in%20a%20loop.html
- https://matplotlib.org/stable/users/explain/colors/colors.html
- https://www.pythonguis.com/tutorials/pyside-plotting-matplotlib/
- https://matplotlib.org/stable/gallery/widgets/polygon_selector_demo.html
- https://docs.astropy.org/en/stable/table/index.html
- https://decovar.dev/blog/2022/02/26/astronomy-databases-tap-adql/
- https://pypi.org/project/trianglesolver/
- https://docs.astropy.org/en/stable/cosmology/io/builtin.html#module-astropy.cosmology._io.table
- https://en.wikipedia.org/wiki/Sagittarius_A*
- https://docs.astropy.org/en/stable/coordinates/skycoord.html
