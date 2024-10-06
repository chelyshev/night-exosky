from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default

import pyvo as vo
from astropy.io.votable import from_table, writeto
from astropy import units as u
import numpy as np

from timebudget import timebudget
timebudget.set_quiet()  # don't show measurements as they happen
timebudget.report_at_exit()  # Generate report when the program exits

@timebudget
def download_nasa_exoplanet_archive_data():

    service = vo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP")

    # https://decovar.dev/blog/2022/02/26/astronomy-databases-tap-adql/
    x_query = """
        SELECT * FROM
        (WITH latestEntries AS
            (SELECT pl_rade, pl_radj, pl_pubdate, pl_name, hostname, ra, dec, sy_plx, sy_dist, sy_gaiamag, st_teff,
            ROW_NUMBER() OVER(
                PARTITION BY pl_name ORDER BY CASE WHEN pl_radj IS NULL THEN 1 ELSE 0 END, pl_pubdate DESC
            ) AS rank
            FROM ps)
        SELECT pl_name, hostname, ra, dec, sy_plx, sy_dist, pl_radj, pl_rade, pl_pubdate, sy_gaiamag, st_teff FROM latestEntries 
        WHERE rank = 1 AND sy_plx is not Null AND sy_gaiamag is not Null AND st_teff is not Null ORDER BY sy_dist )
        """

    resultset = service.search(x_query).to_table()

    votable = from_table(resultset)
    writeto(votable, "nasa_ps_data.vot")
    print("Filtered and downloaded ", len(resultset), "rows from NASA Exoplanet Archive TAP resource")

@timebudget
def download_gaia_data():
    query = "SELECT TOP 1000000 source_id, ra, dec, parallax, phot_g_mean_mag, teff_gspphot from gaiadr3.gaia_source \
    WHERE phot_g_mean_mag is not Null AND \
    teff_gspphot is not Null AND \
    parallax is not Null AND  \
    parallax_error is not Null AND \
    parallax_error < parallax * 0.2 AND \
    phot_g_mean_mag < 19"

    job = Gaia.launch_job_async(query)
    gaia_data = job.get_results()

    gaia_data.write('gaia_data.fits', format='fits', overwrite=True)
    # gaia_data_write.tofile('gaia_data-t2.fits')
    print("Filtered and downloaded: ", len(gaia_data), "rows from ESA Gaia resource")

if __name__ == '__main__':
    download_gaia_data()
    download_nasa_exoplanet_archive_data()
