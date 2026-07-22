from astroquery.vizier import Vizier
import numpy as np
import requests
from astropy.io import votable
from io import BytesIO, StringIO
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd


class GalaxyFinder:
    def __init__(self,ra:float,dec:float,r_arcsec:float,save_path='./'):
        self.ra = self.ra
        self.dec = dec
        self.r_arcsec = r_arcsec
        self.save_path = save_path
        self.r_search = self.r_arcsec + 10 #galaxy morphology


    #===================================
    #         Source Catalog      
    #===================================
    def find_ps(self):
        """
        Cone search for Pan-STARRS catalog
        """
        params = {
            "ra": self.ra,
            "dec": self.dec,
            "radius": self.r_search/3600,  #to degree
            "nDetections.gt": 4,
        }

        response = requests.get(self.BASE_URL, params=params)
        if response.status_code != 200:
            raise RuntimeError(
                f"Pan-STARRS query failed (status {response.status_code})"
            )

        df = pd.read_csv(StringIO(response.text))

        # ---- rename columns ----
        new_columns = {}
        for col in df.columns:
            if col == "raMean":
                new_col = "RA"
            elif col == "decMean":
                new_col = "DEC"
            elif col.endswith("MeanPSFMagErr"):
                new_col = col.replace("MeanPSFMagErr", "_err")
            elif col.endswith("MeanPSFMag"):
                new_col = col.replace("MeanPSFMag", "")
            else:
                new_col = col
            new_columns[col] = new_col

        df = df.rename(columns=new_columns)
        df = df.replace(-999, np.nan)

        self.ps_df = df

    def find_ls(self):
        """
        TAP Search for Legacy Survey Catalog
        """
        import pyvo
        # Connect to the NOIRLab TAP service
        tap = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")

        # Use bounding box instead of cone search
        ra_center, dec_center = self.ra, self.dec
        radius = self.r_search / 3600  # convert arcsec to degrees

        query = f"""
        SELECT TOP 4000
        ra, dec, ls_id, flux_g, flux_ivar_g,
        flux_i, flux_ivar_i, flux_r, flux_ivar_r,
        flux_z, flux_ivar_z
        FROM ls_dr10.tractor
        WHERE ra BETWEEN {ra_center - radius} AND {ra_center + radius}
        AND dec BETWEEN {dec_center - radius} AND {dec_center + radius}
        """

        try:
            # Run query
            results = tap.search(query)
            df = results.to_table().to_pandas()

            # Filter out non-positive flux or ivar (to avoid log errors)
            bands = ['g', 'r', 'i', 'z']
            for band in bands:
                flux = f'flux_{band}'
                ivar = f'flux_ivar_{band}'
                mask = (df[flux] > 0) & (df[ivar] > 0)
                df = df[mask].copy()

                # Convert flux to magnitude and calculate errors
                df[f'{band}'] = 22.5 - 2.5 * np.log10(df[flux])
                df[f'{band}_err'] = 2.5 / np.log(10) * (1 / (np.sqrt(df[ivar]) * df[flux]))

                self.ls_df = df
        except Exception as e:
            print(f"Cannot find LS catalog: {e}")


    #===================================
    #     Galaxy Specified Catalog      
    #===================================

    def find_reglade(self,max_rows=20):
        """
        Cone search REGALADE catalog through VizieR.

        Parameters
        ----------
        max_rows : int
            Maximum number of returned rows

        Returns
        -------
        pandas.DataFrame
        """

        coord = SkyCoord(
            ra=self.ra,
            dec=self.dec,
            unit="deg",
            frame="icrs"
        )


        viz = Vizier(
            columns=["*"],
            row_limit=max_rows
        )


        result = viz.query_region(
            coord,
            radius=self.r_search*u.arcsec,
            catalog="J/A+A/706/A284/regalade"
        )


        if len(result) == 0:
            return pd.DataFrame()


        df = result[0].to_pandas()

        self.reglade_df = df
    



    def find_ned(self,maxrec=20):
        """
        Query NED ConeSearchByPosition API

        Parameters
        ----------
        maxrec : int
            Maximum number of returned objects
        """

        url = (
            "https://ned.ipac.caltech.edu/NED::API/"
            "ConeSearchByPosition"
        )

        params = {
            "RA": f"{self.ra}d",
            "DEC": f"{self.dec}",
            "CSYS": "Equatorial",
            "EQUINOX": "J2000.0",
            "SR": f"{self.r_search/60}",  #convert into arcmin
            "MAXREC": f"{maxrec}"
        }


        response = requests.get(url, params=params, timeout=120)


        if response.status_code != 200:
            print(response.text)
            response.raise_for_status()


        # NED returns VOTable
        vot = votable.parse(
            BytesIO(response.content)
        )

        table = vot.get_first_table().to_table()
        self.ned_df = table.to_pandas()

    
