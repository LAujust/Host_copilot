import numpy as np
import sys, os
from astropy.table import Table
from astropy.io import fits
import requests
from PIL import Image
from io import BytesIO
import matplotlib.pyplot as plt



class Imager:
    def __init__(self,ra,dec,r_arcsec,band='r',save_path=None):
        self.ra = ra
        self.dec = dec
        self.r_arcsec = r_arcsec
        self.band = band
        self.fov = (self.r_arcsec/3600)**2
        self.save_path = save_path
        
        if self.save_path:
            os.makedirs(self.save_path, exist_ok=True)
            
            
    def get_cutout(self,):
        
        #try PS first
        try:
            fits_dir = self.PS_cutout()
        except Exception as e:
            print(f"Pan-STARRS cutout failed: {e}")
            fits_dir = None
            
        if not fits_dir:
            #Fall-back to LS
            try:
                fits_dir = self.LS_cutout()
            except Exception as e:
                print(f"Legacy Survey cutout failed: {e}")
                fits_dir = None
        
        if not fits_dir:
            print("Both Pan-STARRS and Legacy Survey cutouts failed.")
            return None
        else:
            return fits_dir
        
        
        
        
    def PS_cutout(self):
        """Get cutout from Pan-STARRS DR2"""
        
        def _getimages(ra,dec,band):
            
            """Query ps1filenames.py service to get a list of images
            
            ra, dec = position in degrees
            size = image size in pixels (0.25 arcsec/pixel)
            filters = string with filters to include
            Returns a table with the results
            """
            
            service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
            url = f"{service}?ra={ra}&dec={dec}&filters={band}"
            table = Table.read(url, format='ascii')
            return table
        
        
        def _geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
    
            """Get URL for images in the table
            
            ra, dec = position in degrees
            size = extracted image size in pixels (0.25 arcsec/pixel)
            output_size = output (display) image size in pixels (default = size).
                        output_size has no effect for fits format images.
            filters = string with filters to include
            format = data format (options are "jpg", "png" or "fits")
            color = if True, creates a color image (only for jpg or png format).
                    Default is return a list of URLs for single-filter grayscale images.
            Returns a string with the URL
            """
            
            if color and format == "fits":
                raise ValueError("color images are available only for jpg or png formats")
            if format not in ("jpg","png","fits"):
                raise ValueError("format must be one of jpg, png, fits")
            table = _getimages(ra,dec,filters=filters)
            url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
                f"ra={ra}&dec={dec}&size={size}&format={format}")
            if output_size:
                url = url + "&output_size={}".format(output_size)
            # sort filters from red to blue
            flist = ["yzirg".find(x) for x in table['filter']]
            table = table[np.argsort(flist)]
            if color:
                if len(table) > 3:
                    # pick 3 filters
                    table = table[[0,len(table)//2,len(table)-1]]
                for i, param in enumerate(["red","green","blue"]):
                    url = url + "&{}={}".format(param,table['filename'][i])
            else:
                urlbase = url + "&red="
                url = []
                for filename in table['filename']:
                    url.append(urlbase+filename)
            return url
        
        pix_size = int(2*self.r_arcsec/0.25)
        fitsurl = _geturl(self.ra, self.dec, size=pix_size, filters=self.band, format="fits")
        fh = fits.open(fitsurl[0])
        fh.writeto(os.path.join(self.save_path,'ps1_%s_ref.fits'%self.band),overwrite=True)
        print("Saved PS1 cutout to %s"%(os.path.join(self.save_path,'ps1_%s_ref.fits'%self.band)))
        return os.path.join(self.save_path,'ps1_%s_ref.fits'%self.band)
    
    
    def LS_cutout(self):
        from pyvo.dal import sia
        from astropy.utils.data import download_file
        # Connect to the Legacy Survey DR10 SIA service
        FOV = self.r_arcsec/3600 #to degree
        
        DEF_ACCESS_URL = "https://datalab.noirlab.edu/sia/ls_dr10"
        svc = sia.SIAService(DEF_ACCESS_URL)

        # Search for images overlapping the position
        # RA FOV is divided by cos(dec) to correct for spherical projection
        imgTable = svc.search(
            (self.ra, self.dec),
            (FOV / np.cos(self.dec * np.pi / 180), FOV),
            verbosity=2
        ).to_table()

        print(f"Found {len(imgTable)} images overlapping position")
        print(imgTable['obs_bandpass', 'proctype', 'prodtype'])


            # Filter for stacked images in this band
        sel = (
                (imgTable['proctype'].astype(str) == 'Stack') &
                (imgTable['prodtype'].astype(str) == 'image') &
                (np.char.startswith(imgTable['obs_bandpass'].astype(str), self.band))
            )

        if not np.any(sel):
            print(f"No {self.band}-band stacked image found, skipping.")
            return None

        row = imgTable[sel][0]
        url = row['access_url']
        print(f"Downloading {self.band}-band cutout: {url}")

        filename = download_file(url, cache=True, show_progress=True, timeout=120)
        hdu = fits.open(filename)
        hdu.writeto(os.path.join(self.save_path, f"ls_{self.band}.fits"), overwrite=True)
        hdu.close()

        print(f"Saved {self.band}-band cutout to {self.save_path}/ls_{self.band}.fits")
        
        return os.path.join(self.save_path, f"ls_{self.band}.fits")