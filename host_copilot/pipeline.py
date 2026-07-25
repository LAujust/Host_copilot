from .catalog import GalaxyFinder
from .image import Imager
from .utils import *
import astropy.units as u
from ipyaladin import Aladin, EllipseError
from astropy.coordinates import SkyCoord, Angle
from regions import (
        EllipseSkyRegion,
        RegionVisual,
        CircleSkyRegion
    )



class HostPipeline:
    def __init__(self, ra, dec, r_arcsec, zcutout=0.1, quick=True, save_path='./'):
        self.ra = ra
        self.dec = dec
        self.r_arcsec = r_arcsec
        self.zcutout = zcutout
        self.save_path = save_path
        self.quick = quick
        self.pos = SkyCoord(ra, dec, unit='deg', frame='icrs')

        # Initialize GalaxyFinder and Imager
        self.galaxy_finder = GalaxyFinder(ra, dec, r_arcsec, save_path=save_path)
        self.imager = Imager(ra, dec, r_arcsec, save_path=save_path)
        
        
    def filter_and_visualize(self):
        """
        Visualize the results using matplotlib
        """

        fov = (3*self.r_arcsec/3600)**2

        if self.quick:
            if len(self.galaxy_finder.reglade_df)>0 or self.galaxy_finder.reglade_df is not None:
                cat_table = Table.from_pandas(self.galaxy_finder.reglade_df)
                #Only REGLADE catalog
                print(f"")
                cat_table['R1'] = cat_table['R1']*u.arcsec
                cat_table['R2'] = cat_table['R2']*u.arcsec
                cat_table['PA'] = cat_table['PA']*u.deg
                cat_table = cat_table[cat_table['z']<self.zcutout]
                cat_table['sep'] = self.pos.separation(SkyCoord(cat_table['RAJ2000'], cat_table['DEJ2000'], unit='deg')).arcsec * u.arcsec
                
                if len(cat_table) == 0:
                    print("No galaxies found in the REGLADE catalog within the specified redshift cut.")
                    return None, None
                
                else:
                    
                    print(f"Found {len(cat_table)} galaxies in the REGLADE catalog within the specified redshift cut (z<{self.zcutout}) and radius r={self.r_arcsec} arcsec.")
                    print('-'*50)
                    for row in cat_table:
                        print(f"Galaxy: {row['Name']}, z={row['z']:.3f}, sep={row['sep']:.2f}")
                
                    aladin = Aladin(fov=fov, target=self.pos ,survey='CDS/P/PanSTARRS/DR1/color-z-zg-g')
                    
                    #Add Catalog
                    aladin.add_table(cat_table,
                                    shape=EllipseError(
                                    maj_axis="R1",
                                    min_axis="R2",
                                    angle="PA",
                                    default_shape="cross"),
                                    color="cyan",
                                    )


                    #Add Source
                    circle = CircleSkyRegion(
                        center=self.pos, radius=Angle(self.r_arcsec, "arcsec"), visual={"edgecolor": "yellow",'linestyle':'dashed'}
                    )
                    aladin.add_graphic_overlay_from_region([circle])
                    
                    return aladin, cat_table
            else:
                print("No galaxies found in the REGLADE catalog.")
                return None, None
        
        
        
        
    

    def run(self):
        print('=' * 50)
        
        if self.quick:
            print("[QUICK MODE]")
            # Step 1: Find galaxies in the vicinity
            print("Searching for galaxies...")
            self.galaxy_finder.run(quick=self.quick)
            aladin, cat_table = self.filter_and_visualize()
            print('HostPipeline run completed.')
            print('=' * 50)
            return aladin, cat_table

            # # Step 2: Get cutout images
            # print("Retrieving cutout images...")
            # fits_dir = self.imager.get_cutout()
            # if fits_dir:
            #     print(f"Cutout images saved to {fits_dir}")
            # else:
            #     print("Failed to retrieve cutout images.")
            
            

                
            
    
    