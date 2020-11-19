import os, sys, importlib
import rasterio, affine

import geopandas as gpd

from rasterio import features
from affine import Affine
from shapely.geometry import box

class gsm_rasterizer(object):
    def __init__(self, gsm_files, out_folder):
        ''' Used to open, process, and extract GSM coverage files on a country by country basis
        
        :param gsm_files: list of string paths defining the GSM coverage data
        '''
        raw_gsm = {}
        for f in gsm_files:
            gsm_name = os.path.basename(f).replace(".shp", "")
            raw_gsm[gsm_name] = f
        self.gsm_files = raw_gsm
        self.out_folder = out_folder
        if not os.path.exists(self.out_folder):
            os.makedirs(out_folder)
        
    def initial_read_in(self, gsm_files=[]):
        if len(gsm_files) == 0:
            gsm_files = self.gsm_files
        gsm_data = {}
        
        for key, value in gsm_files.items():
            gsm_data[key] = gpd.read_file(value)
            
        self.gsm_data = gsm_data
    
    def get_exact_shape(self, x, boundary):
        if boundary.contains(x):
            return(x)
        else:
            try:
                return(boundary.intersection(x))
            except:
                x = x.buffer(0)
                return(boundary.intersection(x))

            
    def extract_country_vectors(self, iso3, global_data, out_folder = ""):    
        country_data = global_data[global_data['ISO3'] == iso3]
        sample_key = list(self.gsm_data.keys())[0]
        if country_data.crs != self.gsm_data[sample_key].crs:
            country_data = country_data.to_crs(self.gsm_data[sample_key].crs)
        if out_folder == "":
            out_folder = self.out_folder
        self.country_data = country_data
        
        c_shp = box(*country_data.total_bounds)
        self.c_shp = c_shp
        out_gsm_files = []
        
        for key, gsm in self.gsm_data.items():
            out_file = os.path.join(out_folder, "%s.shp" % key)
            out_gsm_files.append(out_file)
            if not os.path.exists(out_file):
                select = gsm[gsm.intersects(c_shp)]    
                if select.shape[0] > 0:
                    # clip geometries for features to actual boundary
                    select['geometry'] = select['geometry'].apply(lambda x: self.get_exact_shape(x, c_shp))
                    select.to_file(out_file)
                else:
                    pass
        self.out_gsm_files = out_gsm_files
        
    def rasterize_gsm_vectors(self):
        for vec_file in self.out_gsm_files:
            if os.path.exists(vec_file):
                gsmD = gpd.read_file(vec_file)
                self.rasterize_gsm(gsmD, vec_file.replace(".shp", ".tif"))
    
    def rasterize_gsm(self, gsm, out_tif, cell_width = 0.01):        
        bounds = gsm.total_bounds
        width_cells  = int(round((bounds[2] - bounds[0]) / cell_width))
        height_cells = int(round(((bounds[3] - bounds[1]) / cell_width)))
        cAffine = affine.Affine(cell_width, 0, bounds[0], 0, cell_width*-1, bounds[3])
        cMeta = {'count':1, 'crs': gsm.crs, 'dtype':'uint8', 'affine':cAffine, 'driver':'GTiff',
                 'transform':cAffine, 
                 'height':height_cells, 'width':width_cells}
        shapes = ((row.geometry,1) for idx, row in gsm.iterrows())
        with rasterio.open(out_tif, 'w', **cMeta) as out:
                burned = features.rasterize(shapes=shapes, fill=0., 
                                            out_shape=(cMeta['height'], cMeta['width']), 
                                            transform=out.transform)
                burned = burned.astype(cMeta['dtype'])
                out.write_band(1, burned)

        
    
    
    
