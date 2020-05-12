import os, sys, importlib, subprocess, logging, json
import rasterio

import geopandas as gpd
import pandas as pd
import numpy as np

from shapely.geometry import Point
from shapely.wkt import loads
from rasterio import features

import vulnerability_mapping as vulmap

import GOSTRocks.rasterMisc as rMisc
import GOSTRocks.misc as misc
import GOSTRocks.osmMisc as osm
import GOSTRocks.Urban.UrbanRaster as urban

# https://mrc-ide.github.io/global-lmic-reports/parameters.html
# https://www.cdc.gov/mmwr/volumes/69/wr/mm6912e2.htm
vul_def = {'0-5'  :0.001,
           '6-10' :0.001,
           '11-15':0.001,
           '16-20':0.002,0
           '21-25':0.005,
           '26-30':0.010,
           '31-35':0.016,
           '36-40':0.023,
           '41-45':0.029,
           '46-50':0.039,
           '51-55':0.058,
           '56-60':0.072,
           '61-65':0.102,
           '66-70':0.117,
           '71-75':0.146,
           '76-80':0.177,
           '81-100':0.180}

def set_up_logging(fileName):
    logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()

    fileHandler = logging.FileHandler("{0}.log".format(fileName))
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
           
class covid_extractor(object):
    ''' Extract and process base information for assessing CoVID vulnerabilities
        1. Urban areas
        2. Population density
        3. Vulnerability mapping
        4. Fraym data
    '''
    def __init__(self, iso3, country_bounds, inG1, inG2, data_folder):
        self.iso3 = iso3
        self.out_folder = data_folder
        self.country_bounds = country_bounds
        
        #Define data layers
        self.adm0_file = os.path.join(data_folder, "adm0.shp")
        if not os.path.exists(self.adm0_file):
            self.country_bounds.to_file(self.adm0_file)
        
        self.adm1_file = os.path.join(data_folder, "adm1.shp")
        self.adm2_file = os.path.join(data_folder, "adm2.shp")
        try:
            self.country_adm1 = inG1.loc[inG1['ISO3'] == iso3].to_crs({'init':'epsg:4326'})
            if not os.path.exists(self.adm1_file):
                self.country_adm1.to_file(self.adm1_file)
        except:
            logging.warning("Could not create admin 1 for %s" % self.iso3)
        try:    
            self.country_adm2 = inG2.loc[inG2['ISO3'] == iso3].to_crs({'init':'epsg:4326'})
            if not os.path.exists(self.adm2_file):
                self.country_adm2.to_file(self.adm2_file)
        except:
            logging.warning("Could not create admin 2 for %s" % self.iso3)
        
        self.country_pop = os.path.join(data_folder, "WP_2020_1km.tif")
        self.urban_pop = os.path.join(data_folder, "WP_2020_1km_urban_pop.tif")
        self.urb_bounds = os.path.join(data_folder, "urban_areas.shp")
        self.hd_urb_bounds = os.path.join(data_folder, "urban_areas_hd.shp")
        self.urban_fishnets = os.path.join(data_folder, "urban_fishnets")
        self.hd_urban_fishnets = os.path.join(data_folder, "hd_urban_fishnets")
        if not os.path.exists(self.urban_fishnets):
            os.makedirs(self.urban_fishnets)
        if not os.path.exists(self.hd_urban_fishnets):
            os.makedirs(self.hd_urban_fishnets)
        self.out_vulnerability = os.path.join(data_folder, "WP2020_vulnerability_map.tif")
        self.out_pbf = os.path.join(data_folder, "osm.osm.pbf")
        self.out_roads = os.path.join(data_folder, "OSM_roads.shp")  
        self.check_data()
        
    def check_data(self, execute=False):
        ''' Check all the input data to determine what needs to be calculated
        '''
        adm_data = False
        if os.path.exists(self.adm0_file) and os.path.exists(self.adm1_file) and os.path.exists(self.adm2_file):
            adm_data = True
        
        vul_data = False
        if os.path.exists(self.out_vulnerability):
            vul_data = True
        
        urb_extents = False
        if (os.path.exists(self.urb_bounds)) or (os.path.exists(self.hd_urb_bounds)) :
            urb_extents = True
            
        urb_fishnets = False
        if (len(os.listdir(self.urban_fishnets)) > 1) or (len(os.listdir(self.hd_urban_fishnets)) > 1):
            urb_fishnets = True
            
        self.data_status = {
            "adm_data":adm_data,
            "vul_data":vul_data,
            "urb_extents":urb_extents,
            "urb_fishnets":urb_fishnets,
        }
        
        if not adm_data:
            logging.warning(f"{self.iso3} admin data: {adm_data}")
        if not vul_data:
            logging.warning(f"{self.iso3} vulnerability data: {vul_data}")
        if not urb_extents:
            logging.warning(f"{self.iso3} urban extents: {urb_extents}")
        if not urb_fishnets:
            logging.warning(f"{self.iso3} urban fishnets: {urb_fishnets}")
                       
        
    def create_urban_data(self, inR, calc_urban=True, calc_hd_urban=True, verbose=False,
                            urb_dens=300, urb_pop=5000, hd_urb_dens=1500, hd_urb_pop=50000):
        ''' Extract urban areas from gridded population data following the EC density methodology
        
        INPUTS:
            iso3 [string] - iso3 code for country of interest
            country_folder [string] - path to output folder 
            country_bounds [geopandas dataframe] - boundary within which to extract
            inR [rasterio object] - population datasets from which country specific pop layer is extracted
            [optional] calc_urban [boolean] - whether to calculate urban (lower density) extents
            [optional] calc_hd_urban [boolean] - whether to calculate hd urban (higher density) extents
            [optional] urb_dens/hd_urb_dens [int] - population density threshold for calculating urban areas IN THE PER PIXEL UNITS OF inR
            [optional] urb_pop/hd_urb_pop [int] - total population threshold for calculating urban areas
            
        RETURNS:
            NA - the function calculates a number of files in the folder country_folder:
                -- WP_2020_1km.tif: population dataaset clipped from inR
                -- urban_areas.shp: vectorization of urban areas calculation
                -- urban_areas_hd.shp: vectorization of high density urban areas calculation
                -- urban_fishnets/URBAN_XX.shp: a 1 km fishnet is created over the five highest population areas in urban_areas.shp 
                -- hd_urban_fishnets/HD_URBAN_XX.shp: a 1 km fishnet is created over the five highest population areas in hd_urban_areas.shp 
        '''

        country_pop = self.country_pop
        urban_pop = self.urban_pop
        urb_bounds = self.urb_bounds
        hd_urb_bounds = self.hd_urb_bounds
        urban_fishnets = self.urban_fishnets
        hd_urban_fishnets = self.hd_urban_fishnets
        iso3 = self.iso3
        country_folder = self.out_folder 
        country_bounds = self.country_bounds

        # Clip pop raster
        if not os.path.exists(country_pop):
            rMisc.clipRaster(inR, country_bounds, country_pop)

        if not os.path.exists(urb_bounds) and calc_urban:
            urbanR = urban.urbanGriddedPop(country_pop)
            urban_vec = urbanR.calculateUrban(densVal = urb_dens, 
                                              totalPopThresh = urb_pop, 
                                              smooth=True,
                                              queen=False,
                                              raster_pop = urban_pop)
            if urban_vec.shape[0] > 0:
                urban_vec.to_file(urb_bounds)

        if not os.path.exists(hd_urb_bounds) and calc_hd_urban:
            urbanR = urban.urbanGriddedPop(country_pop)
            urban_vec = urbanR.calculateUrban(densVal = hd_urb_dens, 
                                              totalPopThresh = hd_urb_pop, 
                                              smooth=True,
                                              queen=True,
                                                raster_pop = urban_pop)
            if urban_vec.shape[0] > 0:
                urban_vec.to_file(hd_urb_bounds)
        #Generate Fishnets
        if not os.path.exists(urban_fishnets):
            os.makedirs(urban_fishnets)

        if not os.path.exists(hd_urban_fishnets):
            os.makedirs(hd_urban_fishnets)

        def create_fishnet(extents_file, out_folder, prefix):
            urban_extents = gpd.read_file(extents_file)
            sel_cities = urban_extents.sort_values(['Pop'], ascending=False).iloc[0:5]
            try:
                sel_cities = misc.project_UTM(sel_cities)
            except:
                sel_cities = sel_cities.to_crs({"init":"epsg:3857"})

            for idx, row in sel_cities.iterrows():
                out_fishnet = os.path.join(out_folder, "%s_%s.shp" % (prefix, row['ID']))
                b = row['geometry'].bounds
                crs_num = sel_cities.crs['init'].split(":")[-1]
                crs_num = int(crs_num)
                misc.createFishnet(out_fishnet, b[0], b[2], b[1], b[3], 1000, 1000, crsNum=crs_num)
                fishnet = gpd.read_file(out_fishnet)
                fishnet = fishnet[fishnet.intersects(row['geometry'])]
                fishnet = fishnet.to_crs({'init':'epsg:4326'})
                fishnet.to_file(out_fishnet)
                if verbose:
                    misc.tPrint("%s: %s" % (prefix, row['ID']))

        if calc_urban:
            create_fishnet(urb_bounds, urban_fishnets, "URBAN")
        if calc_hd_urban:
            create_fishnet(hd_urb_bounds, hd_urban_fishnets, "HD_URBAN")
           
    def calculate_vulnerability(self, pop_folder, pop_files):
        ''' The hospitalization rates listed in the vul_def dictionary (at the top of this script) are
            combined to create a vulnerability layer
            
            INPUTS
                iso3 [string] - iso3 code for country of interest
                country_folder [string] - path to output folder 
                country_bounds [geopandas dataframe] - boundary within which to extract
                pop_files [list of file paths] - list of global demographic rasters which are clipped out with country_bounds
            RETURNS
                NA - creates a file called WP2020_vulnerability_map.tif in country_folder            
        '''
        iso3 = self.iso3
        country_folder = self.out_folder
        country_bounds = self.country_bounds
        out_vulnerability = self.out_vulnerability
        # clip out the temporary vulnerability metrics    
        if not os.path.exists(out_vulnerability):
            wp_files = []
            if not os.path.exists(country_folder):
                os.makedirs(country_folder)
            for pFile in pop_files:
                curR = rasterio.open((os.path.join(pop_folder, pFile)))
                out_file = os.path.join(country_folder, pFile)
                if not os.path.exists(out_file):
                    rMisc.clipRaster(curR, country_bounds, out_file)
                wp_files.append(out_file)
            #Calculate vulnerability
            wp_file_objects  = [vulmap.wp_demographics(os.path.join(country_folder, x)) for x in wp_files]
            vul = vulmap.wp_vulnerability(wp_file_objects, vul_def)
            vul.calculate_vulnerability()
            vul.combine_results(out_vulnerability, '')
            for f in wp_files:
                os.remove(f)

    def extract_osm(self, global_osm_pbf = '/home/public/Data/GLOBAL/OSM/GLOBAL/planet-latest.osm.pbf'):
        ''' Extract a national osm pbf from the global file
            INPUTS
                sel_country [shapely object] - extent to clip
                out_folder [string folder path] - output folder for creating the osm.osm.pbf
        '''
        sel_country = self.country_bounds.unary_union
        out_folder = self.out_folder
        osm_extractor = osm.osmExtraction(osmosisCmd = "/home/wb411133/Code/Osmosis/bin/osmosis",
                                            tempFile = '/home/wb411133/temp/osmosis.sh')
        if not os.path.exists(self.out_pbf):
            cmds = osm_extractor.extractBoundingBox(global_osm_pbf, sel_country, self.out_pbf, execute=False)
            subprocess.check_call(cmds.split(" "))
        if not os.path.exists(self.out_roads):
            roads = osm.convertOSMPBF_DataFrame(self.out_pbf, "lines")
            roads_sel = roads.loc[:,["geometry","highway","Length","osm_id"]]
            roads_sel = roads_sel.loc[[not x is None for x in roads['highway']]]
            roads_sel.to_file(self.out_roads)
        
    def run_zonal(self, file_defs):
        ''' Calculate zonal results for submitted admin and raster
            
            INPUTS
                admin_shapes [geopandas] - features within which to calculate statistics
                rasters [dictionary] - data dictionary containing the raster and the required information
                    { 'HNP_Var1':{
                            'raster_file': 'path_to_raster',
                            'vars':['SUM','MEAN'],
                            'description':'Lorem Ipsum'
                        }
                    }
        '''
        admin_shapes = [self.adm0_file, self.adm1_file, self.adm2_file]
        if os.path.exists(self.urb_bounds):
            admin_shapes.append(self.urb_bounds)
        if os.path.exists(self.hd_urb_bounds):
            admin_shapes.append(self.hd_urb_bounds)
        #Add fishnet files to list
        for folder in [self.urban_fishnets, self.hd_urban_fishnets]:
            for root, dirs, files in os.walk(folder):
                for f in files:
                    if f[-4:] == ".shp":
                        admin_shapes.append(os.path.join(root, f))     
                
        for shp in admin_shapes:
            if os.path.exists(shp) and not "zonal" in shp:
                inD = gpd.read_file(shp)
                out_zonal = shp.replace(".shp", "_zonal.shp")
                if not os.path.exists(out_zonal):
                    for var_name, definition in file_defs.items():
                        if definition['raster_file'] != '_INPUT_':
                            rFile = os.path.join(self.out_folder, definition['raster_file'])
                            if os.path.exists(rFile):
                                # Zonal stats
                                res = rMisc.zonalStats(inD, rFile, minVal=0, reProj=True)
                                res = pd.DataFrame(res, columns=['SUM','MIN','MAX','MEAN'])
                                res.columns = [f"{var_name}_{x}" for x in res.columns]
                                for var in definition['vars']:
                                    inD[f"{var_name}_{var}"] = res[f"{var_name}_{var}"]
                    inD.to_file(out_zonal)

def main():
    # define the input datasets
    global_bounds = "/home/public/Data/GLOBAL/ADMIN/Admin0_Polys.shp"
    global_adm1 = "/home/public/Data/GLOBAL/ADMIN/Admin1_Polys.shp"
    global_adm2 = "/home/public/Data/GLOBAL/ADMIN/Admin2_Polys.shp"
    pop_folder = "/home/public/Data/GLOBAL/Population/WorldPop_PPP_2020/GLOBAL_1km_Demographics"
    output_folder = "/home/wb411133/data/Projects/CoVID"
    population_raster = "/home/public/Data/GLOBAL/Population/WorldPop_PPP_2020/ppp_2020_1km_Aggregated.tif" 
    # Read in raster definitions
    json_file = 'DataDictionary_v2.json'
    with open(json_file, 'r') as j_file:
        file_defs = json.load(j_file, strict=False)
    # Read in the global datasets
    pop_files = os.listdir(pop_folder)
    inG  = gpd.read_file(global_bounds)
    inG1 = gpd.read_file(global_adm1)
    inG2 = gpd.read_file(global_adm2)
    inR = rasterio.open(population_raster)
    
    #set up logging
    set_up_logging("covid_extractor.log")
    
    # Read in data processing file
    data_def_file = "../Notebooks/HNP_Data_Readiness_0no_1yes_Ben.csv"
    inD = pd.read_csv(data_def_file, index_col=0)
    countries = inD['NOTHING'].iloc[1:].values
    nCountries = len(countries)
    idx = 0
    #countries = ['IDN']
    for iso3 in countries:
        # extract national bounds
        misc.tPrint("Processing %s of %s: %s" % (idx, nCountries, iso3))
        idx = idx + 1
        country_folder = os.path.join(output_folder, iso3)
        if not os.path.exists(country_folder):
            os.makedirs(country_folder)
        country_bounds = inG.loc[inG['ISO3'] == iso3]
        country_bounds = country_bounds.to_crs({'init':'epsg:4326'})
        covid_extract = covid_extractor(iso3, country_bounds, inG1, inG2, country_folder)
        if not covid_extract.data_status['vul_data']:
            covid_extract.calculate_vulnerability(pop_folder, pop_files)
        
        if not covid_extract.data_status['urb_extents']:
            covid_extract.create_urban_data(inR)

        covid_extract.run_zonal(file_defs)
        
if __name__ == "__main__":
    main()
    
