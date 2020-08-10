import os, sys, importlib, subprocess, copy, multiprocessing
import rasterio, geohash

import geopandas as gpd
import pandas as pd
import numpy as np

from shapely.geometry import Point
from shapely.wkt import loads
from rasterio import features
from collections import Counter
from multiprocessing import Pool

try:
    from . import vulnerability_mapping as vulmap
    from . import misc
    from . import rasterMisc as rMisc
    from . import UrbanRaster as urban
except:
    import vulnerability_mapping as vulmap 
    import misc
    import rasterMisc as rMisc
    import UrbanRaster as urban
'''
import vulnerability_mapping as vulmap
import misc as misc
import osmMisc as osm
import UrbanRaster as urban
'''
# https://mrc-ide.github.io/global-lmic-reports/parameters.html
# https://www.cdc.gov/mmwr/volumes/69/wr/mm6912e2.htm
vul_def = {'0-5'  :0.001,
           '6-10' :0.001,
           '11-15':0.001,
           '16-20':0.002,
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

hnp_categories = {
    'R10': {
        'Name':'WorldPop',
        'raster_file':'WP_2020_1km.tif',
        'vars':['SUM'],
        'description':'WorldPop total population in 2020'
    },
    'P1': {
        'Name': 'Urban_Population',
        'raster_file':'WP_2020_1km_urban_pop.tif',
        'vars':['SUM'],
        'description':'WorldPop total urban population in 2020'
    },
    'P2': {
        'Name':'Demographics',
        'raster_file':'WP2020_vulnerability_map.tif',
        'vars':['SUM'],
        'description':'Total potential hospitalization load based on WorldPop demographics and published CoVID hospitalization rates'
    },
    'LC': {
        'Name':'Landcover',
        'raster_file':'LC.tif',
        'vars':['C'],
        'unqVals':[11,14,20,30,40,50,60,70,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230],
        'description':'Landcover dataset from Globcover'
    }
} 

def create_fishnet(extents_file, out_folder, prefix, verbose=True):
    ''' Create a 1 km fishnet inside each feature in the input extents_file
    
    INPUT
        extents_file [string] - path to urban extents
        out_folder [string path] - where output shapefiles should be written
        prefix [string] - will be appended to each fidhnet shapefile
    '''
    urban_extents = gpd.read_file(extents_file)
    #sel_cities = urban_extents.sort_values(['Pop'], ascending=False).iloc[0:5]
    sel_cities = urban_extents.sort_values(['Pop'], ascending=False)
    try:
        sel_cities = misc.project_UTM(sel_cities)
    except:
        sel_cities = sel_cities.to_crs({"init":"epsg:3857"})

    for idx, row in sel_cities.iterrows():
        out_fishnet = os.path.join(out_folder, "%s_%s.shp" % (prefix, row['ID']))
        if not os.path.exists(out_fishnet):
            b = row['geometry'].bounds
            crs_num = sel_cities.crs['init'].split(":")[-1]
            crs_num = int(crs_num)
            misc.createFishnet(out_fishnet, b[0], b[2], b[1], b[3], 1000, 1000, crsNum=crs_num)
            fishnet = gpd.read_file(out_fishnet)
            fishnet = fishnet[fishnet.intersects(row['geometry'])]
            fishnet = fishnet.to_crs({'init':'epsg:4326'})
            fishnet['geohash'] = fishnet['geometry'].apply(lambda x: geohash.encode(x.centroid.y, x.centroid.x))
            fishnet.to_file(out_fishnet)
            if verbose:
                misc.tPrint("%s: %s" % (prefix, row['ID']))

def create_urban_data(iso3, country_folder, country_bounds, inR, calc_urban=True, calc_hd_urban=True, verbose=False,
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

    country_pop = os.path.join(country_folder, "WP_2020_1km.tif")
    urban_pop = os.path.join(country_folder, "WP_2020_1km_urban_pop.tif")
    hd_urban_pop = os.path.join(country_folder, "WP_2020_1km_hd_urban_pop.tif")
    urb_bounds = os.path.join(country_folder, "urban_areas.shp")
    hd_urb_bounds = os.path.join(country_folder, "urban_areas_hd.shp")
    urban_fishnets = os.path.join(country_folder, "urban_fishnets")
    hd_urban_fishnets = os.path.join(country_folder, "hd_urban_fishnets")

    # Clip pop raster
    if not os.path.exists(country_pop):
        rMisc.clipRaster(inR, country_bounds, country_pop)

    if not os.path.exists(urban_pop) and calc_urban:
        urbanR = urban.urbanGriddedPop(country_pop)
        urban_vec = urbanR.calculateUrban(densVal = urb_dens, 
                                          totalPopThresh = urb_pop, 
                                          smooth = True,
                                          queen = False,
                                          raster_pop = urban_pop)
        urban_vec.to_file(urb_bounds)

    if not os.path.exists(hd_urban_pop) and calc_hd_urban:
        urbanR = urban.urbanGriddedPop(country_pop)
        urban_vec = urbanR.calculateUrban(densVal = hd_urb_dens, 
                                          totalPopThresh = hd_urb_pop, 
                                          smooth = True,
                                          queen = True,
                                          raster_pop = urban_pop)
        urban_vec.to_file(hd_urb_bounds)
    #Generate Fishnets
    if not os.path.exists(urban_fishnets):
        os.makedirs(urban_fishnets)

    if not os.path.exists(hd_urban_fishnets):
        os.makedirs(hd_urban_fishnets)

    if calc_urban:
        create_fishnet(urb_bounds, urban_fishnets, "URBAN")
    if calc_hd_urban:
        create_fishnet(hd_urb_bounds, hd_urban_fishnets, "HD_URBAN")
           
def calculate_vulnerability(iso3, country_folder, country_bounds, pop_folder, pop_files):
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
    
    # clip out the temporary vulnerability metrics
    out_vulnerability = os.path.join(country_folder, "WP2020_vulnerability_map.tif")
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

def extract_osm(sel_country, out_folder,
                global_osm_pbf = '/home/public/Data/GLOBAL/OSM/GLOBAL/planet-latest.osm.pbf'):
    ''' Extract a national osm pbf from the global file
        INPUTS
            sel_country [shapely object] - extent to clip
            out_folder [string folder path] - output folder for creating the osm.osm.pbf
    '''
     
    osm_extractor = osm.osmExtraction(osmosisCmd = "/home/wb411133/Code/Osmosis/bin/osmosis",
                                        tempFile = '/home/wb411133/temp/osmosis.sh')
    out_pbf = os.path.join(out_folder, "osm.osm.pbf")
    out_roads = os.path.join(out_folder, "OSM_roads.shp")
    if not os.path.exists(out_pbf):
        cmds = osm_extractor.extractBoundingBox(global_osm_pbf, sel_country, out_pbf, execute=False)
        subprocess.check_call(cmds.split(" "))
    if not os.path.exists(out_roads):
        roads = osm.convertOSMPBF_DataFrame(out_pbf, "lines")
        roads_sel = roads.loc[:,["geometry","highway","Length","osm_id"]]
        roads_sel = roads_sel.loc[[not x is None for x in roads['highway']]]
        roads_sel.to_file(out_roads)
        
def run_zonal(admin_shapes, rasters, out_suffix='', iso3=''):
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
            out_suffix [string] - text to append to output zonal file
    '''
    for shp in admin_shapes:
        inD = gpd.read_file(shp)
        out_zonal = shp.replace(".shp", "_zonal%s.csv" % out_suffix)
        misc.tPrint(f"Processed: {iso3} {os.path.basename(shp)}")
        write_out = False
        if not os.path.exists(out_zonal):
            for var_name, definition in rasters.items():
                if os.path.exists(definition['raster_file']):
                    write_out = True
                    if definition['vars'][0] == 'C':
                        uVals = definition['unqVals']
                        res = rMisc.zonalStats(inD, definition['raster_file'], rastType='C', unqVals=uVals, reProj=True)                    
                        res = pd.DataFrame(res, columns=['LC_%s' % x for x in uVals])
                        for column in res.columns:
                            inD[column] = res[column]
                    else:
                        # Zonal stats
                        res = rMisc.zonalStats(inD, definition['raster_file'], minVal=0, reProj=True)
                        res = pd.DataFrame(res, columns=['SUM','MIN','MAX','MEAN'])
                        res.columns = [f"{var_name}_{x}" for x in res.columns]
                        for var in definition['vars']:
                            inD[f"{var_name}_{var}"] = res[f"{var_name}_{var}"]                
            if write_out:
                inD.drop(['geometry'], axis=1, inplace=True)
                pd.DataFrame(inD).to_csv(out_zonal)

def check_zonal(country_folder, remove_bad = False):
    ''' Check the results of the zonal statistics
    '''
    
    stat_files = []
    for root, dirs, files in os.walk(country_folder):
        if not "FINAL" in root:
            for f in files:
                if f[-4:] == ".csv":
                    stat_files.append(os.path.join(root, f))
        
    res = {}
    for stat_file in stat_files:
        if "DHS" in stat_file:
            good_col = "age_final_0_4_househ_SUM"
            bad_col = "R10_SUM"
        else:
            good_col = "R10_SUM"
            bad_col = "age_final_0_4_househ_SUM"
        xx = pd.read_csv(stat_file)
        cols = list(xx.columns)
        val = 0
        if good_col in cols:
            val += 1
        if bad_col in cols:
            val += 10
        res[stat_file] = val
    
    if remove_bad:
        for f, score in res.items():
            if score in [0, 10]:
                os.remove(f)
            
    return(res)
            
def combine_dhs_pop(popRaster, dhs_raster, out_file, factor=100):
    ''' 
    INPUT
        popRaster [rasterio]
        dhs_raster [rasterio]
        out_file [string]
        [optional] factor [int] - number to divide dhs_raster (converts percentage to fraction)
    '''
    inP = popRaster.read()
    dhs = dhs_raster.read()
    if factor != 1:
        dhs = dhs / factor
    dhs_pop = inP * (dhs)
    
    with rasterio.open(out_file, 'w', **popRaster.meta) as outR:
        outR.write(dhs_pop)
            
def summarize_DHS(template, dhs_files, country_folder, iso3):
    ''' combine DHS data with WorldPop population and run zonal stats
    
    INPUT
        template [string] - template raster upon which to base rasterization of DHS
        dhs_files [dictionary] - defines DHS files to process {filename:geopandas}
        country_folder [string to path] - folder to create output
    '''
    # Process DHS data
    inP = rasterio.open(template)
    # get a list of unique columns in the DHS data
    total_columns = 0
    try:
        del(all_columns)
    except:
        pass
    # get a list of all unique columns
    for key, inD in dhs_files.items():
        cur_columns = list(inD.columns.values)
        try:
            all_columns = all_columns + cur_columns
        except:
            all_columns = cur_columns

    col_count = Counter(all_columns)
    unq_columns = [key for key, value in col_count.items() if value == 1] 
    
    dhs_rasters = {}
    for key, inD in dhs_files.items():
        sel_dhs = inD.loc[inD['ISO3'] == iso3]
        if sel_dhs.shape[0] > 0:
            for field in inD.columns:
                if field in unq_columns:
                    out_file = os.path.join(country_folder, f'{key}_{field}.tif')
                    out_file_pop = os.path.join(country_folder, f'{key}_{field}_pop.tif')
                    try:
                        # rasterize the desired field in the inputDHS data                
                        if not os.path.exists(out_file) and not os.path.exists(out_file_pop):
                            rMisc.rasterizeDataFrame(inD, out_file, idField=field, templateRaster = template)

                        #Multiply the rasterized data frame by the population layer
                        if not os.path.exists(out_file_pop):
                            combine_dhs_pop(inP, rasterio.open(out_file), out_file_pop, factor=100)
                        if os.path.exists(out_file):
                            os.remove(out_file)
                        misc.tPrint(f'{iso3}_{key}: {field}')
                        dhs_rasters[f'{key}_{field}'] = {
                            'raster_file': f'{key}_{field}_pop.tif',
                            'vars': ['SUM', 'MEAN'],
                            'description': f'{key}_{field}'
                        }
                    except:
                        misc.tPrint(f"Error processing {key} - {field}")
    return(dhs_rasters)
 
def extract_data(inG, inG1, inG2, inL, inR):
    country_folder = os.path.join(output_folder, iso3)
    adm0_file = os.path.join(country_folder, "adm0.shp")
    adm1_file = os.path.join(country_folder, "adm1.shp")
    adm2_file = os.path.join(country_folder, "adm2.shp")
    lc_file = os.path.join(country_folder, "LC.tif")
           
    if not os.path.exists(country_folder):
        os.makedirs(country_folder)
    country_bounds = inG.loc[inG['ISO3'] == iso3].to_crs({'init':'epsg:4326'})
    if not os.path.exists(adm0_file):
        country_bounds.to_file(adm0_file)
    if not os.path.exists(adm1_file):
        try:
            country_adm1   = inG1.loc[inG1['ISO3'] == iso3].to_crs({'init':'epsg:4326'})
            country_adm1.to_file(adm1_file)
        except:
            misc.tPrint("%s Could not extract ADMIN 1" % iso3)
    if not os.path.exists(adm2_file):
        try:
            country_adm2   = inG2.loc[inG2['ISO3'] == iso3].to_crs({'init':'epsg:4326'})
            country_adm2.to_file(adm2_file)
        except:
            misc.tPrint("%s Could not extract ADMIN 2" % iso3)
    if not os.path.exists(lc_file):
        rMisc.clipRaster(inL, gpd.read_file(adm0_file), lc_file)
        
    calculate_vulnerability(iso3, country_folder, country_bounds, pop_folder, pop_files)
    misc.tPrint("***%s Calculated Vulnerability" % iso3)
    try:
        create_urban_data(iso3, country_folder, country_bounds, inR, calc_urban=False)
        misc.tPrint("***%s Calculated Urban Extents" % iso3)                           
    except:
        misc.tPrint("%s errored on HD clusters" % iso3)
        try:
            create_urban_data(iso3, country_folder, country_bounds, inR, calc_urban=True, calc_hd_urban=False)
        except:
            misc.tPrint("%s errored on all clusters" % iso3)        
    #extract_osm(country_bounds, country_folder)
    #misc.tPrint("***Extracted OSM")
    
 
def run_all(iso3, output_folder, dhs_files):
    country_folder = os.path.join(output_folder, iso3)
    # extract national bounds
    misc.tPrint("Processing %s" % iso3)        
    #summarize DHS
    country_pop = os.path.join(country_folder, "WP_2020_1km.tif")
    dhs_rasters = summarize_DHS(country_pop, dhs_files, country_folder, iso3)
    
    #Run zonal stats
    cur_rasters = copy.deepcopy(hnp_categories)
    for key, values in cur_rasters.items():
        values['raster_file'] = os.path.join(country_folder, values['raster_file'])
        cur_rasters[key] = values
    
    cur_dhs = copy.deepcopy(dhs_rasters)
    for key, values in dhs_rasters.items():
        values['raster_file'] = os.path.join(country_folder, values['raster_file'])
        cur_dhs[key] = values        
    
    all_shps = []
    for root, dirs, files, in os.walk(country_folder):
        for f in files:
            if f[-4:] == ".shp" and not "zonal" in f:
                all_shps.append(os.path.join(root, f))
                
    run_zonal(all_shps, cur_rasters, out_suffix="_BASE", iso3 = iso3)
    misc.tPrint("***%s Calculated Base Zonal" % iso3)
    run_zonal(all_shps, cur_dhs, out_suffix="_DHS", iso3 = iso3)
    misc.tPrint("***%s Calculated DHS Zonal" % iso3)
 
def main():
    # define the input datasets
    global_bounds = "/home/public/Data/GLOBAL/ADMIN/Admin0_Polys.shp"
    global_adm1 = "/home/public/Data/GLOBAL/ADMIN/Admin1_Polys.shp"
    global_adm2 = "/home/public/Data/GLOBAL/ADMIN/Admin2_Polys.shp"
    pop_folder = "/home/public/Data/GLOBAL/Population/WorldPop_PPP_2020/GLOBAL_1km_Demographics"
    output_folder = "/home/wb411133/data/Projects/CoVID"
    population_raster = "/home/public/Data/GLOBAL/Population/WorldPop_PPP_2020/ppp_2020_1km_Aggregated.tif" 
    lcRaster = "/home/public/Data/GLOBAL/LANDCOVER/GLOBCOVER/2015/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif"
    dhs_folder = '/home/public/Data/PROJECTS/CoVID/DHS'
    dhs_files = {}
    for root, dirs, files in os.walk(dhs_folder):
        for f in files:
            if f[-4:] == ".shp":
                dhs_files[f.replace(".shp", "")] = gpd.read_file(os.path.join(os.path.join(root, f)))
    
    # Read in the global datasets
    pop_files = os.listdir(pop_folder)
    inG  = gpd.read_file(global_bounds)
    inG1 = gpd.read_file(global_adm1)
    inG2 = gpd.read_file(global_adm2)
    inR = rasterio.open(population_raster)
    inL = rasterio.open(lcRaster)

    countries = ['ARG', 'PAK', 'ZAF', 'COL', 'ZWE', 'MNG', 'SLE', 'CPV', 'KEN', 'GHA', 'AFG', 'YEM', 'ECU', 'PRY', 'MRT', 'MDV', 'KGZ', 'HTI', 'DJI', 'KHM', 'TJK', 'GMB', 'LKA', 'SEN', 'STP', 'SLV', 'VEN', 'MLI', 'RWA', 'BOL', 'TZA', 'MAR', 'IND', 'IDN', 'SDN', 'AGO', 'BEN', 'BWA', 'BFA', 'BDI', 'CMR', 'CAF', 'TCD', 'COM', 'COG', 'CIV', 'COD', 'SSD', 'ERI', 'ETH', 'GAB', 'GNB', 'GIN', 'LSO', 'LBR', 'MDG', 'MWI', 'MUS', 'MOZ', 'NAM', 'NER', 'NGA', 'SYC', 'SOM', 'SWZ', 'TGO', 'UGA', 'ZMB', 'LCA', 'PHL', 'GTM', 'BGD', 'BRA', 'MEX', 'EGY', 'UKR', 'PER', 'LAO', 'PSE', 'NPL', 'PNG', 'DZA', 'BLR', 'BTN', 'BIH', 'NIC', 'FJI', 'GEO', 'HND', 'JOR', 'MHL', 'MDA', 'MMR', 'MKD', 'PAN', 'WSM', 'SLB', 'TUN', 'TUR', 'URY', 'UZB', 'ALB', 'HRV', 'IRN', 'SRB', 'TTO', 'ATG', 'CHN', 'IRQ']
    countries = ['EGY','PHL','BGD']
    countries = set(countries)
    nCountries = len(countries)
    idx = 0
    all_commands = []
    for iso3 in countries:
        all_commands.append([iso3, output_folder, dhs_files])
        
    with Pool(round(multiprocessing.cpu_count() * 0.8)) as p:
        res = p.starmap(run_all, all_commands)

if __name__ == "__main__":
    main()
    
