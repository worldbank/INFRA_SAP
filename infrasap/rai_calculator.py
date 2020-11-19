import os, sys, time, subprocess, argparse, logging
import rasterio.mask

import geopandas as gpd
import osmnx as ox
import pandas as pd

import GOSTnets as gn
import GOSTnets.load_osm as losm

from . import osm_extractor as osm
from . import rasterMisc as rMisc

from shapely.geometry import box

def extract_rai_network(focal_osm, epsg=3857, rai_buffer=2000):
    ''' extract road network from OSM PBF
    
    INPUTS
    focal_osm [string] - path to osm pbf to be used for calculations
    epsg [int] - epsg code for metres-based measurement codes
    rai_buffer [number]
    
    RETURNS
    [geopandas data frame]
    '''
    G_loader = losm.OSM_to_network(focal_osm)
    G_loader.generateRoadsGDF()
    roadsGPD = G_loader.roadsGPD
    roadsGPD.set_geometry("Wkt", inplace=True)
    roadsGPD = roadsGPD.to_crs({'init':'epsg:%s' % epsg})
    roadsGPD_rai = roadsGPD.copy()
    roadsGPD_rai['Wkt'] = roadsGPD_rai['Wkt'].apply(lambda x: x.buffer(rai_buffer))
    
    roadsGPD_rai['OSMLR'] = roadsGPD_rai['infra_type'].map(osm.OSMLR_Classes)
    def get_num(x):
        try:
            return(int(x))
        except:
            return(5)
    roadsGPD_rai['OSMLR_num'] = roadsGPD_rai['OSMLR'].apply(lambda x: get_num(str(x)[-1]))
    
    return(roadsGPD_rai)
    
def calculate_rai(out_bounds, out_bounds_id, wp_data, roadsGPD_rai, rai_folder):
    ''' Calculate RAI
    
    INPUT
    out_bounds [geopandas dataframe] - admin boundaries within which to summarize RAI
    out_bounds_id [string] - column name in out_bounds containing INDETIFIER
    wp_data [rasterio raster] - gridded population dataset from which to calculate RAI
    roadsGPD_rai [geopandas data frame] - generated from function .extract_rai_network()
    
    RETURNS
    [pandas dataframe]
    '''
    actual_pop = rMisc.zonalStats(out_bounds, wp_data, minVal=0)
    actual_pop = pd.DataFrame(actual_pop, columns=['SUM','MIN','MAX','STDEV'])
    # Generate OSMLR level road maps
    all_res = {}
    for road_count in [1,2,3,4]:
        rai_pop = os.path.join(rai_folder, "WP_RAI_%s.tif" % road_count)
        # Generate pop map within RAI pop
        if not os.path.exists(rai_pop):
            cur_roads = roadsGPD_rai.loc[roadsGPD_rai['OSMLR_num'] <= road_count]
            roads_mask = rasterio.mask.mask(wp_data, [cur_roads.unary_union], invert=False)
            with rasterio.open(rai_pop, 'w', **wp_data.meta) as out:
                out.write(roads_mask[0])
        rai_pop = rMisc.zonalStats(out_bounds, rai_pop, minVal=0)

        rai_pop = pd.DataFrame(rai_pop, columns=['SUM','MIN','MAX','STDEV'])
        rai_pop['ID'] = out_bounds[out_bounds_id]
        rai_pop = rai_pop[['SUM','ID']]
        rai_pop.columns = ['RAI_POP_%s' % road_count, 'ID']
        all_res[road_count] = rai_pop    
        
    final = all_res[1]
    final['RAI_POP_2'] = all_res[2]['RAI_POP_2']
    final['RAI_POP_3'] = all_res[3]['RAI_POP_3']
    final['RAI_POP_4'] = all_res[4]['RAI_POP_4']
    final['POP'] = actual_pop['SUM']
    return(final)
