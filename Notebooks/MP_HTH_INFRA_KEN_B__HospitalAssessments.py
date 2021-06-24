import os, sys, importlib, multiprocessing
import rasterio, affine, gdal

import networkx as nx
import geopandas as gpd
import pandas as pd
import numpy as np
import skimage.graph as graph

from shapely.geometry import Point, shape, box
from shapely.wkt import loads
from shapely.ops import cascaded_union
from rasterio import features

import GOSTnets as gn
import GOSTnets.load_osm as losm

sys.path.append("/home/wb411133/Code/GOST")

import GOSTRocks.rasterMisc as rMisc
import GOSTRocks.misc as misc
import GOSTRocks.Urban.UrbanRaster as urban

from GOSTRocks.misc import tPrint

sys.path.append("../")
import infrasap.market_access as ma

def generate_feature_vectors(network_r_path, inH, threshold, featIdx, verbose=True):
    ''' Generate individual market sheds for each feature in the input dataset
    
    INPUTS
        network_r [string] - path to raster from which to grab index for calculations in MCP
        mcp [skimage.graph.MCP_Geometric] - input graph
        inH [geopandas data frame] - geopandas data frame from which to calculate features
        threshold [list of int] - travel treshold from which to calculate vectors in units of graph
        featIdx [string] - column name in inH to append to output marketshed dataset
        
    RETURNS
        [geopandas dataframe]
    '''
    n = inH.shape[0]
    feat_count = 0
    complete_shapes = []
    
    network_r = rasterio.open(network_r_path)
    traversal_time = network_r.read()[0,:,:]
    mcp = graph.MCP_Geometric(traversal_time)
    thread_id = multiprocessing.current_process().name
    
    for idx, row in inH.iterrows():
        feat_count = feat_count + 1
        if verbose:
            tPrint(f"{thread_id}: {feat_count} of {n}")
        cur_idx = network_r.index(row['geometry'].x, row['geometry'].y)
        if cur_idx[0] > 0 and cur_idx[1] > 0 and cur_idx[0] < network_r.shape[0] and cur_idx[1] < network_r.shape[1]:
            costs, traceback = mcp.find_costs([cur_idx])
            for thresh in threshold:
                within_time = ((costs < thresh) * 1).astype('int16')
                all_shapes = []
                for cShape, value in features.shapes(within_time, transform = network_r.transform):
                    if value == 1.0:
                        all_shapes.append([shape(cShape)])
                complete_shape = cascaded_union([x[0] for x in all_shapes])
                complete_shapes.append([complete_shape, thresh, row[featIdx]])
    final = gpd.GeoDataFrame(complete_shapes, columns=["geometry", "threshold", "IDX"], crs=network_r.crs)
    return(final)
    
def calculate_vectors(network_r_path, inH, out_file):
    thresholds = [1800, 3600, 7200, 14400]
    #polygons = generate_feature_vectors_noRasterio(mcp, transform, inCRS, inH, [30,60,300,600], "MASTER_FCO")
    if not os.path.exists(out_file):
        polygons = generate_feature_vectors(network_r_path, inH, thresholds, "MASTER_FCO")
        polygons.to_file(out_file)
    return(True)
        


if __name__ == '__main__':
    p = multiprocessing.Pool(10)
    step_count = 100
    out_folder = "/home/wb411133/data/Country/KEN"
    travel_folder = os.path.join(out_folder, 'TRAVEL_TIMES')
    v_folder = os.path.join(out_folder, "vector_traveltimes")
    if not os.path.exists(v_folder):
        os.makedirs(v_folder)
    network_map = os.path.join(travel_folder, "road_network.tif")
    all_hospitals = "/home/public/Data/COUNTRY/KEN/HD_INF/merged_hospitals.shp"
    
    network_r = rasterio.open(network_map)
    traversal_time = network_r.read()[0,:,:]
    mcp = graph.MCP_Geometric(traversal_time)
    
    inH = gpd.read_file(all_hospitals)
    bad_facilities = ["VCT Centre","Blood Centre","Facility Type","Radiology Unit","Hospice"]
    inH = inH.loc[~inH['type'].isin(bad_facilities)]
    inH['IDX'] = inH['geometry'].apply(lambda x: network_r.index(x.x, x.y))
    
    inH['x'] = inH['IDX'].apply(lambda x: x[0])
    inH['y'] = inH['IDX'].apply(lambda x: x[1])

    inH = inH[~inH['IDX'].duplicated()]
    
    start_idx = 0
    all_inputs = []
    for end_idx in range(step_count, inH.shape[0], step_count):
        out_file = os.path.join(v_folder, f"sample_dt_60_120__{start_idx}_{end_idx}.shp")
        if not os.path.exists(out_file):
            #all_inputs.append([mcp, network_r.transform, network_r.crs, inH.iloc[start_idx:end_idx,].copy(), out_file])
            all_inputs.append([network_map, inH.iloc[start_idx:end_idx,].copy(), out_file])
        start_idx = end_idx
        
    p.starmap(calculate_vectors, all_inputs)
    