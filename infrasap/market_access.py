import sys, os, importlib
import rasterio

import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox
import GOSTnets as gn
import skimage.graph as graph

from rasterio.mask import mask
from rasterio import features
from shapely.geometry import box, Point
from scipy.ndimage import generic_filter
from pandana.loaders import osm

speed_dict = {
   'residential': 20,  # kmph
   'primary': 40,
   'primary_link':35,
   'motorway':50,
   'motorway_link': 45,
   'trunk': 40,
   'trunk_link':35,
   'secondary': 30,
   'secondary_link':25,
   'tertiary':30,
   'tertiary_link': 25,
   'unclassified':20,
   'living_street':10,
   'service':10
}

def get_speed(x, s_dict):
    try: 
        
        speed = s_dict[x]
    except:
        if type(x) == list:
            try:
                speed = s_dict[x[0]]
            except:
                speed = 5
        else:
            speed=5
    return(speed)

def get_nodes(b, tags):
    nodes = osm.node_query(b[1], b[0], b[3], b[2], tags=tags)
    nodes_geom = [Point(x) for x in zip(nodes['lon'], nodes['lat'])]
    nodes_df = gpd.GeoDataFrame(nodes[['amenity','lat','lon']], geometry=nodes_geom, crs={'init':'epgs:4326'})
    return(nodes_df)
    
def get_roads(b):
    sel_graph = ox.graph_from_bbox(b[3], b[1], b[2], b[0], retain_all=True)
    sel_roads = gn.edge_gdf_from_graph(sel_graph)
    sel_roads['speed'] = sel_roads['highway'].apply(lambda x: get_speed(x, speed_dict))
    return(sel_roads)

def generate_network_raster(inH, sel_roads, min_speed=5, speed_col='speed', resolution=100):   
    ''' Create raster with network travel times from a road network
    
    INPUTS
        inH [rasterio object] - template raster used to define raster shape, resolution, crs, etc.
        sel_roads [geopandas dataframe] - road network to burn into raster
        [optional] min_speed [int] - minimum travel speed for areas without roads
        [optional] speed_col [string] - column in sel_roads that defines the speed in KM/h
        [optional] resolution [int] - resolution of the raster in metres
        
    RETURNS
        [numpy array]
    '''
    # create a copy of inH with value set to slowest walking speed
    distance_data = np.zeros(inH.shape)
    # burn the speeds into the distance_data using the road network 
    sel_roads = sel_roads.sort_values([speed_col])
    shapes = ((row['geometry'], row[speed_col]) for idx, row in sel_roads.iterrows())
    speed_image = features.rasterize(shapes, out_shape=inH.shape, transform=inH.transform, fill=min_speed)
    # convert to a version that claculates the seconds to cross each cell
    traversal_time = resolution / (speed_image * 1000 / (60 * 60)) # km/h --> m/s * resolution of image in metres
    return(traversal_time)

def calculate_travel_time(inH, traversal_time, destinations, out_raster = ''):
    ''' Calculate travel time raster
    
    INPUTS
        inH [rasterio object] - template raster used to identify locations of destinations
        traversal_time [numpy array] - describies per pixel seconds to cross
        destinations [geopandas df] - destinations for nearest calculations
        
    LINKS
        https://scikit-image.org/docs/0.7.0/api/skimage.graph.mcp.html#skimage.graph.mcp.MCP.find_costs
    '''
    # create skimage graph
    mcp = graph.MCP_Geometric(traversal_time)
    cities = list(set([inH.index(x.x, x.y) for x in destinations['geometry']]))
    cities = [x for x in cities if ((x[0] > 0) and (x[1] > 0) and 
                (x[0] <= inH.shape[0]) and (x[1] <= inH.shape[1]))]
    
    costs, traceback = mcp.find_costs(cities)        
    if not out_raster == '':
        meta = inH.meta.copy()
        meta.update(dtype=costs.dtype)
        with rasterio.open(out_raster, 'w', **meta) as out:
            out.write_band(1, costs)
            
    return(costs)
    
    
def get_all_amenities(bounds):
    amenities = ['toilets', 'washroom', 'restroom']
    toilets_tags = '"amenity"~"{}"'.format('|'.join(amenities))
    toilets = get_nodes(inH.bounds, toilets_tags)

    amenities = ['water_points', 'drinking_water', 'pumps', 'water_pumps', 'well']
    water_tags = '"amenity"~"{}"'.format('|'.join(amenities))
    water_points = get_nodes(inH.bounds, water_tags)
        
    amenities = ['supermarket', 'convenience', 'general', 'department_stores', 'wholesale', 'grocery', 'general']
    shp_tags = '"shop"~"{}"'.format('|'.join(amenities))
    shops = get_nodes(inH.bounds, shp_tags)
    