import os, sys, time, subprocess, argparse, logging

import geopandas as gpd
import pandas as pd

from shapely.geometry import box, LineString, Point

def extract_power_plants(power_plants_database, country_bounds):
    pp_d = pd.read_csv(power_plants_database)
    pp_geom = [Point(x) for x in zip(pp_d['longitude'], pp_d['latitude'])]
    pp_d = gpd.GeoDataFrame(pp_d, geometry=pp_geom, crs={'init':'epsg:4326'})
    
    selected_pp = pp_d[pp_d.intersects(country_bounds.unary_union)]
    return(selected_pp)
    
    
def extract_transmission_lines(grid_line_file, out_bounds):
    lines = gpd.read_file(grid_line_file)
    sel_lines = lines[lines.intersects(out_bounds.unary_union)]
    return(sel_lines)