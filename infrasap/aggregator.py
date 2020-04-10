'''
The following module contains a number of functions to aggregate geospatial outputs into tables for InfraSAP analytics.
'''
import geopandas as gpd
import pandas as pd
import rasterio as rio
from rasterio import features
from rasterstats import zonal_stats
import numpy as np

def rasterize_gdf(inD, field, template, outFile=None, nodata=np.nan):
    ''' Convert geopandas geo data frame to raster of equal size/res to template raster
    
    INPUT
    inD [ geopandas data frame ]
    outFile [ string ] - path to save output raster
    field [ string ] - field to rasterize
    template [ string ] - path to template raster
    nodata [ int ] - value for no data
    
    RETURNS
    Raster file. If no outFile is specified, the function returns burned features as numpy array
    '''
    raster_template = rio.open(template)
    # get info from template file
    xRes = raster_template.res[1]
    yRes = raster_template.res[0]
    trans = raster_template.transform
    x_pixels = raster_template.shape[1]
    y_pixels = raster_template.shape[0]
    raster_template.close()
    
    shapes = ((row.geometry,row[field]) for idx, row in inD.iterrows())
    burned = features.rasterize(shapes=shapes, fill=nodata, out_shape=raster_template.shape, transform=trans)
    burned = burned.astype(str(inD[field].dtype))
    
    if outFile:
        with rio.open(
            outFile, 'w', driver = 'GTiff',
            height = y_pixels, width = x_pixels,
            count=1, dtype=str(inD[field].dtype),
            crs=raster_template.crs,
            transform=trans
        ) as new_dataset:
            new_dataset.write_band(1, burned)
    else:
        return burned
    
def pop_weighted_average(target, data_raster, pop_raster, new_field):
    '''Calculate population weighted average from a raster dataset to a target shapefile
    
    INPUT
    target [ GDF ] - Target geo data frame
    data_raster [ path ]
    pop_raster [ path ]
    new_field [ string ] -
    
    RETURNS
    Target GDF with new fields
    '''
    # Load inputs
    pop = rio.open(pop_raster)
    pop_crs = pop.crs.to_string()
    pop_array = pop.read(1, masked=True)
    
    if pop_crs!=target.crs.to_string():
        target = target.to_crs(pop_crs)
    
    # Calculate weights
    zs_sum_pop = pd.DataFrame(zonal_stats(target, pop_array, affine=pop.transform, stats='sum', nodata=pop.nodata)).rename(columns={'sum':'pop_sum'})
    target = target.join(zs_sum_pop)
    pop_sum_array = rasterize_gdf(inD=target, field='pop_sum',template=pop_raster)
    weights = pop_array/pop_sum_array
    
    # Apply weights
    data = rio.open(data_raster).read(1)
    data_weighted = data*weights
    
    zs_sum_data = pd.DataFrame(zonal_stats(target, data_weighted, affine=pop.transform, stats='sum', nodata=pop.nodata)).rename(columns={'sum':new_field})
    
    target = target.join(zs_sum_data)
    return target
    
    
    
