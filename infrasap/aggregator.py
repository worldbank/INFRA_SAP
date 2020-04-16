'''
The following module contains a number of functions to aggregate geospatial outputs into tables for InfraSAP analytics.
'''
import geopandas as gpd
import pandas as pd
import rasterio as rio
from rasterio import features
from rasterstats import zonal_stats
import numpy as np
from shapely.wkt import loads

def rasterize_gdf(inD, field, template, outFile=None, nodata=np.nan):
    ''' Convert geopandas geo data frame to raster of equal size/res to template raster
    
    INPUT
    inD [ geopandas data frame / path ]
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
    if type(target) == str:
        target = gpd.read_file(target)
    
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

def calculate_access_percentages(OD, target, dest_type, rural=False, urban_extents=None thresholds=[0,30,60,120,180,240,300,360,2000]):
    
    if type(target) == str:
        target = gpd.read_file(target).reset_index(drop=True)
    
    if type(OD) == str:
        OD = pd.read_csv(OD, header=[0,1], index_col=0)
        
    origins_geom = OD['origin'][['geometry']].copy()
    origins_geom.loc[:,'geometry'] = origins_geom['geometry'].apply(lambda x: loads(x))
    origins_geom = gpd.GeoDataFrame(origins_geom, crs = 'EPSG:4326', geometry='geometry')
    
    origins_sj = gpd.sjoin(origins_geom, target[['geometry']], how='left', op='intersects')
    OD.loc[:,('origin','target_idx')] = origins_sj['index_right']
    
    if rural:
        origins_sj2 = gpd.sjoin(origins_geom, urban_extents, how='left', op='intersects')
        OD.loc[:,('origin','rural')] = origins_sj2['index_right'].apply(lambda x: 1 if pd.isna(x) else 0)
        # drop all ourban points, still need to check that this works
        OD = OD.loc[OD[('origin','rural')] == 1].copy()        
        
    min_df = OD['origin'].join(pd.DataFrame(OD[dest_type].min(axis=1).apply(lambda x: (x/60)), columns=["tt_min"]))
    min_df = min_df.loc[~pd.isna(min_df.target_idx)]
    min_df.loc[:,"target_idx"] = min_df.target_idx.astype(int)
    min_df.loc[:,'tt_min_cut'] = pd.cut(min_df.tt_min, bins=thresholds)
    
    summary = min_df.groupby(['target_idx','tt_min_cut'])[['grid_code']].sum().unstack()
    summary_pct = summary.apply(lambda x: x/(summary.sum(axis=1)))
    summary_pct.columns = summary_pct.columns.get_level_values(1)
    results = target.join(summary_pct)
    
    return results