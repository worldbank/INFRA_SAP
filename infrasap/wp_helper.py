import sys, os
import rasterio, geojson, json

import numpy as np
import pandas as pd
import geopandas as gpd

from affine import Affine
from rasterio.warp import reproject, Resampling
from scipy.ndimage import generic_filter
from scipy import stats
from rasterio import features
from rasterio.features import rasterize
from shapely.geometry import shape


def resample_wp(in_wp, outfile, factor=10):
    ''' resample the input rasterio object and write to output
    
    INPUT
        in_wp [rasterio] - original 100m resolution wp dataset
        outfile [string path to output file]
        factor [int] - ratio to upscale. Default is 10, meaning res of 100m will become 1000m
    '''
    out_meta = in_wp.meta.copy()
    out_meta.update(width=round(out_meta['width']/factor), height=round(out_meta['height']/factor))
    a = out_meta['transform']
    new_affine = Affine(a[0] * factor, a[1], a[2], a[3], a[4] * factor, a[5])
    out_meta.update(transform=new_affine)
    new_arr = np.empty(shape=(1, out_meta['height'], out_meta['width']))
    reproject(in_wp.read(masked=True), new_arr, src_transform=in_wp.meta['transform'], dst_transform=new_affine, 
                src_crs = in_wp.crs, dst_crs = in_wp.crs, resampling = Resampling.average)
    new_arr = new_arr * (factor * factor)
    new_arr[new_arr < 0] = out_meta['nodata']
    with rasterio.open(outfile, 'w', **out_meta) as out_file:
        out_file.write(new_arr.astype(out_meta['dtype']))            
        
        
def calculateUrban(popRaster, densVal=300, totalPopThresh=5000, smooth=False, verbose=False,
                            raster='', raster_pop=''):
        """
        Generate urban extents from gridded population data through the application of a minimum
            density threshold and a minimum total population threshold
            
        :param popRaster: rasterio object of gridded population data
        :param densVal: integer of the minimum density value to be counted as urban
        :param totalPopThresh: integer minimum total settlement population to ne considered urban
        :param smooth: boolean to run a single modal smoothing function (this should be run when running 
                        on WorldPop as the increased resolution often leads to small holes and funny shapes
        :param verbose: boolean on what messages to receive
        :param raster: string path to create a boolean raster of urban and not. 
                        Empty string is the default and will create no raster
        :param raster_pop: string path to create a raster of the population layer only in the urban areas
                            Empty string is the default and will create no raster
        :returns: GeoPandasDataFrame of the urban extents
        """

        data = popRaster.read()
        urbanData = (data > densVal) * 1
        urbanData = urbanData.astype('int16')
        if verbose:
            print("Read in urban data")
        # Modal filter
        def modal(P):
            mode = stats.mode(P)
            return mode.mode[0]

        if smooth:
            # Run modal filter
            urbanData[0,:,:] = generic_filter(urbanData[0,:,:], modal, (3, 3))
            print("Smoothed urban data")
            
        allFeatures = []
        badFeatures = []
        idx = 0     
        verbose=True
        totalPopThresh = 5000
        # create output array to store urban raster
        urban_raster = urbanData * 0
        for cShape, value in features.shapes(urbanData, transform=popRaster.transform):
            if value == 1:            
                if idx % 100 == 0 and verbose:
                    print("Creating Shape %s" % idx)
                #If the shape is urban, claculate total pop        
                mask = rasterize([(cShape, 0)], out_shape=data[0,:,:].shape,fill=1,transform=popRaster.transform)
                if idx % 100 == 0 and verbose:
                    print("Rasterized feature")
                inData = np.ma.array(data=data, mask=mask.astype(bool))
                curPop = inData.sum() 
                if curPop < 0: # when smoothed, sometimes the pop withh be < 0 because of no data
                    inData = np.ma.array(data=inData, mask=(inData < 0).astype(bool))
                    curPop = inData.sum() 
                if idx % 100 == 0 and verbose:
                    print("Summed Pop")
                if curPop > totalPopThresh:            
                    allFeatures.append([idx, curPop, shape(geojson.loads(json.dumps(cShape)))])
                    urban_raster += (mask^1)
                else:
                    badFeatures.append([idx, curPop, shape(geojson.loads(json.dumps(cShape)))])
            idx = idx + 1
        
        if len(raster):
            out_metadata = popRaster.meta.copy()
            out_metadata['dtype'] = urban_raster.dtype
            out_metadata['nodata'] = 0
            with rasterio.open(raster, 'w', **out_metadata) as rOut:
                rOut.write(urban_raster)
        
        if len(raster_pop):
            out_metadata = popRaster.meta.copy()
            urban_pop = data * urban_raster
            with rasterio.open(raster_pop, 'w', **out_metadata) as rOut:
                rOut.write(urban_pop)
        
        xx = pd.DataFrame(allFeatures, columns=['ID', 'Pop','geometry'])
        xxGeom = gpd.GeoDataFrame(xx, geometry='geometry')
        xxGeom.crs = popRaster.crs
        return(xxGeom)