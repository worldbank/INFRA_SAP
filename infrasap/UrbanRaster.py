#-------------------------------------------------------------------------------
# Calculate urban areas from gridded population data
# Benjamin P Stewart, April 2019
#   Purpose is to create high density urban clusters and urban cluster above minimum
#   density and total population thresholds
#-------------------------------------------------------------------------------

import os, sys, logging, geojson, json, time

import rasterio
import geopandas as gpd
import pandas as pd
import numpy as np

from scipy import stats
from scipy.ndimage import generic_filter
from scipy.sparse.csgraph import connected_components
from rasterio import features
from rasterio.features import rasterize
from shapely.geometry import shape, Polygon

'''prints the time along with the message'''
def tPrint(s):
    print("%s\t%s" % (time.strftime("%H:%M:%S"), s))
    
def calc_areagrid(inRaster):
    """
    # determine area of cells in unprojected raster
    # from https://gis.stackexchange.com/questions/317392/determine-area-of-cell-in-raster-qgis
    Create urban definitions using gridded population data.
        
    :param inRaster: string or rasterio object representing gridded population data 
    :returns: A Numpy ndarray of the same dimensions of the input raster where output units are in square km.
    """

    import numpy.matlib

    if type(inRaster) == str:
        input_raster = rasterio.open(inRaster)
    elif isinstance(inRaster, rasterio.DatasetReader):
        input_raster = inRaster
        # print(input_raster.closed)
        if input_raster.closed == True:
            print('dataset is closed, try to re-open raster')
    else:
        raise(ValueError("Input raster dataset must be a file path or a rasterio object"))

    with input_raster as testif:
        rast = testif.read(1)
        gt = testif.transform
        pix_width = gt[0]
        ulX = gt[2]
        ulY = gt[5]
        rows = testif.height
        cols = testif.width
        lrX = ulX + gt[0] * cols
        lrY = ulY + gt[4] * rows

    lats = np.linspace(ulY,lrY,rows+1)

    a = 6378137
    b = 6356752.3142

    # Degrees to radians
    lats = lats * np.pi/180

    # Intermediate vars
    e = np.sqrt(1-(b/a)**2)
    sinlats = np.sin(lats)
    zm = 1 - e * sinlats
    zp = 1 + e * sinlats

    # Distance between meridians
    #        q = np.diff(longs)/360
    q = pix_width/360

    # Compute areas for each latitude in square km
    areas_to_equator = np.pi * b**2 * ((2*np.arctanh(e*sinlats) / (2*e) + sinlats / (zp*zm))) / 10**6
    areas_between_lats = np.diff(areas_to_equator)
    areas_cells = np.abs(areas_between_lats) * q

    areagrid = np.transpose(np.matlib.repmat(areas_cells,cols,1))

    return areagrid

class urbanGriddedPop(object):
    def __init__(self, inRaster):
        """
        Create urban definitions using gridded population data.
        
        :param inRaster: string or rasterio object representing gridded population data        
        """
        if type(inRaster) == str:
            self.inR = rasterio.open(inRaster)
        elif isinstance(inRaster, rasterio.DatasetReader):
            self.inR = inRaster
        else:
            raise(ValueError("Input raster dataset must be a file path or a rasterio object"))
            
    def calculateUrban(self, densVal=300, totalPopThresh=5000, smooth=False, verbose=False, queen=False,
                        raster='', raster_pop='', print_message=''):
        """
        Generate urban extents from gridded population data through the application of a minimum
            density threshold and a minimum total population threshold
            
        :param densVal: integer of the minimum density value to be counted as urban
        :param totalPopThresh: integer minimum total settlement population to ne considered urban
        :param smooth: boolean to run a single modal smoothing function (this should be run when running 
                        on WorldPop as the increased resolution often leads to small holes and funny shapes
        :param verbose: boolean on what messages to receive
        :param queen: boolean to determine whether to dissolve final shape to connect queen's contiguity
        :param raster: string path to create a boolean raster of urban and not. 
                        Empty string is the default and will create no raster
        :param raster_pop: string path to create a raster of the population layer only in the urban areas
                            Empty string is the default and will create no raster
        :returns: GeoPandasDataFrame of the urban extents
        """

        popRaster = self.inR
        data = popRaster.read()
        urbanData = (data > densVal) * 1
        urbanData = urbanData.astype('int16')
        if verbose:
            tPrint("%s: Read in urban data" % print_message)
        # Modal filter
        def modal(P):
            mode = stats.mode(P)
            return mode.mode[0]
        '''
        if smooth:
            # Run modal filter
            urbanData[0,:,:] = generic_filter(urbanData[0,:,:], modal, (3, 3))
            tPrint("Smoothed urban data")
        '''
        allFeatures = []
        badFeatures = []
        idx = 0     
        # create output array to store urban raster
        urban_raster = urbanData * 0
        for cShape, value in features.shapes(urbanData, transform=popRaster.transform):
            if idx % 1000 == 0 and verbose:
                tPrint("%s: Creating Shape %s" % (print_message, idx))
            if value == 1:            
                if smooth:
                    xx = shape(cShape)
                    xx = Polygon(xx.exterior)
                    cShape = xx.__geo_interface__
                #If the shape is urban, claculate total pop        
                mask = rasterize([(cShape, 0)], out_shape=data[0,:,:].shape,fill=1,transform=popRaster.transform)
                inData = np.ma.array(data=data, mask=mask.astype(bool))
                curPop = inData.sum() 
                if curPop < 0: # when smoothed, sometimes the pop withh be < 0 because of no data
                    inData = np.ma.array(data=inData, mask=(inData < 0).astype(bool))
                    curPop = inData.sum() 
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
        
        if queen:
            # res returns the (width, height) of pixels in the units of its coordinate reference system
            xxGeom['geometry '] = xxGeom.buffer((popRaster.res[0] / 2))
            s = xxGeom['geometry']
            overlap_matrix = s.apply(lambda x: s.intersects(x)).values.astype(int)
            n, ids = connected_components(overlap_matrix)
            xxGeom['group'] = ids
            xxGeom = xxGeom.dissolve(by="group", aggfunc="sum")
        
        return(xxGeom)
                
        