'''
Process OSM files to extract specific data types
'''

import os, sys, time, subprocess, argparse, logging

import geopandas as gpd
import osmnx as ox
import pandas as pd
import networkx as nx

from shapely.geometry import box

# Highway features are reclassified to 4 OSMLR classes for simplification and standardization
#   https://mapzen.com/blog/osmlr-2nd-technical-preview/
OSMLR_Classes = {
"motorway":"OSMLR level 1",
"motorway_link":"OSMLR level 1",
"trunk":"OSMLR level 1",
"trunk_link":"OSMLR level 1",
"primary":"OSMLR level 1",
"primary_link":"OSMLR level 1",

"secondary":"OSMLR level 2",
"secondary_link":"OSMLR level 2",
"tertiary":"OSMLR level 2",
"tertiary_link":"OSMLR level 2",

"unclassified":"OSMLR level 3",
"unclassified_link": "OSMLR level 3",
"residential": "OSMLR level 3",
"residential_link": "OSMLR level 3",

"track": "OSMLR level 4",
"service": "OSMLR level 4"
}

class osmExtraction(object):
    '''
    Download and installation instructions and basic usage can be found here
    https://wiki.openstreetmap.org/wiki/Osmosis    
    '''
    def __init__(self, osmosisCmd = r"C:\WBG\Anaconda\Osmosis\bin\osmosis", tempFile = r"C:\WBG\Anaconda\Osmosis\tempExecution.bat"):
        '''
        osmosisCmd is the path to the file osmosis in the bin folder of the downloaded
            stable version of osmosis (link above)
        '''
        self.osmosisCommand = osmosisCmd
        self.tempFile = tempFile
    
    def extractAmmenities(self, inPbf, outFile, amenityList=["amenity.school"], bounds=[], execute=False):
        ''' Read input osmpbf, extract all buildings and spit out to shapefile
        INPUT
        inPbf [string] - path to input pbf
        outFile [string] - path to output shapefile
        
        OPTIONAL
        amenity [string] - definition of amenity to extract - defaults to schools
        bounds [list of numbers??] - bounds of area to extract, if not defined, covers entire input area
        '''
        if len(amenityList) > 1:
            amenity = ",".join(amenityList)
        else:
            amenity = amenityList[0]
        
        baseCommand = r'{osmCommand} --read-pbf {inPbf} --nkv keyValueList="{amenity}"'.format(
            osmCommand = self.osmosisCommand,
            inPbf = inPbf,             
            amenity=amenity)
        if len(bounds) > 0:
            baseCommand = "{baseCommand} --bounding-box top={top} left={left} bottom={bottom} right={right}".format(
                baseCommand=baseCommand,
                top    = bounds[3],
                left   = bounds[0],
                bottom = bounds[1],
                right  = bounds[2])
        
        baseCommand = '{baseCommand}  --write-pbf {outPbf}'.format(
            baseCommand = baseCommand,
            outPbf = outFile
            )
        if not execute:
            return(baseCommand)
        else:
            with open(self.tempFile, 'w') as outFile:
                outFile.write(baseCommand)
            subprocess.call([self.tempFile], shell=True)
    
    def extractBuildings(self, inPbf, outFile, bounds=[], execute=True):
        ''' Read input osmpbf, extract all buildings and spit out to shapefile
        INPUT
        inPbf [string] - path to input pbf
        outFile [string] - path to output shapefile
        '''
        baseCommand = r"{osmCommand} --read-pbf {inPbf} --tf accept-ways building=*".format(
            osmCommand = self.osmosisCommand,
            inPbf = inPbf, 
            outPbf = outFile)
        if len(bounds) > 0:
            baseCommand = "{baseCommand} --bounding-box top={top} left={left} bottom={bottom} right={right}".format(
                baseCommand=baseCommand,
                top    = bounds[3],
                left   = bounds[0],
                bottom = bounds[1],
                right  = bounds[2])
        baseCommand = "{baseCommand}  --write-pbf {outPbf}".format(
                        baseCommand = baseCommand,
                        outPbf = outFile)
        if not execute:
            return(baseCommand)
        else:
            with open(self.tempFile, 'w') as outFile:
                outFile.write(baseCommand)
            subprocess.call([self.tempFile], shell=True)
    
    def extractBoundingBox(self, inPbf, inShp, outPbf, execute=True):
        ''' Extract the inPbf based on the bounding box of the input shapefile
        INPUT
        inPbf [string] - path to input pbf
        inShp [string] - path to aoi shapefile
        '''
        inD = gpd.read_file(inShp)
        
        baseCommand = r"{osmCommand} --read-pbf {inPbf} --bounding-box top={top} left={left} bottom={bottom} right={right} --write-pbf {outPbf}".format(
            osmCommand = self.osmosisCommand,
            inPbf = inPbf,
            top =  float(inD.total_bounds[3]),
            left = float(inD.total_bounds[0]),
            bottom = float(inD.total_bounds[1]),
            right = float(inD.total_bounds[2]), 
            outPbf = outPbf)
        if not execute:
            return(baseCommand)
        else:
            with open(self.tempFile, 'w') as outFile:
                outFile.write(baseCommand)
            subprocess.call([self.tempFile], shell=True)
            
    
    def extractHighways(self, inPbf, outOSM, values=[1,2,3,4], bounds = [], execute=True):
        ''' Extract highways (all roads) from input osm pbf. Can limit definition of which roads to extract
            by defining the OSMLT class. see osmMisc.OSMLR_Classes for definition of classes
            
        INPUT
            inPbf [string]  - path to input osm pbf file
            outOSM [string] - path to output osm pbf file
            values [list of int] [optional][default = [1,2,3,4] - OSMLR values to extract
            bounds [list of coordinates] - boundary to extract
            execute [boolean] [default = True] - if set to false, command will return the osmosis command, and not execute it
        '''
        highwayVals = []
        for key, value in OSMLR_Classes.items():
            try:
                if int(value.split(" ")[2]) in values:
                    #highwayVals.append("highway=%s" % key)
                    highwayVals.append(key)
            except:
                pass
        allCommands = ",".join(highwayVals)
        baseCommand = r"{osmCmd} --read-pbf {inPbf} --tf accept-ways highway={highwayCommand} --used-node".format(
            osmCmd = self.osmosisCommand, 
            inPbf=inPbf, 
            highwayCommand=allCommands, 
            outPbf=outOSM)
        
        if len(bounds) > 0:
            baseCommand = "{baseCommand} --bounding-box top={top} left={left} bottom={bottom} right={right}".format(
                baseCommand=baseCommand,
                top    = bounds[3],
                left   = bounds[0],
                bottom = bounds[1],
                right  = bounds[2])
                
        baseCommand = "{baseCommand} --write-pbf {outPbf}".format(
                        baseCommand=baseCommand,
                        outPbf=outOSM)
        if not execute:
            return(baseCommand)
        else:
            with open(self.tempFile, 'w') as outFile:
                outFile.write(baseCommand)
            subprocess.call([self.tempFile], shell=True)
