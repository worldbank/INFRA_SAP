'''
Process OSM files to extract specific data types
'''

import os, sys, time, subprocess, argparse, logging
import osmium

import geopandas as gpd
import osmnx as ox
import pandas as pd
import networkx as nx
import shapely.wkb as wkblib

from shapely.geometry import box, LineString, Point
from shapely.ops import transform
from functools import partial

wkbfab = osmium.geom.WKBFactory()


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

class highwayExtractor(osmium.SimpleHandler):
    ''' Extract highways with relevant tags
    '''
    def __init__(self, verbose=False):
        osmium.SimpleHandler.__init__(self)
        self.OSMLR1 = []
        self.OSMLR2 = []
        self.OSMLR3 = []
        self.OSMLR4 = []
        self.OSMLROther = []
        self.tags = ['service','railway',"usage","gauge","electrified"]
    
    def way(self, n):
        try:
            highway = n.tags.get("highway")
            try:
                osmlr_level = OSMLR_Classes[highway]
            except:
                osmlr_level = ""
            wkb = wkbfab.create_linestring(n)
            shp = wkblib.loads(wkb, hex=True)
            res = [n.id, shp, highway]
            for t in self.tags:
                try:
                    res.append(n.tags.get(t))
                except:
                    res.append("")
            
            if osmlr_level == "OSMLR level 1":
                self.OSMLR1.append(res)
            if osmlr_level == "OSMLR level 2":
                self.OSMLR2.append(res)
            if osmlr_level == "OSMLR level 3":
                self.OSMLR3.append(res)
            if osmlr_level == "OSMLR level 4":
                self.OSMLR4.append(res)
            if osmlr_level == "other":
                self.OSMLROther.append(res)
                
        except:
            pass
        
class railwayExtractor(osmium.SimpleHandler):
    ''' Extract railways with relevant tags
    '''
    def __init__(self, verbose=False):
        osmium.SimpleHandler.__init__(self)
        self.railways = []
        self.tags = ['service','railway',"usage","gauge","electrified"]
    
    def way(self, n):
        if n.tags.get("railway"):
            try:
                wkb = wkbfab.create_linestring(n)
                shp = wkblib.loads(wkb, hex=True)
                res = [n.id, shp]
                for t in self.tags:
                    try:
                        res.append(n.tags.get(t))
                    except:
                        res.append("")
                self.railways.append(res)
            except:
                pass
        


# extract directly through osmium
class InfraExtractor(osmium.SimpleHandler):
    """ Extractor for use in osmium SimpleHandler to extract nodes and highways
    """
    def __init__(self, verbose=False):
        ''' Extract nodes representing ports and international airports
        '''
        osmium.SimpleHandler.__init__(self)
        self.verbose = verbose
        self.ports = []
        self.airports = []  
        self.railways = []
        self.rail_stations = []
        self.bridges = []
        self.highways = []
        
    def node(self, n):
        if n.tags.get('aeroway') == 'aerodrome':
            wkb = wkbfab.create_point(n)
            shp = wkblib.loads(wkb, hex=True)
            self.airports.append([n.id, shp, shp.x, shp.y, n.tags.get('aerodrome:type'), n.tags.get('name'), n.tags.get('name:en')])
        if n.tags.get('harbour') == "yes":
            wkb = wkbfab.create_point(n)
            shp = wkblib.loads(wkb, hex=True)
            self.ports.append([n.id, shp, shp.x, shp.y])
        if n.tags.get('railway') == "station":
            wkb = wkbfab.create_point(n)
            shp = wkblib.loads(wkb, hex=True)
            self.rail_stations.append([n.id, shp, n.tags.get("name")])
        if n.tags.get('bridge'):
            wkb = wkbfab.create_point(n)
            shp = wkblib.loads(wkb, hex=True)
            self.bridges.append([n.id, shp.x, shp.y, n.tags.get("name")])
    
    def way(self, n):
        if n.tags.get("railway"):
            try:
                wkb = wkbfab.create_linestring(n)
                shp = wkblib.loads(wkb, hex=True)
                self.railways.append([n.id, shp, n.tags.get("name")])
            except:
                pass
        if n.tags.get('aeroway') == 'aerodrome':
            try:
                wkb = wkbfab.create_multipolygon(n)
            except:
                try:
                    wkb = wkbfab.create_linestring(n)
                except:
                    print(n)
                    wkb = None
            if not wkb is None:
                shp = wkblib.loads(wkb, hex=True)
                self.airports.append([n.id, shp, shp.centroid.x, shp.centroid.y, n.tags.get('aerodrome:type'), n.tags.get('name'), n.tags.get('name:en')])
        if n.tags.get('harbour') == "yes":
            try:
                wkb = wkbfab.create_multipolygon(n)
            except:
                wkb = wkbfab.create_linestring(n)
            shp = wkblib.loads(wkb, hex=True)
            self.ports.append([n.id, shp, shp.centroid.x, shp.centroid.y])
        if n.tags.get("highway"):
            try:
                nodes = [x.ref for x in n.nodes]
                wkb = wkbfab.create_linestring(n)
                shp = wkblib.loads(wkb, hex=True)
                info = [n.id, nodes, shp, n.tags['highway']]
                self.highways.append(info)
            except:
                nodes = [x for x in n.nodes if x.location.valid()]
                if len(nodes) > 1:
                    shp = LineString([Point(x.location.x, x.location.y) for x in nodes])
                    info = [n.id, nodes, shp, n.tags['highway']]
                    self.highways.append(info)
                
def check_international(x):
    try:
        if x['TYPE'].upper() == "International".upper():
            return("International")
        elif "International" in x['Name']:
            return("International")
        elif "International" in x['NameEN']:
            return("International")
        else:
            return(x['TYPE'])
    except:
        return(x['TYPE'])
    
def load_pois(osm_file, exact_bounds):
    h = InfraExtractor()
    h.apply_file(osm_file, locations=True)
    
    airports = pd.DataFrame(h.airports, columns=["OSM_ID", 'geometry', 'x', 'y', "TYPE", "Name", "NameEN"])
    airports_geom = [Point(x) for x in zip(airports['x'], airports['y'])]
    airports.drop(['geometry'], axis=1, inplace=True)
    airports = gpd.GeoDataFrame(airports, geometry=airports_geom, crs={'init':'epsg:4326'})
    airports = airports[airports.intersects(exact_bounds)]
    try:
        airports['TYPE2'] = airports.apply(lambda x: check_international(x), axis=1)
    except:
        airports['TYPE2'] = ""
    
    ports =    pd.DataFrame(h.ports, columns=["OSM_ID",'geometry', 'x', 'y'])
    ports_geom = [Point(x) for x in zip(ports['x'], ports['y'])]
    ports.drop(['geometry'], axis=1, inplace=True)    
    ports =    gpd.GeoDataFrame(ports, geometry=ports_geom, crs={'init':'epsg:4326'})
    ports =    ports[ports.intersects(exact_bounds)]
    
    rails = pd.DataFrame(h.railways, columns=["ID","geometry", "name"])
    rails = gpd.GeoDataFrame(rails, geometry="geometry", crs={'init':'epsg:4326'})
    rails = rails[rails.intersects(exact_bounds)]
    
    stations = pd.DataFrame(h.rail_stations, columns=['ID','geometry','name'])
    stations = gpd.GeoDataFrame(stations, geometry="geometry", crs={'init':'epsg:4326'})
    stations = stations[stations.intersects(exact_bounds)]
    
    highways = pd.DataFrame(h.highways, columns=["ID","nodes","geometry","type"])
    highways = gpd.GeoDataFrame(highways, geometry="geometry", crs={'init':'epsg:4326'})
    highways.drop(['nodes'], axis=1, inplace=True)
    highways['OSMLR'] = highways['type'].map(OSMLR_Classes)
    
    return({
        'airports':airports,
        'ports':ports,
        'railways':rails,
        'rail_stations':stations,
        'highways':highways
    })

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
