###############################################################################
# Summarize OSM roads lengths
# Charles Fox, July 2018
# Purpose: summarize road lengths within features in defined shapefile
###############################################################################

import os, sys, time, subprocess, argparse, logging

import GOSTnets as gn
#from . import OSMNX_POIs

import geopandas as gpd
import osmnx as ox
import pandas as pd
import networkx as nx

from shapely.geometry import box
#import misc


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
        inShp [string or geopandas object] - path to aoi shapefile
        '''
        if type(inShp) == str:
            inD = gpd.read_file(inShp)
        elif type(inShp) == gpd.GeoDataFrame:
            inD = inShp
        else:
            raise(ValueError("inShp needs to be a string or a geopandas object"))
        if inD.crs != {'init':'epsg:4326'}:
            inD = inD.to_crs({'init':'epsg:4326'})
        baseCommand = r"{osmCommand} --read-pbf {inPbf} --bounding-box top={top} left={left} bottom={bottom} right={right} --write-pbf {outPbf}".format(
            osmCommand = self.osmosisCommand,
            inPbf = inPbf,
            top =  float(inD.bounds.maxy),
            left = float(inD.bounds.minx),
            bottom = float(inD.bounds.miny),
            right = float(inD.bounds.maxx), 
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
        

#grid = gpd.read_file(r"Q:\AFRICA\COD\Projects\Shohei_Poverty_Kinshasa\ADMIN\PSUs\bati_ilot_quartier.shp")
#outFolder = r"Q:\AFRICA\COD\Projects\Shohei_Poverty_Kinshasa\ADMIN"
def downloadBaseData(grid, outFolder, amenities=True):
    '''Download OSM based data using OSMNX - roads, schools, hospitals, and churches.
    
    INPUT
    grid [geopandas dataFrame] - area to download in
    outFolder [string] - place to write output
    RETURNS
    dictionary of [geopandas]
    '''
    toReturn = {}
    roadsFile = os.path.join(outFolder, "OSM_Roads.shp")
    if not os.path.exists(roadsFile):
        bbox = box(grid.bounds.minx.min(), grid.bounds.miny.min(), grid.bounds.maxx.max(), grid.bounds.maxy.max())
        #Download road network
        G = ox.graph_from_polygon(bbox, network_type='drive_service')
        roads = gn.edge_gdf_from_graph(G)
        roads['highway'] = roads.highway.astype(str)
        roads['OSMLR'] = roads.highway.map(OSMLR_Classes)
        roads['oneway'] = roads.oneway.astype(int)   
        '''
        for badKeys in ['access', 'bridge','junction', 'lanes','oneway', 'osmid', 'ref', 'service','tunnel','width','stnode','endnode','name']:
            try:
                roads = roads.drop([badKeys],axis=1)
            except:
                print("Could not drop %s" % badKeys)
        '''
        try:        
            roads.to_file(roadsFile)
        except:
            print("Could not write output")
    else:
        roads = pd.read_file
    toReturn['Roads'] = roads
        
    if amenities:    
        #Download Schools
        schools = OSMNX_POIs.AmenityObject('Health', bbox, ['clinic','pharmacy','hospital','health'], "C:/Temp")
        schools = schools.GenerateOSMPOIs()
        schools = schools.RemoveDupes(0.005, roads.crs)
        schools = schools.prepForMA()
        schools.to_csv(os.path.join(outFolder, "OSM_Health.csv"), encoding='utf-8')
        toReturn['Schools'] = schools
        
        #Download Hospitals
        health = OSMNX_POIs.AmenityObject('Education', bbox, ['school','university','secondary school', 'kindergarten', 'college'], "C:/Temp")
        health = health.GenerateOSMPOIs()
        health = health.RemoveDupes(0.005, roads.crs)
        health = health.prepForMA()
        health.to_csv(os.path.join(outFolder, "OSM_Schools.csv"), encoding='utf-8')  
        toReturn['Health'] = health
        
        #Download Churches
        placeOfWorship = OSMNX_POIs.AmenityObject('Churches', bbox, ['place_of_worship'], "C:/Temp")
        placeOfWorship = placeOfWorship.GenerateOSMPOIs()
        placeOfWorship = placeOfWorship.RemoveDupes(0.005, roads.crs)
        placeOfWorship = placeOfWorship.prepForMA()
        placeOfWorship.to_csv(os.path.join(outFolder, "OSM_Churches.csv"), encoding='utf-8') 
        toReturn['placeOfWorship'] = placeOfWorship
        
    return(toReturn)
    
def summarizeOSM(grid, verbose=True, roadsOnly=False):
    ''' Summarizes OSM road length within each feature in the input grid
    
    ---variables---
    grid [GeoDataframe] - each feature will have the length of all roads summarized

    --- To Do ---
    1. The length projection is web mercator - this is not great, and should instead be projected to UTM
    '''
    WGS_84 = {'init' :'epsg:4326'}
    WEB_MERCATOR = {'init' :'epsg:3857'}
    #Extract OSM within the bounding box of the input grid
    if grid.crs != WGS_84:
        grid = grid.to_crs(WGS_84)
    bbox = box(grid.bounds.minx.min(), grid.bounds.miny.min(), grid.bounds.maxx.max(), grid.bounds.maxy.max())
    G = ox.graph_from_polygon(bbox, network_type='drive_service')
    #Limit the OSM grid to important columns, ignore nodes
    nodes,roads = ox.save_load.graph_to_gdfs(G, edges = True)
    if roadsOnly:
        return(roads)
    roads = roads[['geometry','highway','osmid']]
    roads['length'] = roads.length
    roads['highway'] = roads.highway.astype(str)
    roads['OSMLR'] = roads.highway.map(OSMLR_Classes)
    FID_list, OSMLR1_list, OSMLR2_list, OSMLR3_list, OSMLRt_list = [], [], [], [], []
    cnt = 0
    verbose = False
    #Loop through all of the input features to summarize OSM
    for idx, obj in grid.iterrows():
        if idx % 50 == 0 and verbose:
            print ("%s of %s" % (cnt, len(grid.shape[0])))
        roads2 = roads.copy()
        #Find roads that intersect current geometry
        roads2 = roads2[roads2.intersects(obj.geometry)]
        #Find the intersection of the current geometry and the intersecting roads
        roads2['intersecting'] = roads2.geometry.intersection(obj.geometry)    
        roads2 = roads2.set_geometry('intersecting')
        roads2['intersecting_length'] = roads2.length
        #Project the roads to a metre-based CRS
        roads2 = roads2.to_crs(WEB_MERCATOR)
        roads2['intersecting_length'] = roads2.length
        FID_list.append(idx)
        #Summarize total lenght of OSMLR classes
        OSMLR1_list.append(roads2.loc[roads2.OSMLR == 'OSMLR level 1'].intersecting_length.sum())
        OSMLR2_list.append(roads2.loc[roads2.OSMLR == 'OSMLR level 2'].intersecting_length.sum())
        OSMLR3_list.append(roads2.loc[roads2.OSMLR == 'OSMLR level 3'].intersecting_length.sum())
        OSMLRt_list.append(roads2.loc[roads2.OSMLR == 'OSMLR track'].intersecting_length.sum())
        cnt = cnt + 1

    results = pd.DataFrame({'FID':FID_list,
                       'OSMLR level 1': OSMLR1_list,
                       'OSMLR level 2': OSMLR2_list,
                       'OSMLR level 3': OSMLR3_list,
                       'OSMLR track': OSMLRt_list})
    results['totalOSM'] = results['OSMLR level 1'] + results['OSMLR level 2'] + results['OSMLR level 3'] + results['OSMLR track']
    return results

def convertOSMPBF_DataFrame(inOSM, layer):
    ''' Convert an input OSM PBF to a geopandas dataframe
    
    INPUT
    inOSM [string] - path to OSM pbf to convert
    layer [string] - data layer to extract. Select from lines, points, multipolygons, multilinestrings
    
    RETURNS
    [geopandas data frame]
    '''
    import ogr, geojson, json
    from shapely.geometry import shape
    
    driver=ogr.GetDriverByName('OSM')
    data = driver.Open(inOSM)
    layer = data.GetLayer(layer)
    features=[x for x in layer]
    def loadOSMjson(x):
        ''' convert the OSM JSON to something concumable as geopandas
        '''
        geom = shape(geojson.loads(json.dumps(x['geometry'])))
        x['geometry'] = geom
        for key, value in x['properties'].items():
            x[key] = value
        try:
            splitTags = x['other_tags'].split(",")
            for tag in splitTags:
                tSplit = tag.split("=>")
                x[tSplit[0].replace('"', '')] = tSplit[1].replace('"', '')
        except:
            pass
        x.pop("properties", None)
        x.pop("other_tags", None)
        return(x)
    
    allFeats = [loadOSMjson(feat.ExportToJson(as_object=True)) for feat in features]
    inR = pd.DataFrame(allFeats)
    inR = gpd.GeoDataFrame(inR, geometry='geometry')
    inR.crs = {'init': 'epsg:4326'}
    try:
        inR = misc.project_UTM(inR)
    except:
        inR = inR.to_crs({'init':'epsg:3857'})
    inR['Length'] = inR['geometry'].apply(lambda x: x.length)
    inR = inR.to_crs({'init':'epsg:4326'})
    return (inR)
    
if __name__ == "__main__":
    exampleText =  '''
    #Extract OSMLR level 1 and 2 roads from pbf
    python osmMisc.py -i Q:\AFRICA\MRT\INFRA\mauritania-latest_20190103.osm.pbf -o Q:\AFRICA\MRT\INFRA\mauritania-latest_20190103_OSMLR12.osm.pbf -OSMLR 1 2 -OSMExtract
    
    #Extract OSM PBF for bounding box
    python osmMisc.py -i Q:\AFRICA\LBR\INFRA\liberia-latest.osm.pbf -o Q:\AFRICA\LBR\INFRA\liberia-latest_Monrovia.osm.pbf -AOI_file "Q:\AFRICA\LBR\Projects\Monrovia_Imagery\Greater Monrovia admin boundary\Greater_Monrovia_admin_boundary.shp" -BBExtract
    
    #Extract buildings from OSM PBF
    python osmMisc.py -i Q:\AFRICA\LBR\INFRA\liberia-latest_Monrovia.osm.pbf -o Q:\AFRICA\LBR\INFRA\liberia-latest_Monrovia_Buildings.osm.pbf -BuildingExtract
    '''
    parser = argparse.ArgumentParser(description="Generate urban metrics based on a defined grid")
    
    parser.add_argument('-i',       dest="INPUT", action='store', help="Input .osm.pbf")
    parser.add_argument('-o',       dest="OUTPUT", action='store', help="Output .osm.pbf")    
    parser.add_argument('-OSMExtract',  dest="EXTRACT", action='store_true', help="Perform OSMExtract on OSM PBF")    
    parser.add_argument('-OSMLR',   dest='LEVELS', action='store', help="OSMLR levels to extract from the input pbf", nargs='+')
    parser.add_argument('-BBExtract',  dest="BBEXTRACT", action='store_true', help="Perform OSMExtract on OSM PBF")       
    parser.add_argument('-BuildingExtract',  dest="BUILDINGEXTRACT", action='store_true', help="Perform OSMExtract on OSM PBF")       
    parser.add_argument('-AOI_file', dest="AOI", action='store', help="AOI of study area")
    
    args = parser.parse_args()
    print(args)
    #Set logging information
    logging.basicConfig(format='%(asctime)s:%(levelname)s:  %(message)s', level=logging.INFO)
    xx = osmExtraction()
    if args.EXTRACT:
        y = xx.extractOSMLRlvl1(args.INPUT, args.OUTPUT, values=[int(x) for x in args.LEVELS], execute=True)
        print(y)
    if args.BBEXTRACT:
        y = xx.extractBoundingBox(args.INPUT, args.AOI, args.OUTPUT)
        print(y)    
    if args.BUILDINGEXTRACT:
        y = xx.extractBuildings(args.INPUT, args.OUTPUT)
        print(y)       