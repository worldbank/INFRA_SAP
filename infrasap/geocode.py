################################################################################
# Reverse Geocode Shapefile
# Benjamin Stewart, July 2018
# Purpose: Get the name of each feature in a feature class
################################################################################

import sys, os, argparse, logging
import requests, xlrd, unicodedata, time, json
import urllib
import urllib.request
import urllib.parse
import geopandas as gpd 

#This bit makes the function emulate a firefox request...I really don't know if it is necessary
userAgent = "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.2.8) Gecko/20100722 Firefox/3.6.8 GTB7.1 (.NET CLR 3.5.30729)" 
headers = { 'User-Agent' : userAgent }
values = { 's' : 'nothing' }
data = urllib.parse.urlencode(values)

def geocodeDF(inDF, location_cols, username):
    ''' Geocode each name in a pandas data frame
    
    INPUT
        inDF [pandas dataframe]
        location_cols [list of ints] - columns to combine (separated by commas) to create geolocation
    '''
    
    for idx, row in inDF.iterrows():
        location = ','.join(row.iloc[location_cols])
        print(location)


def getJSONresponse(url, data, headers):
    ''' Get JSON response from URL request, return JSON
    '''
    urlReq = urllib.request.urlopen(url).read()
    return(json.loads(urlReq.decode('utf-8')))    

def getLocation (location, data=data, headers=headers):
    ''' perform geocode of the location string
    
    INPUT
        location [string] - actual string to pass to geocoding
        [optional] data [urllib.parse.urlencode] - send to url request to emulate a firefox call ... may not be necessary
        [optional] headers [dictionary] - similar to data, used for emulation
    '''
    location2 = urllib.parse.quote(location) # get the internet version of the location string
    #Perform the geonames.org geolocation    
    geoNamesUrl = "http://api.geonames.org/searchJSON?q=" + location2 + "&maxRows=3&username=petie_stewart&featureClass=P"    
    jsonResponse = getJSONresponse(geoNamesUrl, data, headers)
    ret = dict(lat=-999, lng=-999)
    # If the geonames result has a length
    geoCodeSource = "nothing"
    if (jsonResponse['totalResultsCount'] > 0):
        res = jsonResponse['geonames'][0]
        ret = dict(lat=res['lat'], lng=res['lng'])
        geoCodeSource = "geonames"
    else:
        googleUrl = "http://maps.googleapis.com/maps/api/geocode/json?address=" + location2 + "&sensor=false"
        jsonResponse = getJSONresponse(googleUrl, data, headers)
        if (jsonResponse['status'] == "OK"):
            res = jsonResponse['results'][0]['geometry']['location']
            ret = dict(lat=res['lat'], lng=res['lng'])
            geoCodeSource = "google"
        else:
            print(location + " was not found")
    return({'location':ret, 'source':geoCodeSource})


def getGeogLocation (lat, lng, username):
    ''' Perform geonames reverse geocode
    
    INPUT
        lat/lng [number] - latitude and longitude for search
        username [string] - geonames username for authentication
    '''
    #Perform the geonames.org geolocation    
    geoNamesUrl = "http://api.geonames.org/findNearbyPlaceNameJSON?lat=%s&lng=%s&username=%s" % (lat, lng, username)  
    urlReq = requests.get(geoNamesUrl)
    return(urlReq.json())

def getGeonamesResults(j):
    '''This function extracts results from the GeoNames JSON results
        There could be some error checking here, but it hasn't been necessary so far
    '''
    try:
        res = [
                j['geonames'][0]['name'],
                j['geonames'][0]['fclName'],
                j['geonames'][0]['distance']
                ]
    except:
        res = ['','',-1]
    return res
    
def reverseGeocode(inShapefile, username, outFile=''):
    ''' Search for the name of each feature in a supplied shapefile
    
    INPUT
    inShapefile [string or geopandas geodataframe] - path to input shapefile
    '''
    if isinstance(inShapefile, str):
        inD = gpd.read_file(inShp)
    else:
        inD = inShapefile
    #Re-project to WGS84
    if not inD.crs == {'init': u'epsg:4326'}:
        inD = inD.to_crs({'init': 'epsg:4326'})
    
    inD['Name'] = ''
    inD['NameType'] = ''
    inD['dist'] = 0
    for feat in inD.iterrows():
        lat = feat[1].geometry.centroid.y
        lng = feat[1].geometry.centroid.x
        ###Perform Reverse Geocoding
        geocodeRes = getGeonamesResults(getGeogLocation (lat, lng, username))     
        inD.Name[feat[0]] = geocodeRes[0]
        inD.NameType[feat[0]] = geocodeRes[1]
        inD.dist[feat[0]] = geocodeRes[2]
    if outFile != '':
        ###Write Results
        inD.to_file(outFile)
    return(inD)

if __name__ == '__main__':
    # python reverseGeocode_GeoNames.py -i C:\Temp\temp.shp -o C:\Temp\temp_geocoded.shp -u USERNAME
    #Read in commandline arguments
    parser = argparse.ArgumentParser(description='Determine the name of the nearest named place for each feature in a shapefile')

    parser.add_argument('-i', '--input', dest='IN_SHAPE', action='store', help='Input shape file')
    parser.add_argument('-o', '--output', dest='OUTPUT', action='store', help='Output shapefile')    
    parser.add_argument('-u', '--username', dest='USERNAME', action='store', help='username', default=False)
    
    args = parser.parse_args()
    
    input = args.IN_SHAPE
    output = args.OUTPUT
    username = args.USERNAME
    
    reverseGeocode(input, username, outFile=output)