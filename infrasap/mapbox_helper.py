import json, sys, os, time, argparse, logging
import shapely

import pandas as pd
import numpy as np
import urllib.request as url

### README
### This code has been adpated from the earlier Market Access tools written by Charles Fox.
### Some simple adjustments were made to work with the InfraSAP inputs, but I have not had enough time to re-write this as it's quite complicated.

def CreateODMatrix(input_df, input_df2, lat_name = 'Lat', lon_name = 'Lon', UID = 'ID', 
                    Pop = 'Pop', call_type = 'OSRM', rescue = 0, rescue_num = 0, MB_Toke = '', 
                    sleepTime = 5, osrmHeader = ''):
    '''
    TODO: make the function flexible for MapBox endpoint - add back MB token / different formatting
    '''
#     ffpath = os.path.dirname(infile)
    start = time.time()
    print('\nChosen server: %s\n\nStart time: %s' % (call_type, time.ctime(start)))
#     print('Origins: %s' % infile)
#     print('Destinations: %s\n' % infile_2)

    # Save settings
    save_rate = 5
    def save(returns, j, i, numcalls, rescue_num):
        elapsed_mins = (time.time() - start)/60
        elapsed_secs = (time.time() - start)%60
        total = ((numcalls / float(i)) * (time.time() - start)/60.0)
        remaining = total - elapsed_mins
        print ('\n______________________________________\n')
        print ('\nSave point %s. Running for: %d minutes %d seconds' % (j, elapsed_mins, elapsed_secs))
        print ('\ncalls completed: %d of %d. Est. run time: %d minutes. Time remaining: %d' % (i-1, numcalls, total, remaining))
        print ('\npercentage complete: %d percent' % (((i-1) / float(numcalls)*100)))
        print ('\n______________________________________\n')
        try:
            df = pd.concat(returns)
        except:
            df = returns
        curOutput = os.path.join(ffpath,'temp_file_%d.csv' % rescue_num)
        df.to_csv(curOutput)
        
        
    # Function for calling OSRM server.
    def Call(O_list, D_list, i, O_IDs, D_IDs, header):
        # Convert origins to HTTP request string
        Os = ';'.join(str(coord).replace("'", "").replace(";", "") for coord in O_list)
        # Destinations to HTTP request string
        Ds = ';'.join(str(coord).replace("'", "").replace(";", "") for coord in D_list)
        # Join them together
        data = Os+';'+Ds

        # Define which coords in data string are origins, and which are destinations
        sources = ['%d' % x for x in range(0,len(O_list))]
        sources = ';'.join(str(x).replace("'", "") for x in sources)
        lenth = len(O_list)+len(D_list)
        destinations = ['%d' % x for x in range(len(O_list),lenth)]
        destinations = ';'.join(str(x).replace("'", "") for x in destinations)

        # Build request string
        request = header+data+'?sources='+sources+'&destinations='+destinations+'?&access_token='+MB_Toke
        # Pass request to interweb
        
        try:
            r = url.urlopen(request)
        except:
            print(request)
            time.sleep(5)
            r = url.urlopen(request)
            
        # Error handle
        try:
            # Convert Bytes response to readable Json
            MB_TelTest_json = json.loads(r.read().decode('utf-8'))
            data_block = MB_TelTest_json['durations']
        except:
            data_block = 'null'

        # Build df from JSON
        #sources_label = [str(i['location']) for i in MB_TelTest_json['sources']]
        #dest_label = [str(i['location']) for i in MB_TelTest_json['destinations']]
        sources_label = O_IDs
        dest_label = D_IDs
        chunk = pd.DataFrame(data = data_block, columns = dest_label, index = sources_label)
        # Convert to minutes, stack 2D array to 1D array
        chunk = chunk.stack(level =-1)
        chunk.columns = ['O','D','DIST']
        return(chunk)

    # Generate appropriately split source and destination lists
    def split_and_bundle(in_list,break_size):
        new_list = []
        for i in range (0,(int(max(len(in_list)/break_size,1)))):
            upper = (i+1) * break_size
            lower = (upper - break_size)
            objs = in_list[lower:upper]
            new_list.append(objs)
        if len(in_list) > break_size:
            rem = len(in_list) % break_size
            if rem > 0:
                final = upper+rem
                new_list.append(in_list[upper:final])
        return new_list

    # File Import for sources file
#     input_df = pd.read_csv(infile)
    input_df['source_list'] = input_df[lon_name].map(str).str.cat(input_df[lat_name].map(str), sep = ',')
    input_df['source_list'] = input_df['source_list']+';'
    source_list = input_df['source_list'].values.tolist()
    source_UIDs = input_df[UID].values.tolist()
    #input_df['source_point'] = input_df.apply(lambda x: Point(x[lon_name],x[lat_name]), axis = 1)
    #source_points = input_df['source_point'].tolist()

    # Look to import separate file for destinations; if not, set destinations = sources
#     input_df2 = pd.read_csv(infile_2)
    input_df2['dest_list'] =  input_df2[lon_name].map(str).str.cat(input_df2[lat_name].map(str), sep = ',')
    input_df2['dest_list'] = input_df2['dest_list']+';'
    dest_list = input_df2['dest_list'].values.tolist()
    dest_UIDs = input_df2[UID].values.tolist()           

    if call_type == 'MBT' :
        sources_list = split_and_bundle(source_list, 5)
        dests_list = split_and_bundle(dest_list, 5)
        sources_UIDs = split_and_bundle(source_UIDs, 5)
        dests_UIDs = split_and_bundle(dest_UIDs, 5)
    elif call_type == 'MB'or call_type == 'OSRM':
        sources_list = split_and_bundle(source_list, 12)
        dests_list = split_and_bundle(dest_list, 13)
        sources_UIDs = split_and_bundle(source_UIDs, 12)
        dests_UIDs = split_and_bundle(dest_UIDs, 13)
    else:
        pass
            
    # Run function call across the O-D matrix; output is 'df'
    returns = []
    numcalls = (len(sources_list) * len(dests_list))
    s , d = sources_list, dests_list
    i, j = 1 + (rescue * len(sources_list)), 1 + rescue

    ### Making Calls 
    if call_type == 'Euclid':
        df = EuclidCall(source_list,dest_list,source_points,dest_points)
    else:
        if rescue > 0:
            s = s[rescue:] # possibly rescue -1
            sources_UIDs = sources_UIDs[rescue:]
        print('source list: %s' % len(source_list))
        print('sources list: %s' % len(sources_list))
        print('dest list: %s' % len(dest_list))
        print('dests list: %s' % len(dests_list))
        numcalls_rem = (len(s) * len(d))
        print('\nEstimated remaining calls to chosen server: %d\n' % numcalls_rem)
        print('save points will occur every %d calls\n' % (len(dests_list)))
        if sleepTime > 0:
            time.sleep(sleepTime)
        for O_list in s:
            O_IDs = sources_UIDs[s.index(O_list)]
            for D_list in d:                    
                if sleepTime > 0:
                    time.sleep(sleepTime)
                D_IDs = dests_UIDs[d.index(D_list)]
                if call_type == 'MB':
                    header = 'https://api.mapbox.com/directions-matrix/v1/mapbox/driving/'
                elif call_type == 'MBT':
                    header = 'https://api.mapbox.com/directions-matrix/v1/mapbox/driving-traffic/'
                elif call_type == 'OSRM':
                    header = 'http://router.project-osrm.org/table/v1/driving/'
                    if osrmHeader != '':
                        header = osrmHeader
                try:
                    # prevent server annoyance
                    print('Call to OSRM server number: %d of %s' % (i, numcalls_rem))                
                    returns.append(Call(O_list,D_list,i,O_IDs,D_IDs, header))
                    i += 1
                    j += 1
                except:
                    logging.warning("Error Processing OSRM for i:%s and j:%s" % (i, j))
                    save(returns, j, i, numcalls, rescue_num)
        try:
            df = pd.concat(returns)
        except:
            df = returns

    # re-attach the population of origins and destinations, prep dataframe
    all_matrices = []
    if rescue_num > 0:
        for r in range(0,rescue_num):
            rescued_matrix = pd.read_csv(os.path.join(ffpath,'temp_file_%d.csv' % (r)),header=None)
            rescued_matrix.columns = ['O_UID','D_UID','DIST']
            all_matrices.append(rescued_matrix)
    df = df.reset_index()
    df.columns = ['O_UID','D_UID','DIST']
    all_matrices.append(df)
    new = pd.concat(all_matrices)
    new = new.set_index('O_UID')
    new['DIST'] = new['DIST'].apply(pd.to_numeric)
    popdf = input_df[[UID,Pop]].set_index(UID)
    new['O_POP'] = popdf[Pop]
    new = new.reset_index()
    new = new.set_index('D_UID')
    if dest_list == source_list:
        new['D_POP'] = popdf[Pop]
        new = new.reset_index()
    else:
        popdf_dest = input_df2[[UID,Pop]].set_index(UID)
        new['D_POP'] = popdf_dest[Pop]
        new = new.reset_index()
    new['O_UID'] = new['O_UID'].astype(str)
    new['D_UID'] = new['D_UID'].astype(str)
    new['combo'] = new['O_UID']+'_X_'+new['D_UID']
    new = new.drop_duplicates('combo')
    new = new.drop(['combo'], axis = 1)

    return new