import math,random

# modeled from: https://github.com/gboeing/osmnx/blob/89b5771c6e060dee28fc0f56c3e7ee4c43b0be26/osmnx/projection.py
def project_gdf(gdf):
    # calculate longitude of centroid of union of all geometries in gdf
    avg_lng = gdf["geometry"].unary_union.centroid.x

    # calculate UTM zone from avg longitude to define CRS to project to
    utm_zone = math.floor((avg_lng + 180) / 6) + 1
    
    # testing
    print(f"print utm zone: {utm_zone}")
    
    utm_crs = f"+proj=utm +zone={utm_zone} +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

    # project the GeoDataFrame to the UTM CRS
    gdf_proj = gdf.to_crs(utm_crs)
    
    return gdf_proj
    
#-----------------------------------------------------------------------------
# TEST IF POINT IS IN POLYGON (HELPER ALGORITHM)

# To test if a point is in a shape, a vertical line is drawn from the point to
# the lowest extent. The number of intersections between the vertical line and
# the shape perimeter is determined. If the number of intersections is more than
# zero and is odd, then the test point is in the shape.


# requires list of feature XY coordinates
# requires a coordinates for a test point

def pntInShp(pntLst,testPnt):   
    
    # XY coordinates of the test point...   
    Xr, Yr = testPnt[0], testPnt[1]
    
    ringNo = 0
    
    for ringList in pntLst:

        intersect = 0    # number of intersections
        inPoly = 0
        
        # for each point in feature part...
        for i in range(len(ringList)-1):

            # get X coordinates of segment end points...
            X1,X2 = ringList[i][0], ringList[i+1][0]

            # skip if test point is to the left or right of segment...
            if Xr < X1 and Xr < X2: continue
            if Xr > X1 and Xr > X2: continue
         
            # get Y coordinates of segment end points...
            Y1,Y2 = ringList[i][1], ringList[i+1][1]

            # skip if test point is below segment...
            if Yr < Y1 and Yr < Y2:
                continue

            # if X values are the same, skip...
            if X1 == X2: continue

            # if perimeter segement is vertical, then line must be on it...
            if Y1 == Y2:
                intersect += 1   # count intersection
                continue

            # slope of line between pt1 and pt2...
            m = (Y2-Y1) / (X2-X1)

            # y value on line that has X value equal to test point X coordinate...
            y = m*(Xr-X1)+Y1

            # if test point Y is greater than y on the line between pt1 and pt2...
            if Yr >= y:
                # then vertical line containing test point intersects polygon boundary...
                intersect += 1

        # if more than 1 intersection and total number is odd,
        if intersect > 0 and intersect % 2 != 0:
            # then test point is inside polygon...
            inPoly = 1

        ##        print ringNo,inPoly

        if ringNo == 0 and inPoly == 0: break

        # means that if the point was in the exterior polygon and also inside one of the interior polygons
        # then it is not inside of the polygon 
        if ringNo > 0 and inPoly == 1:
            inPoly = 0
            break

        ringNo += 1

    else:
        inPoly = 1
       
    return inPoly   # 1 if test point is in shape, 0 if not in shape
    
#------------------------------------------------------------------
# IDENTIFY GRID CELLS IN ROW THAT BOUND ROW PIXELS IN FEATURE...

def bounds(studyArea_pntLst,Y):

    # list of intersections between horizontal line and feature...
    intersectLst = []

    for ringLst in studyArea_pntLst:
        
        #print(f"print ringLst: {ringLst}")

        # for each position in feature point list...
        for pos in range(len(ringLst)-1):
            
            # XY coordinates of perimeter segment endpoints...
            X1, Y1 = ringLst[pos][0], ringLst[pos][1]
            X2, Y2 = ringLst[pos+1][0], ringLst[pos+1][1]  
            
            #print(f"print X1, Y1: {X1, Y1}")
            #print(f"print X2, Y2: {X2, Y2}")
            
            #---------------------------------------
            # equation for perimeter segment...

            # if vertical perimeter segment
            if X1 == X2:
                # if test Y not in range of perimeter segment, there is no intersection
                if not(Y1 <= Y <= Y2 or Y2 <= Y <= Y1):
                    continue
                # if Y in range, set x value of intercept
                x = X1

            # if perimeter segment not vertical
            else:
                # perimeter segment slope
                m = (Y2-Y1)/(X2-X1)      
                #print(f"perimeter segment slope: {m}")

                # intercept of segment
                B1 = Y1-(m*X1)
                #print(f"intercept of segment: {B1}")

                # if slope equals zero, there can be no intersection...
                if m == 0:
                    continue

                # calculate x of intercept
                x = (Y - B1) / m
                
                #print(f"print x of intercept: {x}")

                # skip if intersection not on perimeter segment...
                if not(X1 <= x <= X2 or X2 <= x <= X1): continue

            # add intersection X coordinate to list...
            intersectLst.append(x)
        
        
    #print(f"print intersectLst: {intersectLst}")

    # sort list in ascending order...
    intersectLst.sort()

    #print(f"print intersectLst: {intersectLst}")

    #-----------------------------
    # CREATE FINAL RANGE LIST...

    # For each pair of consecutive X values,
    # check if line segment is contained within the shape...


    # final list of ranges (between each set of X values is the inside of shape)...
    rangeLst = []

    # for each pos in intersection list...
    for pos in range(len(intersectLst)-1):

        # X coordinates...
        X1, X2 = intersectLst[pos], intersectLst[pos+1]

        # X coordinate of segment midpoint...
        Xm = X1 + (X2-X1)/2

        # if segment midpoint is in shape, add range to list...
        if pntInShp(studyArea_pntLst,[Xm,Y]) == 1:
            rangeLst.append([X1,X2])

    return rangeLst
    
#-----------------------------
def generate_featPntLst(pntLst, minX, minY, maxX, maxY, cellsize, gridTxtFile):

    # open new text file...
    o_grid = open(gridTxtFile,'w')

    # add blank lines for header lines (to be added in later)...
    for x in range(100):
        o_grid.write("\n")

    # total row count...
    nrows = 0

    featPntLst = []   # list of coordinates of points inside feature 

    # set current Y to maxY...
    Y = maxY

    # column count...
    ncols = 0

    Ycounter = 0
    Xcounter = 0

    # for each row...
    while Y >= minY:
        
        #print(f"Ycounter is: {Ycounter}")
        Ycounter += 1
        
        #if Ycounter > 10:
            #break
        
        nrows += 1   # count row

        line = ""    # text file line

        # get range of X values within which lie the inside of feature...
        rangeLst = bounds(pntLst,Y)
        
        #print(f"rangeLst is: {rangeLst}")

        # testing
        #break
        
        # set X to leftmost column...
        X = minX
        
        # column count...
        ncols = 0
     
        # for each column..
        while X <= maxX:
            
            #print(f"Xcounter is: {Xcounter}")
            Xcounter += 1
            #if Xcounter > 10:
                #break

            ncols +=1   # count column

            # for each range of X values...
            for x1,x2 in rangeLst:
                # if current X is within range, point is in shape...
                if x1 <= X <= x2:
                    line += "1 "    # "1" indicates feature pixel
                    featPntLst.append([X,Y])   # at point to feature point list
                    break
            # if current X is not within an range of interior X values
            else:
                line += "0 "    # "0" indicates non-feature pixel

            # X coordinate of next column
            X += cellsize

        line += "\n"    # add return to line
        #print(line)
        o_grid.write(line)   # write line to txt file

        # Y coordinate for next row...
        Y -= cellsize

    # close grid...
    o_grid.close()

    #--------------------------
    # Add header lines to grid text file...

    # coordinates of lower left corner of grid (Note: this is the corner
    # of the pixel, not the center)
    ll_Y = maxY - (nrows-1)*cellsize - cellsize/2     # Y coordinate
    ll_X = minX-cellsize/2                            # X coordinate

    # open existing text file in edit mode...
    o_grid = open(gridTxtFile,'r+') 

    # write header lines to grid text file...
    o_grid.write("ncols %s\n" % ncols)
    o_grid.write("nrows %s\n" % nrows)
    o_grid.write("xllcorner %s\n" % ll_X)
    o_grid.write("yllcorner %s\n" % ll_Y)
    o_grid.write("cellsize %s\n" % cellsize)
    o_grid.write("NODATA_value %s\n" % -9999)
    o_grid.close()
    
    return featPntLst


#----------------------------------------------------------------------
# PROXIMITY INDEX
def proximity(featPntLst,centerX,centerY,r):
   
    Xc,Yc = centerX, centerY

    inPix = 0 
    sumD = 0   # sum of distances

    # for each feature point...
    count = 0
    for X,Y in featPntLst:
        count = count + 1
        # distance to center...
        d = ((X-Xc)**2 + (Y-Yc)**2)**.5
        
        #if count < 5:
            #print(f"print X: {X} and Y :{Y}")
            #print(f"print distance: {d}")
        
        sumD += d

        # if distance < EAC radius, then pixel is in
        if d < r:
            inPix += 1  # count pixel

    # calculate average distance
    D_to_Center = sumD / len(featPntLst)

    return D_to_Center, inPix
    
#---------------------------------------------------------------------------

## SECTION 3: COHESION FUNCTION

# The average distance between all pairs of points in the shape. Estimated
# using 30 samples of 1000 points...

def interpointDistance(ptList): # requires list of XY coordinates of points in shape

    # number of points in shape...
    numPts = len(ptList)

    samplSize = 1000    # number of points in sample
    samples = 30        # number of samples
        
    avgD = 0    # average interpoint distance

    # run specified number of samples...
    for t in range(samples):

        total_D = 0     # cumulative distance
        cnt = 0         # number of calculations
        
        # select a random sample of shape points...
        sampLst = random.sample(ptList, samplSize)

        # for each point in sample, calculate distance between it and
        # every other point...
        for pt1 in sampLst:
            
            # get coordinates of a point...   
            X1 = pt1[0]
            Y1 = pt1[1]

            # calculate distance to all other points in sample...
            for pt2 in sampLst:

                # skip pt2 if it is the same as pt1...
                if pt1 == pt2: continue 

                # get coord. of point from sample...
                X2 = pt2[0]
                Y2 = pt2[1]

                # calculate distance to pt1...
                dist = ((X1-X2)**2 + (Y1-Y2)**2)**.5
                
                # cumulative interpoint distance...
                total_D += dist

                # number of calculations...
                cnt += 1
                    
        # average interpoint distance...
        avgD += (total_D / cnt) / samples
         
    return avgD
    
## SECTION 9: THE SPIN INDEX (RELATIVE MOMENT OF INERTIA)...
# Change header ZD
# orig: def spin(XYLst,centroidXY):
def spin(XYLst,centroidX, centroidY):

    # XY coordiantes of shape centroid...
    Xc, Yc = centroidX, centroidY

    # sum of distance squared...
    sum_dsqr = 0

    # count of number of pixels...
    cnt = 0

    # for each shape point...
    for pt in XYLst:

        # XY coordiante of point...
        X, Y = pt[0], pt[1]

        # distance to center squared...
        dsqr = (X-Xc)**2+(Y-Yc)**2
        
        # sum of squared distances...
        sum_dsqr += dsqr

        # count of points...
        cnt += 1

    return sum_dsqr / cnt
    
#---------------------------------------------------------------------------
# GET LIST OF POINTS EVENLY SPACED ALONG PERIMETER
# Vertices not preserved...

def PerimeterPnts(coordLst, numPnts):

    perimPntLst_allRings = []

    totalPerim = 0    

    ringPerimLst = []

    for ringLst in coordLst:

        # perimeter length...
        perim = 0

        # for each perimeter segment, add lengths...
        for pos in range(len(ringLst)-1):
            x1,y1 = ringLst[pos][0], ringLst[pos][1]
            x2,y2 = ringLst[pos+1][0], ringLst[pos+1][1]
            d = ((x1-x2)**2 + (y1-y2)**2)**.5
            perim += d

        ringPerimLst.append(perim)
        totalPerim += perim

    spacingLst = []

    for perim in ringPerimLst:
        rNumPnts = int((perim / totalPerim) * numPnts)
        if rNumPnts == 0: 
            d = rNumPnts = 0
        else:
            d = perim / rNumPnts
            if d < .5: d = .5
        spacingLst.append([d, rNumPnts])

    
    # for each ring in polygon...
    for ringNo in range(len(coordLst)):

        ringLst = coordLst[ringNo]
        d = spacingLst[ringNo][0]
        numPnts = spacingLst[ringNo][1]
        
        if numPnts == 0: continue
        
        #---------------------------------
        # GENERATE POINTS ALONG PERIMETER OF CURRENT RING...

        # list of selected points...
        perimPntLst = []

        # position in vertex list...
        pos = done = 0

        # endpoint coordinates of first perimeter segment...
        X1,Y1 = ringLst[pos][0], ringLst[pos][1]
        X2,Y2 = ringLst[pos+1][0], ringLst[pos+1][1]

        # for each point desired...
        #for pntNum in xrange(numPnts+15):
        while True:
                   
            # determine the min and max X values for the segment...
            xLst = [X1,X2]
            xLst.sort()
            xMin = xLst[0]
            xMax = xLst[-1]

            # determine the min and max Y values for the segment...
            yLst = [Y1,Y2]
            yLst.sort()
            yMin = yLst[0]
            yMax = yLst[-1]

            # 1 = valid point found; 0 = no valid point found...
            pnt1 = pnt2 = 0
            
            # if segment slope is nearly vertical...
            if abs(X1 - X2) <= abs(Y1-Y2)/100:
                x1 = x2 = X1    
                y1 = Y1 - d      
                y2 = Y1 + d
                
            # if segment slope is horizontal...
            elif Y2-Y1 == 0:
                y1 = y2 = Y1
                x1 = X1 - d
                x2 = X1 + d

            # if segment is not nearly vertical, calculate slope and y-intercept
            else:
            
                m = (Y2-Y1) / (X2-X1)    # slope
                b = Y1-(m*X1)            # y-intercept
               
                #---------------------------
                # find point on line that is distance d from (X1,Y1)...

                # do not understand this original code:
                # coefficients in the quadratic equation...
                # A = 1+m**2
                # B = 2*m*b - 2*X1 - 2*Y1*m
                # C = X1**2 - d**2 + Y1**2 - 2*Y1*b + b**2

                # calculate x intercepts using quadratic equation...
                # x1 = (-B + (B**2-4*A*C)**.5) / (2*A)
                # y1 = m*x1+b

                # x2 = (-B - (B**2-4*A*C)**.5) / (2*A)
                # y2 = m*x2+b
                
                # finding the point on line that is distance d from (X1,Y1) this way
                d_line = ((X2-X1)**2+(Y2-Y1)**2)**.5
                d_ratio = d/d_line
                
                # line direction of ring can be ordered in either direction depending if it is the outer ring, or an interior ring
                # this is the reason for there being two sets of coordinates
                x1,y1 = ((1-d_ratio)*X1+d_ratio*X2),((1-d_ratio)*Y1+d_ratio*Y2)
                x2,y2 = ((1-d_ratio)*X1-d_ratio*X2),((1-d_ratio)*Y1-d_ratio*Y2)

            # check if 1st point is on the segment...
            # print(f"print x1: {x1}")
            # print(f"print type x1: {type(x1)}")
            # print(f"print xMin: {xMin}")
            # print(f"print type xMin: {type(xMin)}")
            # print(f"print xMax: {xMax}")
            # print(f"print type xMax: {type(xMax)}")
            
            if xMin <= x1 <= xMax:
                if yMin <= y1 <= yMax:
                    # check if point is a vertex...
                    if ((x1-X2)**2 + (y1-Y2)**2)**.5 < .00001:
                        pos += 1
                        # if position is last vertex, analysis is finished...
                        if pos >= len(ringLst)-1:
                            break
                        X2,Y2 = ringLst[pos+1][0], ringLst[pos+1][1]
                    pnt1 = 1

            # check if 2nd point is on the segment...
            if xMin <= x2 <= xMax:
                if yMin <= y2 <= yMax:
                    # check if point is a vertex...
                    if ((x2-X2)**2 + (y2-Y2)**2)**.5 < .00001:
                        pos += 1
                        # if position is last vertex, analysis is finished...
                        if pos >= len(ringLst)-1:
                            break
                        X2,Y2 = ringLst[pos+1][0], ringLst[pos+1][1] 
                    pnt2 = 1
            
                    
            #---------------------------
            
            dNext = d  # additional distance needed (set to full separation distance initially)...

            # if neither point is on line segment, move along successive segments...        
            while pnt1 == pnt2 == 0:

                # calculate additional distance needed to meet spacing requirement...   
                dNext = dNext-((X1-X2)**2 + (Y1-Y2)**2)**.5

                # move position to next vertex (to get next segment)...
                pos += 1

                # if position is last vertex, analysis is finished...
                if pos >= len(ringLst)-1:
                    break

                # if close to vertex, take vertex as the point...
                if dNext < .0001:
                    x1,y1 = X2,Y2
                    X2,Y2 = ringLst[pos+1][0], ringLst[pos+1][1]
                    pnt1 = 1
                    break
            
                # coordinates of next perimeter segment...
                X1,Y1 = ringLst[pos][0], ringLst[pos][1]
                X2,Y2 = ringLst[pos+1][0], ringLst[pos+1][1]

                # determine the min and max X values for the segment...
                xLst = [X1,X2]
                xLst.sort()
                xMin = xLst[0]
                xMax = xLst[-1]

                # determine the min and max Y values for the segment...
                yLst = [Y1,Y2]
                yLst.sort()
                yMin = yLst[0]
                yMax = yLst[-1]

                # 1 = valid point found; 0 = no valid point found...
                pnt1 = pnt2 = 0

                # if segment slope is nearly vertical...
                if abs(X1 - X2) <= abs(Y1-Y2)/100:
                    x1 = x2 = X1
                    y1 = Y1 - dNext
                    y2 = Y1 + dNext

                # if segment is not nearly vertical, calculate slope and y-intercept
                else:                
               
                    m = (Y2-Y1) / (X2-X1)
                    b = Y1-(m*X1)

                    #---------------------------
                    # find point on line that is distance d from (X1,Y1)...

                    # do not understand this original code:
                    # coefficients in the quadratic equation...
                    # A = 1+m**2
                    # B = 2*m*b - 2*X1 - 2*Y1*m
                    # C = X1**2 - dNext**2 + Y1**2 - 2*Y1*b + b**2
                  
                    # calculate x intercepts using quadratic equation...
                    # x1 = (-B + (B**2-4*A*C)**.5) / (2*A)
                    # y1 = m*x1+b

                    # x2 = (-B - (B**2-4*A*C)**.5) / (2*A)
                    # y2 = m*x2+b
                    
                    # finding the point on line that is distance d from (X1,Y1) this way
                    d_line = ((X2-X1)**2+(Y2-Y1)**2)**.5
                    d_ratio = dNext/d_line
                    
                    x1,y1 = ((1-d_ratio)*X1+d_ratio*X2),((1-d_ratio)*Y1+d_ratio*Y2)
                    x2,y2 = ((1-d_ratio)*X1-d_ratio*X2),((1-d_ratio)*Y1-d_ratio*Y2)

                # print(f"print2 x1: {x1}")
                # print(f"print2 type x1: {type(x1)}")
                # print(f"print2 xMin: {xMin}")
                # print(f"print2 type xMin: {type(xMin)}")
                # print(f"print2 xMax: {xMax}")
                # print(f"print2 type xMax: {type(xMax)}")
                
                if not isinstance(x1,complex):
                    # check if 1st point is on the segment...
                    if xMin <= x1 <= xMax:
                        if yMin <= y1 <= yMax:
                            pnt1 = 1  
                            
                if not isinstance(x2,complex):
                    # check if 2nd point is on the segment...
                    if xMin <= x2 <= xMax:
                        if yMin <= y2 <= yMax:
                            pnt2 = 1
                        
            # if position is last vertex, analysis is finished...
            if pos >= len(ringLst)-1:
                break

            # if point 1 is valid, and not already in list, add to list
            if pnt1 == 1:
                if [x1,y1] not in perimPntLst:
                    perimPntLst.append([x1,y1])

            # if point 2 is valid, and not already in list, add to list
            elif pnt2 == 1:
                if [x2,y2] not in perimPntLst:
                    perimPntLst.append([x2,y2])

            # set 1st endpoint to last perimeter point...
            X1,Y1 = perimPntLst[-1][0], perimPntLst[-1][1]

        # add last vertex point to perimeter point list...
        perimPntLst.append(ringLst[-1])

        perimPntLst_allRings.append(perimPntLst)

    return perimPntLst_allRings
    
# SECTION 13: CALCULATE PIXEL DISTANCE TO PERIPHERY (for depth and girth indices)...

def pt_distToEdge(featPntLst,perimPntLst):

    pt_dToE = []    # list for pixel distances to nearest edge pixel
    perimPnts = []

    for ringLst in perimPntLst:
        perimPnts.extend(ringLst)

    # for each pixel in shape...
    for X1,Y1 in featPntLst:

        dLst = []   # list of distances from shape pixel to all edge pixels

        # for each edge pixel
        for X2,Y2 in perimPnts:

            # calculate distance between shape and perimeter...
            d = ((X1-X2)**2 + (Y1-Y2)**2)**.5

            # add distance to list...
            dLst.append(d)

        # sort distance list to get minimum distance...
        dLst.sort()

        # add minimum distance to list...
        pt_dToE.append(dLst[0])

    return pt_dToE
    
## SECTION 15: THE GIRTH INDEX (MAX RADIUS OF INSCRIBED CIRCLE)
def girth(pt_dToE):

    # copy pixel distance to edge list...
    sortLst = pt_dToE[:]

    # sort list to identify maximum distance
    sortLst.sort()
    shp_Girth = sortLst[-1] # max distance to edge...

    # identify position of inmost point in pt_dToE list...
    ptPos = pt_dToE.index(shp_Girth)

    return shp_Girth

## SECTION 14: THE DEPTH INDEX (DISTANCE TO EDGE OF SHAPE)...
def depth(pt_dToE):

    depth = 0
    
    cnt = len(pt_dToE)

    # for each distance in list...
    for d in pt_dToE:
        
        # calculate average of distances...
        depth += d / cnt

    return depth

# requires list with XY of proximate center and a list with XY for perimeter points...
def dispersion(center,perimPntLst):

    sumD = 0
    dList = []

    # XY coordinates of center...
    Xc, Yc = center[0], center[1]

    # for point in shape perimeter...
    for X1,Y1 in perimPntLst:

        # double check if need to change to great-circle distance
        # calculate distance to proximate center...
        d = ((X1-Xc)**2+(Y1-Yc)**2)**.5

        # add distance to list...
        dList.append(d)

        # cumulative distance...
        sumD += d

    # number of points in list...
    numPts = len(dList)

    # average distance from center to perimeter...
    avgD = sumD / numPts

    # calculate std deviation between distances...
    V = 0

    for d in dList:
        V += abs(d - avgD)  # cumulative std deviation (difference between distance and average distance)


    # calculate dispersion index
##    avgVi = V / (2 * numPts)    # std deviation / 2n
    avgVi = V / (numPts)
    dispersionIndex = (avgD - avgVi)/avgD

    return dispersionIndex, avgD
    
# HELPER ALGORITHM (CONVEX HULL): RIGHT-HAND TURN TEST

# Tests if 3 consecutive points form a right hand turn...

def right_turn(pts):

    ## right_turn algorithm borrowed from BoundingContainers.py by Dan Patterson.

    # for three points p, q, r...
    # calculate the determinant of a special matrix with three 2D points.
    # The sign, "-" or "+", determines the side, right or left,
    # respectively, on which the point r lies, when measured against
    # a directed vector from p to q.  Use Sarrus' Rule to calculate
    # the determinant. (could also use the Numeric package...) 
    
    p,q,r = pts[0], pts[1], pts[2]  # for points p, q, and r

    # Sarrus' Rule (for computing determinants)...
    sum1 = q[0]*r[1] + p[0]*q[1] + r[0]*p[1]
    sum2 = q[0]*p[1] + r[0]*q[1] + p[0]*r[1]

    # if sum1 < sum2, then p,q,r form a right-hand turn...
    if sum1 < sum2: rightTurn = 1
    
    # else p,q,r do not form a left-hand turn
    else: rightTurn = 0
    
    return rightTurn
    
# SECTION 17: CREATE CONVEX HULL (FOR USE IN THE DETOUR INDEX)...

    ## The ConvexHull algorithm was borrowed from BoundingContainers.py written by:

        #  Dan Patterson
        #  Dept of Geography and Environmental Studies
        #  Carleton University, Ottawa, Canada
        #  Dan_Patterson@carleton.ca
        #
        # Date  Jan 25 2006

    ## Starting with a point known to be on the hull -- leftmost point -- test
    ## every 3 consecutive points to see if they form right-hand turns. For any set
    ## of 3 consecutive points that does not form a right-hand turn, 
    ## then the middle point is not part of the hull - it is not included in the hull list.
    ## The points of the upper hull are derived in the 1st section, points in the
    ## lower hull are derived in the 2nd section. Both sets of points are combined to form
    ## the complete hull list.

def ConvexHull(vertexLst):

    hullPoints_allParts = []

    pntLst = vertexLst[:]

    pntLst.sort()  # sort vertex points in ascending order

    #------------------------------------------------
    # Build upper half of the hull (leftmost point to rightmost point in hull).

    upperLst = [pntLst[0], pntLst[1]]   # add 1st two points to upperLst

    # for all points except 1st two...
    for p in pntLst[2:]:
        upperLst.append(p)     # add point to upper list

        # check if last 3 points in the lowerLst form a right-hand turn
        # (must be more than 2 points in lowerLst)...
        while len(upperLst) > 2 and right_turn(upperLst[-3:]) == 0:
            # delete the middle of the 3 points if they do not form a right-hand turn...
            del upperLst[-2]
    
    # Build lower half of the hull (rightmost point to leftmost point in hull).
    pntLst.reverse()  # sort pntLst in descending order
    lowerLst = [pntLst[0], pntLst[1]]
    
    for p in pntLst[2:]:
        lowerLst.append(p)  # add point to lower list

        # check if last 3 points in the lowerLst form a right-hand turn
        # (must be more than 2 points in lowerLst)...
        while len(lowerLst) > 2 and right_turn(lowerLst[-3:]) == 0:
            # delete the middle of the 3 points if they do not form a right-hand turn... 
            del lowerLst[-2]  
    
    # remove duplicates...
    del lowerLst[0]
    del lowerLst[-1]

    # combine upperLst and lowerLst...
    hullPoints = upperLst+lowerLst

    #-----------------------------------------
    # get total perimeter of convex hull...

    hullPerim = 0    # length of hull perimeter

    # for each point in hull...
    for pos in range(len(hullPoints)):

        # starting point of segment...
        X1, Y1 = hullPoints[pos][0], hullPoints[pos][1]

        # endpoint of segment
        if pos+1 < len(hullPoints)-1:
            X2, Y2 = hullPoints[pos+1][0], hullPoints[pos+1][1]   # next point in list
        else:
            X2, Y2 = hullPoints[0][0], hullPoints[0][1]   # 1st point in list (completes final segment)

        # add distance of current hull segment
        hullPerim += ((X1-X2)**2 + (Y1-Y2)**2)**.5
    
    return hullPerim
    
# SECTION 18: THE RANGE INDEX (RADIUS OF MINIMUM CIRCUMSCRIBED CIRCLE)...

# for every pair of points in shape, calculate distance - this
# distance is the minimum diameter of a circle that circumscribes the shape

# use the points in the exterior polygon shape
def Range(pntLst):

    dMax = 0    # maximum distance between 2 pts on the perimeter

    # for a point in the convex hull...
    for X1,Y1 in pntLst:

        # for another point in the convex hull...
        for X2,Y2 in pntLst:
     
            # distance between point 1 and point 2
            d = ((X1-X2)**2+(Y1-Y2)**2)**.5

            # if distance is greater than max distance
            if d > dMax:
                # set maximum distance
                dMax = d
                # keep XY coordinates of both points
                finalPtLst = [[X1,Y1],[X2,Y2]]

    # return maximum distance between 2 points
    return dMax
    
# DETERMINE IF A GIVEN LINE IS FULLY CONTAINED WITHIN A SHAPE...

def lineInPoly(line1,perimPntLst):

    poly_pntLst = []

    for ringLst in perimPntLst:
        poly_pntLst.extend(ringLst)


    # coordinates of line endpoints...
    X1,Y1 = float(line1[0][0]),float(line1[0][1])
    X2,Y2 = float(line1[1][0]),float(line1[1][1])

    # number of vertices in shape...
    numVertices = len(poly_pntLst)

    # initialize variables...
    intersect_cnt = lineInPoly = 0

    # for each position in vertex list...
    for pos in range(len(poly_pntLst)-1):

        # vertex coordinates...            
        x1,y1 = poly_pntLst[pos][0], poly_pntLst[pos][1]

        # if pos is not last item in list, get coordinates of next point...
        if pos < len(poly_pntLst)-2:
            x2,y2 = poly_pntLst[pos+1][0], poly_pntLst[pos+1][1]

        # if pos is last item in list, get coordinates of 1st point (final line segment)
        else:
            x2,y2 = poly_pntLst[0][0], poly_pntLst[0][1]

        
        # if both lines are the same, then test line is the same as the perimeter line segment...
        # and is inside the shape...
        if [x1,y1] in [[X1,Y1],[X2,Y2]] and [x2,y2] in [[X1,Y1],[X2,Y2]]:
            lineInPoly = 1          # no intersection
            break
        
        #---------------------------------------
        # equation for 1st line...

        # slope of line segment...
        if X1 == X2: M1 = 99999999999   # vertical line
        else: M1 = (Y2-Y1)/(X2-X1)      # calculate slope

        # intercept of segment
        B1 = Y1-(M1*X1)

        #---------------------------------------
        # equation of line 2...

        # slope of line segment
        if x1 == x2: m2 = 99999999999
        else: m2 = (y2-y1)/(x2-x1)

        # intercept of line segment...
        b2 = y1-m2*x1

        #---------------------------------------
        # test for intersection...

        # if equal slopes, lines will never intersect so skip...
        
        if M1 != m2:       
            # calculate XY of intersection...                 
            X = (b2-B1)/(M1-m2)
            Y = M1*X + B1

            # distance between intersection and endpoint...
            D1 = ((X-X1)**2+(Y-Y1)**2)**.5      # end pnt 1
            D2 = ((X-X2)**2+(Y-Y2)**2)**.5      # end pnt 2
           
            # if almost no distance, intersection is an endpoint, do not count...
            if D1 < 1 or D2 < 1:
                continue

            # if intersection is not an endpoint, check if it is on the line segment...
            else:

                # sort X coordinates for each segment in ascending order...
                XLst,xLst = [X1,X2],[x1,x2]
                XLst.sort()
                xLst.sort()

                # if intersection X is in X range of both lines, then count it...
                if XLst[0] <= X <= XLst[1] and xLst[0] <= X <= xLst[1]:              
                    intersect_cnt += 1

    # if at least one intersection found, line is not fully contained within the shape...
    if intersect_cnt > 0:
        lineInPoly = 0
    # if no intersections found, make sure line is within shape (test line midpoint)...
    else:

        # calculate coordinates of test line midpoint...
        Xm = X1 + (X2-X1)/2
        Ym = Y1 + (Y2-Y1)/2

        # if line mid is in the shape, then entire line must be in shape...
        if pntInShp(perimPntLst,[Xm,Ym]) == 1:
            lineInPoly = 1
        
    # value = 0 indicates line not fully contained within shape
    # value = 1 indicates line fully contained within shape
    return lineInPoly
    
