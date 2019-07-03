# -*- coding: utf-8 -*-
"""
Proximity Optimization function.
@author Ry Arnold
@author Grady Lynn
"""

from collections import deque
import numpy as np
from scipy.spatial import distance
from heapq import heappush, heappop
import math

def proximity(xy_dims, featureList, scale=20):
    assert len(xy_dims)==2 and xy_dims[0]>=0 and xy_dims[1]>=0 , "Incorrectly formatted chip dimmensions. Dims must be non-negative x and y dimensions."
    
    def isInChip(point):
        return point[0] >= 0 and point[0] < xy_dims[0] and point[1] >= 0 and point[1] < xy_dims[1]
    
    ans = np.full(xy_dims, np.inf)
    global countChanged, c
#    visited = np.full(xy_dims, 0)
    
    def getQuadrant(origin, point):
        if point == origin: return 0
        if point[0] > origin[0] and point[1] >= origin[1] and point[0] - origin[0] > point[1] - origin[1]: return 1
        if point[0] > origin[0] and point[1] > origin[1] and point[0] - origin[0] <= point[1] - origin[1]: return 2
        if point[0] <= origin[0] and point[1] > origin[1] and -(point[0] - origin[0]) < point[1] - origin[1]: return 3
        if point[0] < origin[0] and point[1] > origin[1] and -(point[0] - origin[0]) >= point[1] - origin[1]: return 4
        if point[0] < origin[0] and point[1] <= origin[1] and point[0] - origin[0] < point[1] - origin[1]: return 5
        if point[0] < origin[0] and point[1] < origin[1] and point[0] - origin[0] >= point[1] - origin[1]: return 6
        if point[0] >= origin[0] and point[1] < origin[1] and -(point[0] - origin[0]) > point[1] - origin[1]: return 7
        if point[0] > origin[0] and point[1] < origin[1] and -(point[0] - origin[0]) <= point[1] - origin[1]: return 8 
        
    up = []
    down = []
    right = []
    left = []
    for x in range(xy_dims[0]):
        up.append((x, xy_dims[1] - 1))
        down.append((x, 0))
    for y in range(xy_dims[1]):
        right.append((xy_dims[0] - 1, y))
        left.append((0, y))
    
    def getNearestSides(outsidePoint):
        if outsidePoint[0] < 0:
            if outsidePoint[1] < 0: return left + down[1:]
            elif outsidePoint[1] >= xy_dims[1]: return left + up[1:]
            else: return left
        elif outsidePoint[0] >= xy_dims[0]:
            if outsidePoint[1] < 0: return right[1:] + down
            elif outsidePoint[1] >= xy_dims[1]: return right[:-1] + up
            else: return right
        else:
            if outsidePoint[1] < 0: return down
            else: return up
    
    def solvePoint(x):
        global countChanged, c
        
        feature = x[1][0]
        point = x[1][1]
        d = x[0]
        
        c += 1
        
        if not isInChip(point): return
        
#        visited[point[0], point[1]] += 1
        
        if d < ans[point[0], point[1]]:
            ans[point[0], point[1]] = d
            countChanged += 1
        else: return
        
        quad = getQuadrant(feature, point)
        
        if quad == 0:
            #Adjacents
            heappush(q, (scale, (feature, (point[0] + 1, point[1]))) ) # Quadrant 1
            heappush(q, (scale, (feature, (point[0], point[1] + 1)))) # Quadrant 3
            heappush(q, (scale, (feature, (point[0] - 1, point[1])))) # Quadrant 5
            heappush(q, (scale, (feature, (point[0], point[1] - 1)))) # Quadrant 7
            #Diagonals
            heappush(q, (math.sqrt(2)*scale, (feature, (point[0] + 1, point[1] + 1)))) # Quadrant 2
            heappush(q, (math.sqrt(2)*scale, (feature, (point[0] - 1, point[1] + 1)))) # Quadrant 4
            heappush(q, (math.sqrt(2)*scale, (feature, (point[0] - 1, point[1] - 1)))) # Quadrant 6
            heappush(q, (math.sqrt(2)*scale, (feature, (point[0] + 1, point[1] - 1)))) # Quadrant 8
            return
        elif quad == 1:
            pushFirst = (feature, (point[0] + 1, point[1]))
            pushSecond = (feature, (point[0] + 1, point[1] + 1))
        elif quad == 2:
            pushFirst = (feature, (point[0], point[1] + 1))
            pushSecond = (feature, (point[0] + 1, point[1] + 1))
        elif quad == 3:
            pushFirst = (feature, (point[0], point[1] + 1))
            pushSecond = (feature, (point[0] - 1, point[1] + 1))
        elif quad == 4:
            pushFirst = (feature, (point[0] - 1, point[1]))
            pushSecond = (feature, (point[0] - 1, point[1] + 1))
        elif quad == 5:
            pushFirst = (feature, (point[0] - 1, point[1]))
            pushSecond = (feature, (point[0] - 1, point[1] - 1))
        elif quad == 6:
            pushFirst = (feature, (point[0], point[1] - 1))
            pushSecond = (feature, (point[0] - 1, point[1] - 1))
        elif quad == 7:
            pushFirst = (feature, (point[0], point[1] - 1))
            pushSecond = (feature, (point[0] + 1, point[1] - 1))
        elif quad == 8:
            pushFirst = (feature, (point[0] + 1, point[1]))
            pushSecond = (feature, (point[0] + 1, point[1] - 1))

        heappush(q, (distance.euclidean(pushFirst[0], pushFirst[1]) * scale, pushFirst))
        heappush(q, (distance.euclidean(pushSecond[0], pushSecond[1]) * scale, pushSecond))
        return
        
    q = []
    countChanged = 0
    c = 0
    for feature in featureList:
        if isInChip(feature):
            heappush(q, (0, (feature, feature)) )
        else:
            for edgePoint in getNearestSides(feature):
                heappush(q, (distance.euclidean(feature, edgePoint)*scale, (feature, edgePoint)))
        
    while countChanged < xy_dims[0] * xy_dims[1]:
        x = heappop(q)
        solvePoint(x)

#    print(visited)
    print(c)
    return ans


# This don't work as well as the above method...
#def computeProximity(xy_dims, featureList, scale=20):
#    assert len(xy_dims)==2 and xy_dims[0]>=0 and xy_dims[1]>=0 , "Incorrectly formatted chip dimmensions. Dims must be non-negative x and y dimensions."
#    
#    ans = np.full(xy_dims, np.inf)
#    
#    def isInChip(point):
#        return point[0] >= 0 and point[0] < xy_dims[0] and point[1] >= 0 and point[1] < xy_dims[1]
#    
#    def nearestPointInChip(feature):
#        if isInChip(feature): return feature
#        x = min(max(0, feature[0]), xy_dims[0])
#        y = min(max(0, feature[1]), xy_dims[1])
#        return (x, y)
#    
#    def solveQuadrant(feature, deque, quadrant):
#        point = deque.popleft()
#        if not isInChip(point): return
#        
#        d = distance.euclidean(feature, point) * scale
#        if d < ans[point[0], point[1]]:
#            ans[point[0], point[1]] = d
#        else: return
#        
#        if quadrant == 1:
#            pushFirst = (point[0] + 1, point[1]) # 
#            pushSecond = (point[0], point[1] + 1) # 
#        elif quadrant == 2:
#            pushFirst = (point[0], point[1] + 1)
#            pushSecond = (point[0] - 1, point[1])
#        elif quadrant == 3:
#            pushFirst = (point[0] - 1, point[1])
#            pushSecond = (point[0], point[1] - 1)
#        elif quadrant == 4:
#            pushFirst = (point[0], point[1] - 1)
#            pushSecond = (point[0] + 1, point[1])
#            
#        if len(deque) > 0 and deque[-1] == pushFirst:
#            deque.append(pushSecond)
#            return
#        else:
#            deque.append(pushFirst)
#            deque.append(pushSecond)
#            return
#        
#    q = deque()
#    for feature in featureList:
#        startOnChip = nearestPointInChip(feature)
#        ans[startOnChip[0], startOnChip[1]] = distance.euclidean(feature, startOnChip)
#        
#        q.append((startOnChip[0] + 1, startOnChip[1]))
#        while len(q) > 0:
#            solveQuadrant(feature, q, 1)
#        
#        q.append((startOnChip[0], startOnChip[1] + 1))
#        while len(q) > 0:
#            solveQuadrant(feature, q, 2)
#           
#        q.append((startOnChip[0] - 1, startOnChip[1]))    
#        while len(q) > 0:
#            solveQuadrant(feature, q, 3)
#            
#        q.append((startOnChip[0], startOnChip[1] - 1))    
#        while len(q) > 0:
#            solveQuadrant(feature, q, 4)
#            
#    return ans