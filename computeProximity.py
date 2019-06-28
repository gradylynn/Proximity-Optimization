# -*- coding: utf-8 -*-
"""
Proximity Optimization function.
@author Ry Arnold
@author Grady Lynn
"""

from collections import deque
import numpy as np
from scipy.spatial import distance


def proximity(xy_dims, featureList, scale=20):
    assert len(xy_dims)==2 and xy_dims[0]>=0 and xy_dims[1]>=0 , "Incorrectly formatted chip dimmensions. Dims must be non-negative x and y dimensions."
    
    ans = np.full(xy_dims, np.inf)
    
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
    
    def isInChip(point):
        return point[0] >= 0 and point[0] < xy_dims[0] and point[1] >= 0 and point[1] < xy_dims[1]
        
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
            if outsidePoint[1] < 0: return left + down
            elif outsidePoint[1] >= xy_dims[1]: return left[::-1] + up
            else: return left
        elif outsidePoint[0] >= xy_dims[0]:
            if outsidePoint[1] < 0: return right + down[::-1]
            elif outsidePoint[1] >= xy_dims[1]: return right[::-1] + up[::-1]
            else: return right
        else:
            if outsidePoint[1] < 0: return down
            else: return up
    
    def solvePoint(feature, point):
        if not isInChip(point): return
        
        d = distance.euclidean(feature, point) * scale
        if d < ans[point[0], point[1]]:
            ans[point[0], point[1]] = d
        else: return
        
        quad = getQuadrant(feature, point)
        
        if quad == 0:
            #Diagonals
            q.append((feature, (point[0] + 1, point[1] + 1))) # Quadrant 2
            q.append((feature, (point[0] - 1, point[1] + 1))) # Quadrant 4
            q.append((feature, (point[0] - 1, point[1] - 1))) # Quadrant 6
            q.append((feature, (point[0] + 1, point[1] - 1))) # Quadrant 8
            #Adjacents
            q.append((feature, (point[0] + 1, point[1]))) # Quadrant 1
            q.append((feature, (point[0], point[1] + 1))) # Quadrant 3
            q.append((feature, (point[0] - 1, point[1]))) # Quadrant 5
            q.append((feature, (point[0], point[1] - 1))) # Quadrant 7
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

        if q[-1] == pushFirst and len(q) > 0:
            q.append(pushSecond)
            return
        else:
            q.append(pushFirst)
            q.append(pushSecond)
            return
  
    q = deque()
    for feature in featureList:
        if isInChip(feature):
            q.append((feature, feature))
        else:
            for edgePoint in getNearestSides(feature):
                q.append((feature, edgePoint))
        
    while len(q) > 0:
        x = q.popleft()
        solvePoint(x[0], x[1])

    return ans


# This don't work as well as the above method...
def computeProximity(xy_dims, featureList, scale=20):
    assert len(xy_dims)==2 and xy_dims[0]>=0 and xy_dims[1]>=0 , "Incorrectly formatted chip dimmensions. Dims must be non-negative x and y dimensions."
    
    ans = np.full(xy_dims, np.inf)
    
    def isInChip(point):
        return point[0] >= 0 and point[0] < xy_dims[0] and point[1] >= 0 and point[1] < xy_dims[1]
    
    def nearestPointInChip(feature):
        if isInChip(feature): return feature
        x = min(max(0, feature[0]), xy_dims[0])
        y = min(max(0, feature[1]), xy_dims[1])
        return (x, y)
    
    def solveQuadrant(feature, deque, quadrant):
        point = deque.popleft()
        if not isInChip(point): return
        
        d = distance.euclidean(feature, point) * scale
        if d < ans[point[0], point[1]]:
            ans[point[0], point[1]] = d
        else: return
        
        if quadrant == 1:
            pushFirst = (point[0] + 1, point[1]) # 
            pushSecond = (point[0], point[1] + 1) # 
        elif quadrant == 2:
            pushFirst = (point[0], point[1] + 1)
            pushSecond = (point[0] - 1, point[1])
        elif quadrant == 3:
            pushFirst = (point[0] - 1, point[1])
            pushSecond = (point[0], point[1] - 1)
        elif quadrant == 4:
            pushFirst = (point[0], point[1] - 1)
            pushSecond = (point[0] + 1, point[1])
            
        if len(deque) > 0 and deque[-1] == pushFirst:
            deque.append(pushSecond)
            return
        else:
            deque.append(pushFirst)
            deque.append(pushSecond)
            return
        
    q = deque()
    for feature in featureList:
        startOnChip = nearestPointInChip(feature)
        ans[startOnChip[0], startOnChip[1]] = distance.euclidean(feature, startOnChip)
        
        q.append((startOnChip[0] + 1, startOnChip[1]))
        while len(q) > 0:
            solveQuadrant(feature, q, 1)
        
        q.append((startOnChip[0], startOnChip[1] + 1))
        while len(q) > 0:
            solveQuadrant(feature, q, 2)
           
        q.append((startOnChip[0] - 1, startOnChip[1]))    
        while len(q) > 0:
            solveQuadrant(feature, q, 3)
            
        q.append((startOnChip[0], startOnChip[1] - 1))    
        while len(q) > 0:
            solveQuadrant(feature, q, 4)
            
    return ans