# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from collections import deque
import numpy as np
from scipy.spatial import distance

def computeProximity(xy_dims, featureList, scale=20):
    assert len(xy_dims)==2 and xy_dims[0]>=0 and xy_dims[1]>=0 , "Incorrectly formatted chip dimmensions. Dims must be non-negative x and y dimensions."
    
    ans = np.full(xy_dims, np.inf)
    
    def isInChip(point):
        return point[0] >= 0 and point[0] < xy_dims[0] and point[1] >= 0 and point[1] < xy_dims[1]
    
    def nearestPointInChip(feature):
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
            pushFirst = (point[0] + 1, point[1])
            pushSecond = (point[0], point[1] + 1)
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
        if isInChip(feature): q.append(feature)
        else: q.append(nearestPointInChip(feature))
        
        while len(q) > 0:
            solveQuadrant(feature, q, 1)
        
        q.append((feature[0] - 1, feature[1]))
        while len(q) > 0:
            solveQuadrant(feature, q, 2)
            
        q.append((feature[0], feature[1] - 1))
        while len(q) > 0:
            solveQuadrant(feature, q, 3)
            
        q.append((feature[0] + 1, feature[1] - 1))
        while len(q) > 0:
            solveQuadrant(feature, q, 4)
            
    return ans