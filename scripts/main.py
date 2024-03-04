import fiona
import shapely
from math import sqrt,inf
import numpy as np
from shapely.geometry import shape, mapping, LineString

with fiona.open("data/meandry_leva.shp") as levy_breh:
    linestring_left=shape(levy_breh[0]['geometry'])
        
        
with fiona.open("data/meandry_prava.shp") as pravy_breh:
    linestring_right=shape(pravy_breh[0]['geometry'])
    
seg_linestring_left=shapely.segmentize(linestring_left,inf)
seg_linestring_right=shapely.segmentize(linestring_right,inf)

array_left=np.array(seg_linestring_left.coords)
array_right=np.array(seg_linestring_right.coords)

class NearestNeighborsCouple (object):
    def __init__(self, p_point, query_point, distance):
        self.p=p_point
        self.q=query_point
        self.distance=distance
    
    @classmethod
    def get_nn_couple (cls, p_points, q_point):
        current_min_dist = inf
        for p in p_points:
            dist = np.linalg.norm(p - q_point)
            if dist<current_min_dist:
                current_min_dist=dist
                nn_couple=NearestNeighborsCouple(p,q_point,dist)
        return nn_couple

def get_nn_list(p_points, q_points):
    nn=[]
    for point in q_points:
        nnc=NearestNeighborsCouple.get_nn_couple(p_points,point)
        nn.append(nnc)
    return nn

nn_lr=get_nn_list(array_right,array_left)
nn_rl=get_nn_list(array_left,array_right)
pricky=np.array([])
print (len(nn_lr))
print (len(nn_rl))
for item_lr in nn_lr:

    for item_rl in nn_rl:
        #print (item_rl.p)
        if all(item_rl.p==item_lr.q) and all(item_rl.q==item_lr.p):
            pricka=LineString([item_lr.p,item_lr.q])
            #print (pricka)
            pricky=np.append(pricky, pricka)
print (len(pricky))

schema = {
    'geometry': 'LineString',
    'properties': {'id': 'int'},
}
# Write a new Shapefile
with fiona.open('data/vysledek.shp', 'w', 'ESRI Shapefile', schema) as c:
    ## If there are multiple geometries, put the "for" loop here
    x=0
    for item in pricky:
        c.write({
            'geometry': mapping(item),
            'properties': {'id': x},
        })
        x+=1