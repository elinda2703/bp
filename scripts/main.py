import fiona
import shapely
from math import sqrt,inf
import numpy as np
from shapely.geometry import shape, mapping, LineString

with fiona.open("data/vzorecek_leva.shp") as levy_breh:
    linestring_left=shape(levy_breh[0]['geometry'])
        
        
with fiona.open("data/vzorecek_prava.shp") as pravy_breh:
    linestring_right=shape(pravy_breh[0]['geometry'])
    
seg_linestring_left=shapely.segmentize(linestring_left,25)
seg_linestring_right=shapely.segmentize(linestring_right,25)

array_left=np.array(seg_linestring_left.coords)
array_left=np.flipud(array_left)
print(len(array_left))
array_right=np.array(seg_linestring_right.coords)
print(len(array_right))

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

"""
def get_splits
def delete_used
def delete_unpairable



"""


#pricky=[]


def get_all_pricky(array_left,array_right,left_start,right_start,left_end,right_end,pricky=None):
    if pricky is None:
        pricky = {}
    nn_lr=get_nn_list(array_right[right_start:right_end],array_left[left_start:left_end])
    nn_rl=get_nn_list(array_left[left_start:left_end],array_right[right_start:right_end])

    left_splits=[]
    right_splits=[]
    #print (len(nn_lr))
    #print (len(nn_rl))
    
    for item_lr in nn_lr:
        #print(item_lr.q)
        for item_rl in nn_rl:
            #print (item_rl.p)
            if all(item_rl.p==item_lr.q) and all(item_rl.q==item_lr.p):
                left_index = np.where(np.all(array_left == item_lr.q, axis=1))[0][0]
                right_index = np.where(np.all(array_right == item_lr.p, axis=1))[0][0]
                pricka=LineString([item_lr.p,item_lr.q])
                #print (pricka)
                #pricky=np.append(pricky, pricka)
                pricky[left_index]=pricka
                left_splits.append(left_index)
                right_splits.append(right_index)
    print(len(pricky))
    #asi checknout krizeni    
    left_splits.append(left_end)
    right_splits.append(right_end)
    
    lower_index_left=None
    lower_index_right=None
        
    for left_split,right_split in zip(left_splits,right_splits):
        if lower_index_left==None and lower_index_right==None:
            lower_index_left=left_start
            lower_index_right=right_end
        else:
            lower_index_left=upper_index_left+1
            lower_index_right=upper_index_right+1
            
        upper_index_left=left_split
        upper_index_right=right_split    
        if upper_index_left-lower_index_left >=2 and upper_index_right-lower_index_right >=2:
            get_all_pricky(array_left,array_right,lower_index_left,lower_index_right,upper_index_left,upper_index_right,pricky)
    return pricky


pricky_final=get_all_pricky(array_left,array_right,0,0,len(array_left),len(array_right))        
"""
print(len(pricky))
for pricka in pricky:
    kolize_right=shapely.intersection(pricka,linestring_right)
    if type(kolize_right)==shapely.geometry.multipoint.MultiPoint:
        pricky.remove(pricka)
    else:
        kolize_left=shapely.intersection(pricka,linestring_left)
        if type(kolize_left)==shapely.geometry.multipoint.MultiPoint:
            pricky.remove(pricka)
"""
sorted_pricky_final = dict(sorted(pricky_final.items()))

final_array=np.array(sorted_pricky_final.values())

print(final_array)

#print(pricky_final)    
schema = {
    'geometry': 'LineString',
    'properties': {'id': 'int'},
}
# Write a new Shapefile
with fiona.open('data/vysledek.shp', 'w', 'ESRI Shapefile', schema) as c:
    ## If there are multiple geometries, put the "for" loop here

    for key,value in sorted_pricky_final.items():
        c.write({
            'geometry': mapping(value),
            'properties': {'id': int(key)},
        })