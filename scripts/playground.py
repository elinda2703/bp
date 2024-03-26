import fiona
import time
import shapely
from math import sqrt,inf
import numpy as np
from shapely.geometry import shape, mapping, LineString
from matplotlib import pyplot as plt 

start=time.time()

with fiona.open("data/meandry_leva.shp") as levy_breh:
    linestring_left=shape(levy_breh[0]['geometry'])
        
        
with fiona.open("data/meandry_prava.shp") as pravy_breh:
    linestring_right=shape(pravy_breh[0]['geometry'])
    
    
seg_linestring_left=shapely.segmentize(linestring_left,25)
seg_linestring_right=shapely.segmentize(linestring_right,25)


array_left=np.array(seg_linestring_left.coords)
print(len(array_left))
array_right=np.array(seg_linestring_right.coords)
array_right=np.flipud(array_right)
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


def get_all_pricky(array_left,array_right,left_start,right_start,left_end,right_end,pricky=None):
    if pricky is None:
        pricky = {}
    #nn_lr=get_nn_list(array_right[right_start:right_end],array_left[left_start:left_end])
    #nn_rl=get_nn_list(array_left[left_start:left_end],array_right[right_start:right_end])

    left_splits=[]
    right_splits=[]
    #print (len(nn_lr))
    #print (len(nn_rl))
    
    for vertex0 in array_left[left_start:left_end]:
        nnc=NearestNeighborsCouple.get_nn_couple(array_right[right_start:right_end],vertex0)
        for vertex1 in array_left[left_start:left_end]:
            if all(vertex1!=vertex0):
                dist=np.linalg.norm(vertex1 - nnc.p)
                if dist<nnc.distance:
                    break
                else:
                    continue
        if dist>nnc.distance:                
            left_index = np.where(np.all(array_left == vertex0, axis=1))[0][0]
            right_index = np.where(np.all(array_right == nnc.p, axis=1))[0][0]
            pricka=[vertex0,nnc.p]                    
            pricky[left_index]=pricka
            left_splits.append(left_index)
            right_splits.append(right_index)
            
        
    #print(len(pricky))
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

sorted_pricky_final = dict(sorted(pricky_final.items()))

sorted_pricky_list=list(sorted_pricky_final.values())

final_array=np.array(sorted_pricky_list)

widths = np.linalg.norm(np.diff(final_array, axis=1), axis=-1)

banks=final_array.swapaxes(0,1)

skoro_kilometraz=np.linalg.norm(np.diff(banks, axis=1), axis=-1)

distances_between_transects=np.mean(skoro_kilometraz,axis=0)/1000

cum_km=np.cumsum(distances_between_transects)

cum_km0=np.insert(cum_km,0,0)

end=time.time()

plt.title("Matplotlib demo") 
plt.xlabel("kilometrage (km)") 
plt.ylabel("river width (m)") 
plt.plot(cum_km0,widths)
plt.grid() 
plt.show()


total_time = end - start
print("\n"+ str(total_time))
#print(pricky_final)    
schema = {
    'geometry': 'LineString',
}
# Write a new Shapefile
with fiona.open('data/vysledek.shp', 'w', 'ESRI Shapefile', schema) as c:
    ## If there are multiple geometries, put the "for" loop here

    for pricka in final_array:
        pricka_linestring=LineString(pricka)
        c.write({
            'geometry': mapping(pricka_linestring),
        })

with fiona.open('data/novy_brehy.shp', 'w', 'ESRI Shapefile', schema) as c:
    ## If there are multiple geometries, put the "for" loop here

    for bank in banks:
        bank_linestring=LineString(bank)
        c.write({
            'geometry': mapping(bank_linestring),
        })