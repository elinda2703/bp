import fiona
import time
import shapely
from math import sqrt,inf
import numpy as np
from shapely.geometry import shape, mapping, LineString
from matplotlib import pyplot as plt
from sklearn.neighbors import NearestNeighbors

start=time.time()

with fiona.open("data/meandry_leva.shp") as levy_breh:
    linestring_left=shape(levy_breh[0]['geometry'])
        
        
with fiona.open("data/meandry_prava.shp") as pravy_breh:
    linestring_right=shape(pravy_breh[0]['geometry'])
    
"""with fiona.open("data/vzorecek_ostrovy.shp") as ostrovy:
    linestring_islands=shape(ostrovy[2]['geometry'])"""
    
seg_linestring_left=shapely.segmentize(linestring_left,5)
seg_linestring_right=shapely.segmentize(linestring_right,5)
#seg_linestring_islands=shapely.segmentize(linestring_islands,25)

array_left=np.array(seg_linestring_left.coords)
print(len(array_left))
array_right=np.array(seg_linestring_right.coords)
#array_right=np.flipud(array_right)
print(len(array_right))
#array_islands=np.array(seg_linestring_islands.coords)
#https://stackoverflow.com/questions/12369484/searching-for-k-nearest-points
neigh = NearestNeighbors(n_neighbors=1)

def get_unique_shortest(distances,pop_indices):
    unique_pop_indices = np.unique(pop_indices)
    closest_query_indices=[]
    for pop_index in unique_pop_indices:
        indices_pop_indices = np.where(pop_indices == pop_index)[0]
        closest_query_index = indices_pop_indices[np.argmin(distances[indices_pop_indices])]
        closest_query_indices.append(closest_query_index)
    return unique_pop_indices,closest_query_indices


pricky_lr=[]
pricky_rl=[]

neigh.fit(array_right)
distances_lr,indices_r=neigh.kneighbors(array_left)

unique_r_indices,closest_l_indices=get_unique_shortest(distances_lr,indices_r)
   
neigh.fit(array_left)
distances_rl,indices_l=neigh.kneighbors(array_right)

unique_l_indices,closest_r_indices=get_unique_shortest(distances_lr,indices_r)

for left_index1,right_index1 in zip(closest_l_indices,unique_r_indices):
    for right_index2,left_index2 in zip(closest_r_indices,unique_l_indices):
        if right_index1==right_index2 and left_index2==left_index1:
            pricka=[array_right[right_index1],array_left[left_index1]]
            pricky_rl.append(pricka)



    
schema = {
    'geometry': 'LineString',
}
# Write a new Shapefile
"""
with fiona.open('data/vysledek_lr_meandry.shp', 'w', 'ESRI Shapefile', schema) as c:
    ## If there are multiple geometries, put the "for" loop here

    for pricka in pricky_lr:
        pricka_linestring=LineString(pricka)
        c.write({
            'geometry': mapping(pricka_linestring),
        })
"""

with fiona.open('data/vysledek_pokusicek.shp', 'w', 'ESRI Shapefile', schema) as c:
    ## If there are multiple geometries, put the "for" loop here

    for pricka in pricky_rl:
        pricka_linestring=LineString(pricka)
        c.write({
            'geometry': mapping(pricka_linestring),
        })
