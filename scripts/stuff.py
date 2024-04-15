import fiona
import time
import shapely
from math import sqrt,inf
import numpy as np
from shapely.geometry import shape, mapping, LineString, Point
from matplotlib import pyplot as plt
from sklearn.neighbors import NearestNeighbors

start=time.time()

with fiona.open("data/meandry_leva.shp") as levy_breh:
    linestring_left=shape(levy_breh[0]['geometry'])
        
        
with fiona.open("data/meandry_prava.shp") as pravy_breh:
    linestring_right=shape(pravy_breh[0]['geometry'])
    
islands=[]

with fiona.open("data/orlice_25_ostrovy.shp") as ostrovy:
    for feature in ostrovy:
        linestring_islands=shape(feature['geometry'])
        islands.append(linestring_islands)

def segmenty_array(vstup):
    segmentized=shapely.segmentize(vstup,25)
    arrayed=np.array(segmentized.coords)
    return arrayed
    
       
seg_islands=list(map(segmenty_array,islands))
    
seg_linestring_left=shapely.segmentize(linestring_left,25)
seg_linestring_right=shapely.segmentize(linestring_right,25)
#seg_linestring_islands=shapely.segmentize(linestring_islands,25)

array_left=np.array(seg_linestring_left.coords)
print(len(array_left))


schema_point = {
    'geometry': 'Point',
}

with fiona.open('data/seg_levy.shp', 'w', 'ESRI Shapefile', schema_point) as c:
    ## If there are multiple geometries, put the "for" loop here

    for breh in array_left:
        breh_point=Point(breh)
        c.write({
            'geometry': mapping(breh_point),
        })

array_right=np.array(seg_linestring_right.coords)
print(len(array_right))

with fiona.open('data/seg_pravy.shp', 'w', 'ESRI Shapefile', schema_point) as c:
    ## If there are multiple geometries, put the "for" loop here

    for breh in array_right:
        breh_point=Point(breh)
        c.write({
            'geometry': mapping(breh_point),
        })
        
        
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

def get_all(array_left,array_right,left_start,right_start,left_end,right_end,pricky=None):
    
    if pricky==None:
        pricky=[]
        
    left_splits=[]
    right_splits=[]
    
    pricky_rl=[]
    pricky_lr=[]

    neigh.fit(array_right[right_start:right_end])
    distances_lr,indices_r=neigh.kneighbors(array_left[left_start:left_end])

    unique_r_indices,closest_l_indices=get_unique_shortest(distances_lr,indices_r)
    
    neigh.fit(array_left[left_start:left_end])
    distances_rl,indices_l=neigh.kneighbors(array_right[right_start:right_end])

    unique_l_indices,closest_r_indices=get_unique_shortest(distances_rl,indices_l)

    for left_index1,right_index1 in zip(closest_l_indices,unique_r_indices):
        for right_index2,left_index2 in zip(closest_r_indices,unique_l_indices):
            if right_index1==right_index2 and left_index2==left_index1:
                left_index=left_index1+left_start
                right_index=right_index1+right_start
                left_splits.append(left_index)
                right_splits.append(right_index)
                pricka=[array_right[right_index],array_left[left_index]]
                pricky.append(pricka)
    """            
    for right_index,left_index in zip(closest_r_indices,unique_l_indices):
        left_splits.append(left_index)
        right_splits.append(right_index)
        pricka=[array_right[right_index],array_left[left_index]]
        pricky_rl.append(pricka)
        
    for left_index,right_index in zip(closest_l_indices,unique_r_indices):
        left_splits.append(left_index)
        right_splits.append(right_index)
        pricka=[array_right[right_index],array_left[left_index]]
        pricky_lr.append(pricka)
    """    
            
    lower_id_left=left_splits[0]
                
    for left_split in left_splits[1:]:
        lr_lower=closest_l_indices.index(lower_id_left)
        lr_upper=closest_l_indices.index(left_split)
        #np where :) a smrdis
        rl_lower=unique_l_indices.tolist().index(lower_id_left)
        rl_upper=unique_l_indices.tolist().index(left_split)
        
        number_of_pricky_lr=lr_upper-lr_lower
        number_of_pricky_rl=rl_upper-rl_lower
        
        lower_id_left=left_split
        
        if number_of_pricky_lr!=number_of_pricky_rl:
            winner=max(number_of_pricky_rl,number_of_pricky_lr)
            if winner==number_of_pricky_lr:
                for left_index,right_index in zip(closest_l_indices[lr_lower:lr_upper],unique_r_indices[lr_lower:lr_upper]):
                    pricka=[array_right[right_index],array_left[left_index]]
                    pricky_lr.append(pricka)
            else:
                for right_index,left_index in zip(closest_r_indices[rl_lower:rl_upper],unique_l_indices[rl_lower:rl_upper]):
                    pricka=[array_right[right_index],array_left[left_index]]
                    pricky_rl.append(pricka)
    
    """       
    left_splits.append(left_end)
    right_splits.append(right_end)
    
    lower_index_left=None
    lower_index_right=None
        
    for left_split,right_split in zip(left_splits,right_splits):
        if lower_index_left==None and lower_index_right==None:
            lower_index_left=left_start
            lower_index_right=right_start
        else:
            lower_index_left=upper_index_left+1
            lower_index_right=upper_index_right+1
            
        upper_index_left=left_split
        upper_index_right=right_split    
        if upper_index_left-lower_index_left >=10 and upper_index_right-lower_index_right >=10:
            get_all(array_left,array_right,lower_index_left,lower_index_right,upper_index_left,upper_index_right,pricky)"""
    return pricky_rl,pricky_lr,pricky

pricky_final_rl,pricky_final_lr,pricky_final=get_all(array_left,array_right,0,0,len(array_left),len(array_right))
    
schema_line = {
    'geometry': 'LineString',
}
# Write a new Shapefile

with fiona.open('data/vysledek_lr.shp', 'w', 'ESRI Shapefile', schema_line) as c:
    ## If there are multiple geometries, put the "for" loop here

    for pricka in pricky_final_lr:
        pricka_linestring=LineString(pricka)
        c.write({
            'geometry': mapping(pricka_linestring),
        })


with fiona.open('data/vysledek_rl.shp', 'w', 'ESRI Shapefile', schema_line) as c:
    ## If there are multiple geometries, put the "for" loop here

    for pricka in pricky_final_rl:
        pricka_linestring=LineString(pricka)
        c.write({
            'geometry': mapping(pricka_linestring),
        })

with fiona.open('data/vysledek_pokusicek2.shp', 'w', 'ESRI Shapefile', schema_line) as c:
    ## If there are multiple geometries, put the "for" loop here

    for pricka in pricky_final:
        pricka_linestring=LineString(pricka)
        c.write({
            'geometry': mapping(pricka_linestring),
        })