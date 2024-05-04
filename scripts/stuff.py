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
    segmentized=shapely.segmentize(vstup,5)
    arrayed=np.array(segmentized.coords)
    return arrayed
    
       
seg_islands=list(map(segmenty_array,islands))
    
seg_linestring_left=shapely.segmentize(linestring_left,5)
seg_linestring_right=shapely.segmentize(linestring_right,5)
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
right_bank_plus_islands=np.concatenate((array_right,*seg_islands),axis=0)
left_bank_plus_islands=np.concatenate((array_left,*seg_islands),axis=0)        

with fiona.open('data/seg_pravy_plus_islands.shp', 'w', 'ESRI Shapefile', schema_point) as c:
    ## If there are multiple geometries, put the "for" loop here

    for breh in right_bank_plus_islands:
        breh_point=Point(breh)
        c.write({
            'geometry': mapping(breh_point),
        })
neigh = NearestNeighbors(n_neighbors=1)


def get_unique_shortest(distances,pop_indices):
    unique_pop_indices = np.unique(pop_indices)
    closest_query_indices=[]
    for pop_index in unique_pop_indices:
        indices_pop_indices = np.where(pop_indices == pop_index)[0]
        closest_query_index = indices_pop_indices[np.argmin(distances[indices_pop_indices])]
        closest_query_indices.append(closest_query_index)
    closest_query_indices=np.array(closest_query_indices)
    return unique_pop_indices,closest_query_indices

def get_duplicates(array1,array2):
    all_values=np.concatenate((array1,array2),axis=0)
    unique_values,unique_count=np.unique(all_values,return_counts=True,axis=0)
    duplicates=unique_values[unique_count==2]
    return duplicates
    
def get_indices_onesided(queries,*population):
    population_merged=np.vstack(population)
    neigh.fit(population_merged) #population = right bank + all islands
    distances,indices_to=neigh.kneighbors(queries) #query = left bank

    unique_to_indices,closest_from_indices=get_unique_shortest(distances,indices_to) #unique_r_indices includes islands
    
    return closest_from_indices,unique_to_indices

def get_filter(to_bank,unique_to_indices):    
    without_islands_filter=unique_to_indices<len(to_bank)
    islands_only_filter=unique_to_indices>=len(to_bank)
    return without_islands_filter,islands_only_filter

def filter_by_islands(from_indices,to_indices,island_filter):

    final_from_indices=from_indices[island_filter] #to be compared
    final_to_indices=to_indices[island_filter] #to be compared
        
    return final_from_indices,final_to_indices


    
    

def get_all(array_left,array_right,seg_islands):
    
    
    
    unpacked_islands=np.concatenate((seg_islands),axis=0)

    
    from_left,to_right_plus_islands_unique=get_indices_onesided(array_left,array_right,unpacked_islands)
         
    without_islands_filter,islands_only_filter=get_filter(array_right,to_right_plus_islands_unique)
    left_indices_without_islands,right_indices_without_islands=filter_by_islands(from_left,to_right_plus_islands_unique,without_islands_filter) #to be compared with switched sides
    left_indices_to_islands_only,island_candidates_for_left=filter_by_islands(from_left,to_right_plus_islands_unique,islands_only_filter) #to be compared
    island_candidates_for_left_adjusted=island_candidates_for_left-len(array_right)
    island_candidates_for_left_coords=unpacked_islands[island_candidates_for_left_adjusted]
    island_indices_to_left,left_indices_for_islands=get_indices_onesided(island_candidates_for_left_coords,array_left) #to be compared
    
    
    """
    unique_r_indices_without_islands=unique_r_indices[unique_r_indices<len(array_right)] #to be compared
    closest_l_indices_without_islands=closest_l_indices[unique_r_indices<len(array_right)] #to be compared

    island_candidates_for_left=unique_r_indices[unique_r_indices>=len(array_right)] #to be compared
    closest_l_indices_only_islands=closest_l_indices[unique_r_indices>=len(array_right)] #to be compared
    
    
    
    
    #left_to_islands=np.stack((array_left[closest_l_indices_only_islands],unpacked_islands[island_candidates_for_left_adjusted]),axis=1)

    
    neigh.fit(array_left) #population = left bank only
    distances_islands_to_left,left_indices_for_islands=neigh.kneighbors(unpacked_islands[island_candidates_for_left_adjusted]) #query = only previously chosen island vertices
       
    
    unique_left_indices_to_islands_only,closest_island_indices_to_left=get_unique_shortest(distances_islands_to_left,left_indices_for_islands) #to be compared
    """
    #islands_to_left=np.stack((array_left[unique_left_indices_to_islands_only],unpacked_islands[island_candidates_for_left_adjusted][closest_island_indices_to_left]),axis=1)
    left_bank_plus_islands=np.concatenate((array_left,unpacked_islands),axis=0)
    neigh.fit(left_bank_plus_islands) #population = left bank + all islands
    distances_rl,indices_l=neigh.kneighbors(array_right) #query = right bank

    unique_l_indices,closest_r_indices=get_unique_shortest(distances_rl,indices_l) #unique_l_indices includes islands
    
    unique_l_indices_without_islands=unique_l_indices[unique_l_indices<len(array_left)] #to be compared
    closest_r_indices_without_islands=closest_r_indices[unique_l_indices<len(array_left)] #to be compared
    


    island_candidates_for_right=unique_l_indices[unique_l_indices>=len(array_left)]
    closest_r_indices_only_islands=closest_r_indices[unique_l_indices>=len(array_left)] #to be compared
    
    neigh.fit(array_right) #population = right bank only
    distances_islands_to_right,right_indices_for_islands=neigh.kneighbors(left_bank_plus_islands[island_candidates_for_right]) #query = only previously chosen island vertices   
    
    unique_right_indices_to_islands_only,closest_island_indices_to_right=get_unique_shortest(distances_islands_to_right,right_indices_for_islands) #to be compared

    """
    all_lr=np.stack((closest_l_indices_without_islands,unique_r_indices_without_islands),axis=1)
    all_rl=np.stack((unique_l_indices_without_islands,closest_r_indices_without_islands),axis=1)
    
    duplicate_transects=get_duplicates(all_lr,all_rl)
    left_indices,right_indices=np.swapaxes(duplicate_transects,0,1)
    transects=np.stack((array_left[left_indices],array_right[right_indices]),axis=1)"""
    
    all_left_to_islands=np.stack((left_indices_to_islands_only,island_candidates_for_left_adjusted),axis=1)
    all_islands_to_left=np.stack((left_indices_for_islands,island_candidates_for_left_adjusted[island_indices_to_left]),axis=1)
    
    duplicate_left_islands=get_duplicates(all_left_to_islands,all_islands_to_left)
    left_indices_is,island_indices_l=np.swapaxes(duplicate_left_islands,0,1)
    transects_left_islands=np.stack((array_left[left_indices_is],unpacked_islands[island_indices_l]),axis=1)
    
    """
    for left_index1,right_index1 in zip(closest_l_indices,unique_r_indices):
        for right_index2,left_index2 in zip(closest_r_indices,unique_l_indices):
            if right_index1==right_index2 and left_index2==left_index1:
                left_index=left_index1+left_start
                right_index=right_index1+right_start
                left_splits.append(left_index)
                right_splits.append(right_index)
                pricka=[array_right[right_index],array_left[left_index]]
                pricky.append(pricka)"""
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
    return transects_left_islands

trasects_left_islands_final=get_all(array_left,array_right,seg_islands)


    
schema_line = {
    'geometry': 'LineString',
}
# Write a new Shapefile
"""
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
"""
with fiona.open('data/vysledek_pokusicek3.shp', 'w', 'ESRI Shapefile', schema_line) as c:
    ## If there are multiple geometries, put the "for" loop here

    for pricka in left_islands_final:
        pricka_linestring=LineString(pricka)
        c.write({
            'geometry': mapping(pricka_linestring),
        })

with fiona.open('data/vysledek_pokusicek4.shp', 'w', 'ESRI Shapefile', schema_line) as c:
    ## If there are multiple geometries, put the "for" loop here

    for pricka in islands_to_left_final:
        pricka_linestring=LineString(pricka)
        c.write({
            'geometry': mapping(pricka_linestring),
        })

with fiona.open('data/vysledek_pokusicek5.shp', 'w', 'ESRI Shapefile', schema_line) as c:
    ## If there are multiple geometries, put the "for" loop here

    for pricka in trasects_left_islands_final:
        pricka_linestring=LineString(pricka)
        c.write({
            'geometry': mapping(pricka_linestring),
        })