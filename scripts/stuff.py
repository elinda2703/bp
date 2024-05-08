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

schema_line = {
    'geometry': 'LineString',
}

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
    unique_pop_indices_arr = np.unique(pop_indices)
    closest_query_indices=[]
    unique_pop_indices=[]
    for pop_index in unique_pop_indices_arr:
        indices_pop_indices = np.where(pop_indices == pop_index)[0]
        min_id=np.argmin(distances[indices_pop_indices])
        min_dist=distances[indices_pop_indices][min_id]
        if min_dist > 60:  #hodne prozatimni reseni
            continue
        closest_query_index = indices_pop_indices[min_id]
        closest_query_indices.append(closest_query_index)
        unique_pop_indices.append(pop_index)
    closest_query_indices=np.array(closest_query_indices)
    unique_pop_indices=np.array(unique_pop_indices)
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

def get_one_side(from_coords,to_coords,island_coords):
    from_unique_unfiltered,to_unique_unflitered=get_indices_onesided(from_coords,to_coords,island_coords)
         
    without_islands_filter,islands_only_filter=get_filter(to_coords,to_unique_unflitered)
    from_indices_without_islands,to_indices_without_islands=filter_by_islands(from_unique_unfiltered,to_unique_unflitered,without_islands_filter) #to be compared with switched sides
    from_indices_to_islands_only1,island_candidates=filter_by_islands(from_unique_unfiltered,to_unique_unflitered,islands_only_filter) #to be compared
    island_candidates_adjusted1=island_candidates-len(to_coords)
    island_candidates_coords=island_coords[island_candidates_adjusted1]
    island_indices,from_indices_for_islands2=get_indices_onesided(island_candidates_coords,from_coords)
    island_indices_final2=island_candidates_adjusted1[island_indices]
    
    return from_indices_without_islands,to_indices_without_islands,from_indices_to_islands_only1,island_candidates_adjusted1,from_indices_for_islands2,island_indices_final2

def get_transects(a1,a2,b1,b2):
    pass

def no_islands(from_coords,to_coords):
    pass

def export_line(arr,file):
    with fiona.open(f"data/{file}.shp", 'w', 'ESRI Shapefile', schema_line) as c:
        for line in arr:
            linestring=LineString(line)
            c.write({
                'geometry': mapping(linestring),
            })
   

def get_all(array_left,array_right,seg_islands):
    transect_dict={}
    unpacked_islands=np.concatenate((seg_islands),axis=0)
    banks_merged=np.concatenate((array_left,array_right),axis=0)
    

    
    islands_il,left_il=get_indices_onesided(unpacked_islands,array_left)
    islands_ir,right_ir=get_indices_onesided(unpacked_islands,array_right)
    
    #transect_dict["transects_li"]=np.stack((array_left[left_li],unpacked_islands[islands_li]),axis=1)
    transect_dict["transects_il"]=np.stack((array_left[left_il],unpacked_islands[islands_il]),axis=1)
    
    #transect_dict["transects_ri"]=np.stack((array_right[right_ri],unpacked_islands[islands_ri]),axis=1)
    transect_dict["transects_ir"]=np.stack((array_right[right_ir],unpacked_islands[islands_ir]),axis=1)
    
    lower_int=0
    upper_int=0
    
    left_isl_ints=[]
    right_isl_ints=[]
  
    for isl in seg_islands:
        upper_int+=len(isl)
        
        sel_isl_left=np.logical_and(islands_il>=lower_int,islands_il<upper_int)
        sel_isl_right=np.logical_and(islands_ir>=lower_int,islands_ir<upper_int)
        
        left_min_id=np.argmin(left_il[sel_isl_left])
        left_min=left_il[sel_isl_left][left_min_id]
        left_max_id=np.argmax(left_il[sel_isl_left])
        left_max=left_il[sel_isl_left][left_max_id]        
        right_min_id=np.argmin(right_ir[sel_isl_right])
        right_min=right_ir[sel_isl_right][right_min_id]
        right_max_id=np.argmin(right_ir[sel_isl_right])
        right_max=right_ir[sel_isl_right][right_max_id]
        
        left_isl_ints.append(np.arange(left_min,left_max))
        right_isl_ints.append(np.arange(right_min,right_max))
        
        lower_int+=len(isl)
        
    used_for_islands_left=np.concatenate((left_isl_ints),axis=0)
    used_for_islands_right=np.concatenate((right_isl_ints),axis=0)
    
    mask_left = np.ones_like(array_left, dtype=bool)
    mask_left[used_for_islands_left] = False
    array_left_without_islands = array_left[mask_left].reshape(-1,2)
    
    mask_right = np.ones_like(array_right, dtype=bool)
    mask_right[used_for_islands_right] = False
    array_right_without_islands = array_right[mask_right].reshape(-1,2)
    
    left_lr,right_lr=get_indices_onesided(array_left_without_islands,array_right_without_islands)
    right_rl,left_rl=get_indices_onesided(array_right_without_islands,array_left_without_islands)
    
    transect_dict["transects_lr"]=np.stack((array_left_without_islands[left_lr],array_right_without_islands[right_lr]),axis=1)
    transect_dict["transects_rl"]=np.stack((array_left_without_islands[left_rl],array_right_without_islands[right_rl]),axis=1)
      
    all_left_to_right=np.stack((left_lr,right_lr),axis=1)
    all_right_to_left=np.stack((left_rl,right_rl),axis=1)
    
    duplicate_bank_transects=get_duplicates(all_left_to_right,all_right_to_left)
    left_indices,right_indices=np.swapaxes(duplicate_bank_transects,0,1)
    transect_dict["transects_banks"]=np.stack((array_left_without_islands[left_indices],array_right_without_islands[right_indices]),axis=1)
    """
    #all_left_to_islands=np.stack((left_li,islands_li),axis=1)
    #all_islands_to_left=np.stack((left_il,islands_il),axis=1)

    duplicate_left_islands_transects=get_duplicates(all_left_to_islands,all_islands_to_left)
    left_indices_islands,islands_indices_left=np.swapaxes(duplicate_left_islands_transects,0,1)
    transect_dict["transects_left_islands"]=np.stack((array_left[left_indices_islands],unpacked_islands[islands_indices_left]),axis=1)
    """
    #all_right_to_islands=np.stack((right_ri,islands_ri),axis=1)
    #all_islands_to_right=np.stack((right_ir,islands_ir),axis=1)
    """
    duplicate_right_islands_transects=get_duplicates(all_right_to_islands,all_islands_to_right)
    right_indices_islands,islands_indices_right=np.swapaxes(duplicate_right_islands_transects,0,1)
    transect_dict["transects_right_islands"]=np.stack((array_right[right_indices_islands],unpacked_islands[islands_indices_right]),axis=1)
    
    
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
    return transect_dict

transects_final=get_all(array_left,array_right,seg_islands)

for key,value in transects_final.items():
    export_line(value,key)