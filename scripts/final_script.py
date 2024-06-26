import fiona
import time
import shapely
from math import sqrt,inf
import numpy as np
from shapely.geometry import shape, mapping, LineString, Point
from matplotlib import pyplot as plt
from sklearn.neighbors import NearestNeighbors
import pandas as pd

schema_line = {
    'geometry': 'LineString',
}

schema_point = {
    'geometry': 'Point',
}

transect_dict={}

neigh = NearestNeighbors(n_neighbors=1)

def convert_shp_to_linestrings(file):
    linestring_list=[]
    with fiona.open(f"data/{file}_clip.shp") as banks:
        for feature in banks:
            linestring_feature=shape(feature['geometry'])
            linestring_list.append(linestring_feature)
    return linestring_list

def segmentize_array(linestring,step):
    seg_linestring=shapely.segmentize(linestring,step)
    arr=np.array(seg_linestring.coords)
    return arr

def export_line(arr,file):
    with fiona.open(f"results/{file}.shp", 'w', 'ESRI Shapefile', schema_line) as c:
        for line in arr:
            linestring=LineString(line)
            c.write({
                'geometry': mapping(linestring),
            })

def export_points(arr,file):
    with fiona.open(f"results/{file}.shp", 'w', 'ESRI Shapefile', schema_point) as c:
        for point in arr:
            point_final=Point(point)
            c.write({
                'geometry': mapping(point_final),
            })
            
def eliminate_incorrect(entire_from_bank,chosen_from_indices,entire_to_bank,chosen_to_indices):
    from_bank_vectors=entire_from_bank[1:] - entire_from_bank[:-1]
    
    from_bank_vectors_plus_last=np.vstack((from_bank_vectors,from_bank_vectors[-1]))
    from_bank_vectors_plus_first=np.vstack((from_bank_vectors[0],from_bank_vectors))
    
    transect_vectors=entire_to_bank[chosen_to_indices]-entire_from_bank[chosen_from_indices]
    
    cross_product1=np.cross(from_bank_vectors_plus_last[chosen_from_indices],transect_vectors)
    cross_product2=np.cross(from_bank_vectors_plus_first[chosen_from_indices],transect_vectors)
    
    return cross_product1,cross_product2

def mask_for_right(cross_product1,cross_product2):
    mask=(cross_product1 > 0) & (cross_product2 > 0)
    return mask

def mask_for_left_and_islands(cross_product1,cross_product2):
    mask=(cross_product1 < 0) & (cross_product2 < 0)
    return mask

def get_duplicates(array1,array2):
    all_values=np.concatenate((array1,array2),axis=0)
    unique_values,unique_count=np.unique(all_values,return_counts=True,axis=0)
    duplicates=unique_values[unique_count==2]
    return duplicates

def get_unique_shortest(distances,pop_indices):
    unique_pop_indices = np.unique(pop_indices)
    closest_query_indices=[]
    for pop_index in unique_pop_indices:
        indices_pop_indices = np.where(pop_indices == pop_index)[0]
        min_id=np.argmin(distances[indices_pop_indices])
        closest_query_index = indices_pop_indices[min_id]
        closest_query_indices.append(closest_query_index)
    closest_query_indices=np.array(closest_query_indices)
    return unique_pop_indices,closest_query_indices

def get_indices_onesided(queries,population,unique):
    neigh.fit(population)
    distances,indices_to=neigh.kneighbors(queries)
    distances=distances.reshape(-1,)
    indices_to=indices_to.reshape(-1,)
    
    if unique==True:
        unique_to_indices,closest_from_indices=get_unique_shortest(distances,indices_to) #unique_r_indices includes islands
        return closest_from_indices,unique_to_indices
    else:
        return distances,indices_to

def get_island_intervals_on_banks(islands_list,bank):
    id_list_intervals=[]
    id_list_bank=[]
    islands_indices=[]
    distances_mean=[]
    for island in islands_list:
        distances,chosen_bank_indices=get_indices_onesided(island,bank,unique=False)
        mask_ingredients=eliminate_incorrect(island,np.arange(len(island)),bank,chosen_bank_indices)
        mask=mask_for_left_and_islands(*mask_ingredients)
        chosen_bank_indices_masked=chosen_bank_indices[mask]
        island_indices_masked=np.arange(len(chosen_bank_indices))[mask]
        distance_mean=np.mean(distances[mask])
        distances_mean.append(distance_mean)
        id_list_bank.append(chosen_bank_indices_masked)
        islands_indices.append(island_indices_masked)
        min_id=np.argmin(chosen_bank_indices_masked)
        minimum=chosen_bank_indices_masked[min_id]
        max_id=np.argmax(chosen_bank_indices_masked)
        maximum=chosen_bank_indices_masked[max_id]
        id_interval=np.arange(minimum,maximum+1)
        id_list_intervals.append(id_interval)
        
    distances_mean_sorted, id_list_intervals_sorted, islands_indices_sorted,id_list_bank_sorted,islands_list_sorted = zip(*sorted(zip(distances_mean, id_list_intervals, islands_indices,id_list_bank,islands_list)))
    
    used_bank_intervals=None
    used_island_intervals=None
    cumulative_islands_vertices_count=0
    counter=0
    transect_list=[]
    for bank_interval,island_indices,bank_ids in zip(id_list_intervals_sorted,islands_indices_sorted,id_list_bank_sorted):
        if counter==0:
            used_bank_intervals=bank_interval
            island_rank_interval=np.full_like(bank_interval,counter)
            island_rank_ids=np.full_like(bank_ids,counter)
            used_island_ids=[island_indices]
            used_bank_ids=[bank_ids]
        else:
            doubles,double_banks,double_used=np.intersect1d(bank_interval,np.vstack((used_bank_intervals)),return_indices=True)
            if doubles.size != 0:
                assigned_islands_to_intervals=island_rank_interval[double_used]
                assigned_islands_to_ids=assigned_islands_to_intervals[np.isin(bank_interval[double_banks],bank_ids)]
                
                confronted_islands=np.unique(assigned_islands_to_intervals)
                
                for island_id in confronted_islands:
                    from_islands,to_islands=get_indices_onesided(islands_list_sorted[counter][island_indices][np.isin(bank_ids,doubles[np.where(assigned_islands_to_intervals==island_id)])],islands_list_sorted[island_id],unique=True)
                    transects=np.stack((islands_list_sorted[counter][island_indices][np.isin(bank_ids,doubles[np.where(assigned_islands_to_intervals==island_id)])][from_islands],islands_list_sorted[island_id][to_islands]),axis=1)
                    transect_list.append(transects)
                np.put(island_rank_interval,double_used,counter)
                int_without_doubles=np.delete(bank_interval,double_banks,0)
                used_bank_intervals=np.concatenate((used_bank_intervals,int_without_doubles))
                island_rank_add=np.full_like(int_without_doubles,counter)
                island_rank_interval=np.concatenate((island_rank_interval,island_rank_add))
            else:   
                used_bank_intervals=np.concatenate((used_bank_intervals,bank_interval))
                island_rank_add=np.full_like(bank_interval,counter)
                island_rank_interval=np.concatenate((island_rank_interval,island_rank_add))
                island_rank_ids_add=np.full_like(bank_ids,counter)
                island_rank_ids=np.concatenate((island_rank_ids,island_rank_ids_add))                  
        counter+=1        
                
    return transect_list                 
            
"""          
def get_island_intervals_on_banks_modified(islands_list,bank):
    id_list_intervals=[]
    id_list_bank=[]
    islands_indices=[]
    chosen_distances=[]
    island_tags=[]
    counter=0
    for island in islands_list:
        distances,chosen_bank_indices=get_indices_onesided(island,bank,unique=False)
        mask_ingredients=eliminate_incorrect(island,np.arange(len(island)),bank,chosen_bank_indices)
        mask=mask_for_left_and_islands(*mask_ingredients)
        chosen_bank_indices_masked=chosen_bank_indices[mask]
        island_indices_masked=np.arange(len(chosen_bank_indices))[mask]
        distances_masked=distances[mask]
        island_tags.append(np.full_like(distances_masked,counter))
        chosen_distances.append(distances_masked)
        id_list_bank.append(chosen_bank_indices_masked)
        islands_indices.append(island_indices_masked)
        counter+=1
    counter=0    
    for bank_ids in id_list_bank:
        intersect,current_bank_double_ids,rest_bank_double_ids=np.intersect1d(bank_ids,np.concatenate(id_list_bank[counter+1:]),return_indices=True)
        if intersect.size != 0:
           confronted_islands=np.concatenate(island_tags[counter+1:])[rest_bank_double_ids] 
           confronted_islands_unique=np.unique(np.concatenate(confronted_islands))
           for island_tag in confronted_islands_unique:
                rest_distances_mean=np.mean(distances_masked[island_tag][np. isin(bank_ids,intersect[np.where(island_tags[counter+1:][rest_bank_double_ids]==island_tag)])])
                current_distances_mean=np.mean(distances_masked[counter])"""
                
               
def islands_to_bank_intersection_check(island_list,bank_vertices,linestring_islands):
    cumsum_to_adjust=[0]
    lengths_cumsum=0
    for line in island_list[0:-1]:
        lengths_cumsum+=len(line)
        cumsum_to_adjust.append(lengths_cumsum)
    all_island_vertices_length=len(np.concatenate(island_list))
    final_from_ids_adjusted=[]
    final_to_ids_adjusted=[]
    to_ids_without_intersection_adjusted=[]
    counter_big=0
    used_intervals=[]
    
    final_linestrings=[]
    for island_origin in island_list:
        linestring_list=[]
        from_ids_with_intersections=[]
        from_ids_without_intersection=[]
        transects_without_intersection=[]
        to_indices_adjusted=[]
        to_ids_without_intersection_adjusted=[]
        distances,bank_indices=get_indices_onesided(island_origin,bank_vertices,unique=False)
        used_interval=np.arange(np.amin(bank_indices),np.amax(bank_indices)+1)
        used_intervals.append(used_interval)
        island_indices=np.arange(len(island_origin))
        mask_preparation=eliminate_incorrect(island_origin,island_indices,bank_vertices,bank_indices)
        mask=mask_for_left_and_islands(*mask_preparation)
        to_indices_oriented=bank_indices[mask]
        island_indices_oriented=island_indices[mask]
        
        to_indices_adjusted=to_indices_oriented+all_island_vertices_length
        island_indices_adjusted=island_indices_oriented+cumsum_to_adjust[counter_big]
        
        
        
        for to_id,from_id in zip(to_indices_oriented,island_indices_oriented):
            temporary_linestring=LineString((bank_vertices[to_id],island_origin[from_id]))
            linestring_list.append(temporary_linestring)
        counter_small=0
        for island_to_be_checked_ls,island_to_be_checked_arr in zip(linestring_islands,island_list):
            if counter_small==counter_big:
                counter_small+=1
                continue
            if transects_without_intersection and from_ids_without_intersection:
                linestring_list=transects_without_intersection
                island_indices_oriented=from_ids_without_intersection
                to_indices_adjusted=to_ids_without_intersection_adjusted
                transects_without_intersection=[]
                from_ids_without_intersection=[]
                from_ids_with_intersections=[]
                to_ids_without_intersection_adjusted=[]
                

            for linestring,from_id,to_id in zip(linestring_list,island_indices_oriented,to_indices_adjusted):
                intersection=shapely.intersection(linestring,island_to_be_checked_ls)
                if intersection.is_empty==False:
                    from_ids_with_intersections.append(from_id)
                else:
                    transects_without_intersection.append(linestring)
                    from_ids_without_intersection.append(from_id)
                    
                    to_ids_without_intersection_adjusted.append(to_id)
                    
                    
                        
            if from_ids_with_intersections:
                distances2,islands_to=get_indices_onesided(island_origin[from_ids_with_intersections],island_to_be_checked_arr,unique=False)
                islands_from=np.array(from_ids_with_intersections)
                for from_id,to_id in zip(islands_from,islands_to):
                    temporary_linestring=LineString((island_to_be_checked_arr[to_id],island_origin[from_id]))
                    transects_without_intersection.append(temporary_linestring)
                    from_ids_without_intersection.append(from_id)
                    to_ids_without_intersection_adjusted.append(to_id+cumsum_to_adjust[counter_small])
            counter_small+=1
        from_ids_adjusted=np.array(from_ids_without_intersection)+cumsum_to_adjust[counter_big]
        final_from_ids_adjusted.append(from_ids_adjusted)
        final_to_ids_adjusted.append(np.array(to_ids_without_intersection_adjusted))
        final_linestrings.append(transects_without_intersection)            
        counter_big+=1
    final_intervals=np.unique(np.concatenate(used_intervals))                 
    return final_linestrings,final_intervals,final_from_ids_adjusted,final_to_ids_adjusted            
        
        
                   
                
        

def islands_to_bank(islands_vertices, bank_vertices, name):
    distances,bank_indices=get_indices_onesided(islands_vertices,bank_vertices,unique=False)

    mask_preparation=eliminate_incorrect(islands_vertices,np.arange(len(islands_vertices)),bank_vertices,bank_indices)

    mask=mask_for_left_and_islands(*mask_preparation)

    transect_dict[f"{name}"]=np.stack((bank_vertices[bank_indices][mask],islands_vertices[mask]),axis=1)
    
    return bank_indices
    
    
linestring_left=convert_shp_to_linestrings("left")
linestring_right=convert_shp_to_linestrings("right")
linestring_islands=convert_shp_to_linestrings("islands")

array_left=segmentize_array(linestring_left[0],20)
array_right=segmentize_array(linestring_right[0],20)

arrays_islands=[]

for island in linestring_islands:
    array_island=segmentize_array(island,20)
    arrays_islands.append(array_island)
unpacked_islands=np.vstack(arrays_islands)

def final_islands(arrays_islands,array_bank,linestring_islands,name):    
    
    linestringos,intervalos,from_ids,to_ids=islands_to_bank_intersection_check(arrays_islands,array_bank,linestring_islands)
    unpacked_linestringos=sum(linestringos,[])
    unpacked_from_ids=np.concatenate(from_ids)
    unpacked_to_ids=np.concatenate(to_ids)
    transect_list=[]
    for transect in unpacked_linestringos:
        linestringos_array=segmentize_array(transect,inf)
        transect_list.append(linestringos_array)
    to_vertices,from_vertices=np.array(transect_list).swapaxes(0,1)

    widths=np.linalg.norm(to_vertices-from_vertices, axis=1)


    idk,indices_to_use=get_unique_shortest(widths,to_vertices)

    final_array=np.stack((to_vertices[indices_to_use],from_vertices[indices_to_use]),axis=1)
    
    
    export_line(final_array,f"{name}")
    return intervalos,unpacked_from_ids[indices_to_use],unpacked_to_ids[indices_to_use]
    
intervals_right,from_ids_right,to_ids_right=final_islands(arrays_islands,array_right,linestring_islands,"islands_to_right")
intervals_left,from_ids_left,to_ids_left=final_islands(arrays_islands,array_left,linestring_islands,"islands_to_left")

array_left_without_used=np.delete(array_left,intervals_left,axis=0)
array_right_without_used=np.delete(array_right,intervals_right,axis=0)

left_lr,right_lr=get_indices_onesided(array_left_without_used,array_right_without_used,unique=True)
right_rl,left_rl=get_indices_onesided(array_right_without_used,array_left_without_used,unique=True)

transect_dict["all_lr"]=np.stack((array_left_without_used[left_lr],array_right_without_used[right_lr]),axis=1)
transect_dict["all_rl"]=np.stack((array_left_without_used[left_rl],array_right_without_used[right_rl]),axis=1)


islands_plus_right=np.vstack((unpacked_islands,array_right))

transect_dict["adjusted_right"]=np.stack((islands_plus_right[to_ids_right],islands_plus_right[from_ids_right]),axis=1)

cumsum_to_adjust=[0]
lengths_cumsum=0
for line in arrays_islands:
        lengths_cumsum+=len(line)
        cumsum_to_adjust.append(lengths_cumsum)

panda_data=pd.DataFrame({'from_ids': from_ids_right, 'to_ids': to_ids_right, "width":np.linalg.norm(islands_plus_right[to_ids_right]-islands_plus_right[from_ids_right],axis=1)})


cumsum_to_adjust.append(panda_data['to_ids'].max() + 1)

# Use pd.cut to bin the data into specified intervals
panda_data['interval'] = pd.cut(panda_data['to_ids'], bins=cumsum_to_adjust, right=False)

# Group by the intervals
grouped = panda_data.groupby('interval')

for key, item in grouped:
    print(grouped.get_group(key), "\n\n")

      
#chosen_bank_indices_left=islands_to_bank(unpacked_islands,array_left,"transects_il")
#chosen_bank_indices_right=islands_to_bank(unpacked_islands,array_right,"transects_ir")


    
    
    

#all_lall=get_indices_onesided(array_left[chosen_bank_indices_left],right_plus_islands,unique=False)

#transect_dict["transects_lall"]=np.stack((array_left[chosen_bank_indices_left],right_plus_islands[all_lall]),axis=1)      
    


for key,value in transect_dict.items():
    export_line(value,key)
