import fiona
import time
import shapely
from math import sqrt,inf
import numpy as np
from shapely.geometry import shape, mapping, LineString, Point
from matplotlib import pyplot as plt
from sklearn.neighbors import NearestNeighbors

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
    with fiona.open(f"data/{file}.shp") as banks:
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
        id_interval=np.arange(minimum,maximum)
        id_list_intervals.append(id_interval)
        
    distances_mean_sorted, id_list_intervals_sorted, islands_indices_sorted,id_list_bank_sorted,islands_list_sorted = zip(*sorted(zip(distances_mean, id_list_intervals, islands_indices,id_list_bank,islands_list)))
    
    used_bank_intervals=None
    used_island_intervals=None
    cumulative_islands_vertices_count=0
    counter=0
    
    for bank_interval,island_indices,bank_ids in zip(id_list_intervals_sorted,islands_indices_sorted,id_list_bank_sorted):
        if counter==0:
            used_bank_intervals=bank_interval
            island_rank=np.full_like(bank_interval,counter)
            used_island_ids=[island_indices]
            used_bank_ids=[bank_ids]
        else:
            doubles,double_banks,double_used=np.intersect1d(bank_interval,np.vstack((used_bank_intervals)),return_indices=True)
            if doubles.size != 0:
                confronted_islands=np.unique(island_rank[double_used])
                np.put(island_rank,double_banks,counter)
                int_without_doubles=np.delete(bank_interval,double_banks,0)
                used_bank_intervals=np.concatenate((used_bank_intervals,int_without_doubles))
                island_rank_add=np.full_like(int_without_doubles,counter)
                island_rank=np.concatenate((island_rank,island_rank_add))
                for island_id in confronted_islands:
                    from_islands,to_islands=get_indices_onesided(islands_list_sorted[counter][island_indices][np.isin(bank_ids,doubles)],islands_list_sorted[island_id],unique=True)
            else:   
                used_bank_intervals=np.concatenate((used_bank_intervals,bank_interval))
                island_rank_add=np.full_like(bank_interval,counter)
                island_rank=np.concatenate((island_rank,island_rank_add))                  
        counter+=1        
                
                
            
          
            
                
        

def islands_to_bank(islands_vertices, bank_vertices, name):
    islands_indices,bank_indices=get_indices_onesided(islands_vertices,bank_vertices,unique=True)

    mask_preparation=eliminate_incorrect(islands_vertices,islands_indices,bank_vertices,bank_indices)

    mask=mask_for_left_and_islands(*mask_preparation)

    transect_dict[f"{name}"]=np.stack((bank_vertices[bank_indices][mask],islands_vertices[islands_indices][mask]),axis=1)
    
    return bank_indices
    
    
linestring_left=convert_shp_to_linestrings("left")
linestring_right=convert_shp_to_linestrings("right")
linestring_islands=convert_shp_to_linestrings("islands")

array_left=segmentize_array(linestring_left[0],10)
array_right=segmentize_array(linestring_right[0],10)

arrays_islands=[]

for island in linestring_islands:
    array_island=segmentize_array(island,10)
    arrays_islands.append(array_island)
    
unpacked_islands=np.vstack(arrays_islands)

get_island_intervals_on_banks(arrays_islands,array_left)      
#chosen_bank_indices_left=islands_to_bank(unpacked_islands,array_left,"transects_il")
#chosen_bank_indices_right=islands_to_bank(unpacked_islands,array_right,"transects_ir")

right_plus_islands=np.vstack((array_right,unpacked_islands))

#all_lall=get_indices_onesided(array_left[chosen_bank_indices_left],right_plus_islands,unique=False)

#transect_dict["transects_lall"]=np.stack((array_left[chosen_bank_indices_left],right_plus_islands[all_lall]),axis=1)      
    
islands_to_bank(unpacked_islands,array_right,"transects_ir")
islands_to_bank(unpacked_islands,array_left,"transects_il")

for key,value in transect_dict.items():
    export_line(value,key)
