import fiona
import time
import shapely
from math import sqrt,inf
import numpy as np
from shapely.geometry import shape, mapping, LineString, Point
from matplotlib import pyplot as plt
from sklearn.neighbors import NearestNeighbors
import pandas as pd
from scipy.optimize import linear_sum_assignment
schema_line = {
    'geometry': 'LineString',
}

schema_point = {
    'geometry': 'Point',
}

schema_multipolygon = {
    'geometry': 'MultiPolygon',
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
            
def get_cross_products(entire_from_bank,chosen_from_indices,entire_to_bank,chosen_to_indices):
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

def mask_for_left(cross_product1,cross_product2):
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
    return closest_query_indices,unique_pop_indices

def get_indices_onesided(queries,population):
    neigh.fit(population)
    distances,indices_to=neigh.kneighbors(queries)
    distances=distances.reshape(-1,)
    indices_to=indices_to.reshape(-1,)
    return distances,indices_to

def bank_to_bank (from_bank,to_bank):
    distances,indices_to=get_indices_onesided(from_bank,to_bank)
    closest_from_indices,unique_to_indices=get_unique_shortest(distances,indices_to)
    return closest_from_indices,unique_to_indices


def find_optimal_unique_closest_vertices(left_bank, right_bank):
    # Calculate the distance matrix
    distances = np.linalg.norm(left_bank[:, np.newaxis] - right_bank, axis=2)
    
    # Use the linear sum assignment (Hungarian algorithm) to find the optimal assignment
    row_ind, col_ind = linear_sum_assignment(distances)

    # Return the indices of the closest vertices on the right bank and their corresponding distances
    return col_ind, row_ind

def get_gaps(left,right):
    if len(right) != len(left):
        print("hmmmm")
    
    left_gaps=np.linalg.norm((left[1:]-left[:-1]),axis=1)
    right_gaps=np.linalg.norm((right[1:]-right[:-1]),axis=1)
    all_gaps=np.concatenate((left_gaps,right_gaps))
    return all_gaps

    

linestring_left=convert_shp_to_linestrings("left")
linestring_right=convert_shp_to_linestrings("right")
linestring_islands=convert_shp_to_linestrings("islands")

array_left=segmentize_array(linestring_left[0],20)
array_right=segmentize_array(linestring_right[0],20)

island_polygons=shapely.polygonize(linestring_islands)

multipolygon_islands=shapely.multipolygons(shapely.get_parts(island_polygons))

with fiona.open(f"results/islands_polygons.shp", 'w', 'ESRI Shapefile', schema_multipolygon) as c:
    c.write({
        'geometry': mapping(multipolygon_islands),
    })



left_lr,right_lr=bank_to_bank(array_left,array_right)
right_rl,left_rl=bank_to_bank(array_right,array_left)

mask_lr=mask_for_left(*get_cross_products(array_left,left_lr,array_right,right_lr))
mask_rl=mask_for_right(*get_cross_products(array_right,right_rl,array_left,left_rl))

left_lr_masked=left_lr[mask_lr]
right_lr_masked=right_lr[mask_lr]

left_rl_masked=left_rl[mask_rl]
right_rl_masked=right_rl[mask_rl]

all_left=np.concatenate((left_lr_masked,left_rl_masked),axis=None)
all_right=np.concatenate((right_lr_masked,right_rl_masked),axis=None)

all_left_sorted=np.sort(np.unique(all_left))
all_right_sorted=np.sort(np.unique(all_right))

print(len(all_left_sorted))
print(len(all_right_sorted))


lr_ids=np.stack((left_lr[mask_lr],right_lr[mask_lr]),axis=1)
rl_ids=np.stack((left_rl[mask_rl],right_rl[mask_rl]),axis=1)

duplicates=get_duplicates(lr_ids,rl_ids)
left_duplicate_ids,right_duplicate_ids=np.swapaxes(duplicates,0,1)



start_left=left_duplicate_ids[0]
start_right=right_duplicate_ids[0]

left_indices=[np.array([start_left])]
right_indices=[np.array([start_right])]


for left_id,right_id in zip(left_duplicate_ids[1:],right_duplicate_ids[1:]):
    end_left=left_id
    end_right=right_id
    
    left_filtered=all_left_sorted[np.logical_and(all_left_sorted>start_left,all_left_sorted<end_left)]
    right_filtered=all_right_sorted[np.logical_and(all_right_sorted>start_right,all_right_sorted<end_right)]
    
    if len(left_filtered) > 0 or len(right_filtered) > 0:
        ids_right,ids_left=find_optimal_unique_closest_vertices(array_left[left_filtered],array_right[right_filtered])
        
        left_bpm_with_ends=np.array([start_left,*left_filtered[ids_left],end_left])
        right_bpm_with_ends=np.array([start_right,*right_filtered[ids_right],end_right])
        
        left_lr_filtered=left_lr_masked[np.logical_and(left_lr_masked>=start_left,left_lr_masked<=end_left)]
        right_lr_filtered=right_lr_masked[np.logical_and(right_lr_masked>=start_right,right_lr_masked<=end_right)]
        
        left_rl_filtered=left_rl_masked[np.logical_and(left_rl_masked>=start_left,left_rl_masked<=end_left)]
        right_rl_filtered=right_rl_masked[np.logical_and(right_rl_masked>=start_right,right_rl_masked<=end_right)]
        
        if len(left_lr_filtered) == len(right_lr_filtered) and len(left_rl_filtered) == len(right_rl_filtered):
        
            gaps_bpm=get_gaps(array_left[left_bpm_with_ends],array_right[right_bpm_with_ends])
            gaps_lr=get_gaps(array_left[left_lr_filtered],array_right[right_lr_filtered])
            gaps_rl=get_gaps(array_left[left_rl_filtered],array_right[right_rl_filtered])
            
            max_gap_bpm=np.max(gaps_bpm)
            max_gap_lr=np.max(gaps_lr)
            max_gap_rl=np.max(gaps_rl)
            
            
            smallest_gap=min(max_gap_bpm,max_gap_lr,max_gap_rl)
            
            if smallest_gap == max_gap_bpm:
                left_indices.append(left_filtered[ids_left])
                right_indices.append(right_filtered[ids_right])
            elif smallest_gap == max_gap_lr:
                left_indices.append(left_lr_filtered[1:-1])
                right_indices.append(right_lr_filtered[1:-1])
            else:
                left_indices.append(left_rl_filtered[1:-1])
                right_indices.append(right_rl_filtered[1:-1])  
        
        
        
    left_indices.append(np.array([end_left]))
    right_indices.append(np.array([end_right]))    
    
    start_left=left_id
    start_right=right_id
    
final_left=np.concatenate(left_indices)
final_right=np.concatenate(right_indices)

final_banks=np.stack((array_left[final_left],array_right[final_right]),axis=0)
final_transects=banks=final_banks.swapaxes(0,1)


skoro_kilometraz=np.linalg.norm(np.diff(final_banks, axis=1), axis=-1)

distances_between_transects=np.mean(skoro_kilometraz,axis=0)/1000

cum_km=np.cumsum(distances_between_transects)

cum_km0=np.insert(cum_km,0,0)


effective_widths=[]
for transect in final_transects:
    transect_linestring=LineString((transect))
    effective_transect=shapely.difference(transect_linestring,multipolygon_islands)
    effective_widths.append(effective_transect.length)


plt.title("River width") 
plt.xlabel("kilometrage (km)") 
plt.ylabel("river width (m)") 
plt.plot(cum_km0,effective_widths)
plt.grid() 
plt.show()


transect_dict["combination"]=np.stack((array_left[final_left],array_right[final_right]),axis=1)

transect_dict["left_to_right"]=np.stack((array_left[left_lr_masked],array_right[right_lr_masked]),axis=1)

transect_dict["right_to_left"]=np.stack((array_left[left_rl_masked],array_right[right_rl_masked]),axis=1)

transect_dict["duplos"]=np.stack((array_left[left_duplicate_ids],array_right[right_duplicate_ids]),axis=1)




for key,value in transect_dict.items():
    export_line(value,key)
