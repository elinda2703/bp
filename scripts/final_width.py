import fiona
import os
import shapely
import math
import numpy as np
from shapely.geometry import shape, mapping, LineString, Point
from sklearn.neighbors import NearestNeighbors

from scipy.optimize import linear_sum_assignment
schema_line = {
    'geometry': 'LineString',
    'properties': {'width': 'float','kilometrag':'float'},
    
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
    with fiona.open(f"important/miss/{file}.shp") as banks:
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
    #vektorove souciny odpovidajich segmentu brehu a pricek
    
    #vektory po sobe jdoucich vrcholu celeho brehu ve smeru toku
    from_bank_vectors=entire_from_bank[1:] - entire_from_bank[:-1]
     
    #pridani posledniho vektoru posledniho segmentu brehu na konec
    from_bank_vectors_plus_last=np.vstack((from_bank_vectors,from_bank_vectors[-1]))
    
    #pridani posledniho vektoru prvniho segmentu brehu na zacatek 
    from_bank_vectors_plus_first=np.vstack((from_bank_vectors[0],from_bank_vectors)) 
    
    #smerove vektory pricek od vychoziho brehu
    transect_vectors=entire_to_bank[chosen_to_indices]-entire_from_bank[chosen_from_indices]
     
    #vektorovy soucin
    cross_product1=np.cross(from_bank_vectors_plus_last[chosen_from_indices],transect_vectors)
    
    #druhy vektorovy soucin tak, abychom meli vektorovy soucin pro oba segmenty, meyi kterymi pricka vznikla 
    cross_product2=np.cross(from_bank_vectors_plus_first[chosen_from_indices],transect_vectors) 
    
    return cross_product1,cross_product2

def mask_for_right(cross_product1,cross_product2):
    # maska pro pricky vychazejici z praveho brehu, jestli opravdu smeruji doleva
    mask=(cross_product1 > 0) & (cross_product2 > 0)
    return mask

def mask_for_left(cross_product1,cross_product2):
    # dtto pro levy breh
    mask=(cross_product1 < 0) & (cross_product2 < 0)
    return mask

def get_duplicates(array1,array2):
    # najde pricky, ktere se vyskutiju v obou stech, numpy nema 2d intersect, tak byla zvoleno np.unique
    all_values=np.concatenate((array1,array2),axis=0)
    unique_values,unique_count=np.unique(all_values,return_counts=True,axis=0)
    duplicates=unique_values[unique_count==2]
    return duplicates

def get_unique_shortest(distances,pop_indices):
    # odstraneni shluku
    unique_pop_indices = np.unique(pop_indices)
    closest_query_indices=[]
    # iterace pres unikatni indexy vybrane na cilovem brehu a vybrani nejkratsi pricky, resp. jejiho indexu
    for pop_index in unique_pop_indices:
        indices_pop_indices = np.where(pop_indices == pop_index)[0]
        min_id=np.argmin(distances[indices_pop_indices])
        closest_query_index = indices_pop_indices[min_id]
        closest_query_indices.append(closest_query_index)
    closest_query_indices=np.array(closest_query_indices)
    return closest_query_indices,unique_pop_indices

def get_indices_onesided(queries,population):
    # hledani nejblizsiho souseda (knihovna sci-kit), nasledny reshape
    neigh.fit(population)
    distances,indices_to=neigh.kneighbors(queries)
    distances=distances.reshape(-1,)
    indices_to=indices_to.reshape(-1,)
    return distances,indices_to

def bank_to_bank (from_bank,to_bank):
    # nejblizsi soused a odstrneni shluku, vrati indexy vrcholu na vychozim a cilovem brehu, dulezite poradi
    distances,indices_to=get_indices_onesided(from_bank,to_bank)
    closest_from_indices,unique_to_indices=get_unique_shortest(distances,indices_to)
    return closest_from_indices,unique_to_indices


def find_optimal_unique_closest_vertices(left_bank, right_bank):

    distances = np.linalg.norm(left_bank[:, np.newaxis] - right_bank, axis=2)
    row_ind, col_ind = linear_sum_assignment(distances)
    return col_ind, row_ind


def get_gaps(left,right):
    # spocita delku mezer pro vsechny dvojice vrcholu pricek na jednom brehu, pouzito pro useky mezi "dvojitymi prickami"
    left_gaps=np.linalg.norm((left[1:]-left[:-1]),axis=1)
    right_gaps=np.linalg.norm((right[1:]-right[:-1]),axis=1)
    all_gaps=np.concatenate((left_gaps,right_gaps))
    return all_gaps

path_island = 'important/miss/islands.shp'
isFile = os.path.isfile(path_island)    

linestring_left=convert_shp_to_linestrings("left")
linestring_right=convert_shp_to_linestrings("right")
#linestring_center=convert_shp_to_linestrings("line")
if isFile==True:
    linestring_islands=convert_shp_to_linestrings("islands")
    island_polygons=shapely.polygonize(linestring_islands)
    multipolygon_islands=shapely.multipolygons(shapely.get_parts(island_polygons))
else:
    multipolygon_islands=None

array_left=segmentize_array(linestring_left[0],200)
array_right=segmentize_array(linestring_right[0],200)

left_lr,right_lr=bank_to_bank(array_left,array_right)
mask_lr=mask_for_left(*get_cross_products(array_left,left_lr,array_right,right_lr))
left_lr_masked=left_lr[mask_lr]
right_lr_masked=right_lr[mask_lr]

right_rl,left_rl=bank_to_bank(array_right,array_left)
mask_rl=mask_for_right(*get_cross_products(array_right,right_rl,array_left,left_rl))
left_rl_masked=left_rl[mask_rl]
right_rl_masked=right_rl[mask_rl]

lr_ids=np.stack((left_lr[mask_lr],right_lr[mask_lr]),axis=1)
rl_ids=np.stack((left_rl[mask_rl],right_rl[mask_rl]),axis=1)

duplicates=get_duplicates(lr_ids,rl_ids)
left_duplicate_ids,right_duplicate_ids=np.swapaxes(duplicates,0,1)

all_left=np.concatenate((left_lr_masked,left_rl_masked),axis=None)
all_right=np.concatenate((right_lr_masked,right_rl_masked),axis=None)

all_left_sorted=np.sort(np.unique(all_left))
all_right_sorted=np.sort(np.unique(all_right))

print(len(all_left_sorted))
print(len(all_right_sorted))

left_duplicate_ids_plus_end=np.concatenate((left_duplicate_ids,[all_left_sorted[-1]]),axis=0)
right_duplicate_ids_plus_end=np.concatenate((right_duplicate_ids,[all_right_sorted[-1]]),axis=0)

start_left=0
start_right=0

left_indices=[]
right_indices=[]


for left_id,right_id in zip(left_duplicate_ids_plus_end,right_duplicate_ids_plus_end):
    end_left=left_id
    end_right=right_id
    
    left_filtered=all_left_sorted[np.logical_and(all_left_sorted>start_left,all_left_sorted<end_left)]
    right_filtered=all_right_sorted[np.logical_and(all_right_sorted>start_right,all_right_sorted<end_right)]
    
    if len(left_filtered) > 0 and len(right_filtered) > 0:
        ids_right,ids_left=find_optimal_unique_closest_vertices(array_left[left_filtered],array_right[right_filtered])
        
        left_bpm_with_ends=np.array([start_left,*left_filtered[ids_left],end_left])
        right_bpm_with_ends=np.array([start_right,*right_filtered[ids_right],end_right])
        
        filter_lr = np.where((right_lr_masked>=start_right) & (right_lr_masked<=end_right))[0]
        
        left_lr_filtered=left_lr_masked[filter_lr]
        right_lr_filtered=right_lr_masked[filter_lr]
        
        filter_check_lr=(left_lr_filtered>=start_left) & (left_lr_filtered<=end_left)
        
        left_lr_final=left_lr_filtered[filter_check_lr]
        right_lr_final=right_lr_filtered[filter_check_lr]
        
        filter_rl = np.where((left_rl_masked>=start_left) & (left_rl_masked<=end_left))[0]
        
        left_rl_filtered=left_rl_masked[filter_rl]
        right_rl_filtered=right_rl_masked[filter_rl]
        
        filter_check_rl=(right_rl_filtered>=start_right) & (right_rl_filtered<=end_right)
        
        left_rl_final=left_rl_filtered[filter_check_rl]
        right_rl_final=right_rl_filtered[filter_check_rl]
        
        
        if len(left_lr_final) == len(right_lr_final) and len(left_rl_final) == len(right_rl_final):
            
            gaps_lr=np.flip(np.sort(get_gaps(array_left[left_lr_final],array_right[right_lr_final])))
            gaps_rl=np.flip(np.sort(get_gaps(array_left[left_rl_final],array_right[right_rl_final])))
            gaps_bpm=np.flip(np.sort(get_gaps(array_left[left_bpm_with_ends],array_right[right_bpm_with_ends])))
            
            if np.array_equal(gaps_lr,gaps_bpm) or np.array_equal(gaps_rl,gaps_bpm):
                if np.array_equal(gaps_lr,gaps_rl):
                    left_indices.append(left_filtered[ids_left])
                    right_indices.append(right_filtered[ids_right])

                else:    
                    for gap_lr,gap_rl in zip(gaps_lr,gaps_rl):
                        if gap_rl == gap_lr:
                            continue
                        elif gap_rl < gap_lr:
                                left_indices.append(left_rl_filtered[1:-1])
                                right_indices.append(right_rl_filtered[1:-1])
                                break
                        else:
                                left_indices.append(left_lr_filtered[1:-1])
                                right_indices.append(right_lr_filtered[1:-1])
                                break
                        
            else:
                filt=np.array([0,1,2])
                for gap_lr,gap_rl,gap_bpm in zip(gaps_lr,gaps_rl,gaps_bpm):
                    arr=np.array([gap_lr,gap_rl,gap_bpm])
                    min_gap_length=np.min(arr[filt])
                    min_gap_where=np.where(arr == min_gap_length)
                    filt_where=np.intersect1d(min_gap_where,filt)
                    if len(filt_where)==3:
                        continue
                    elif len(filt_where)==2:
                        filt=filt_where
                        continue
                    else:
                        if filt_where[0]==0:
                            left_indices.append(left_lr_filtered[1:-1])
                            right_indices.append(right_lr_filtered[1:-1])
                        elif filt_where[0]==1:
                            left_indices.append(left_rl_filtered[1:-1])
                            right_indices.append(right_rl_filtered[1:-1])
                        else:
                            left_indices.append(left_filtered[ids_left])
                            right_indices.append(right_filtered[ids_right])
                        break
                            
        
        
    if end_left!=start_left and end_right!=start_right:    
        left_indices.append(np.array([end_left]))
        right_indices.append(np.array([end_right]))    
    
    start_left=left_id
    start_right=right_id
    
final_left=np.concatenate(left_indices)
final_right=np.concatenate(right_indices)

final_banks=np.stack((array_left[final_left],array_right[final_right]),axis=0)
final_transects=final_banks.swapaxes(0,1)

effective_widths=[]
effective_transects=[]
kilometrages=[]

if multipolygon_islands !=None:
    for transect in final_transects:
        transect_linestring=LineString((transect))
        effective_transect=shapely.difference(transect_linestring,multipolygon_islands)
        #center_intersect=shapely.intersection(effective_transect,linestring_center[0])
        #kilometrage=linestring_center[0].project(center_intersect)
        #if math.isnan(kilometrage):
            #continue
        #kilometrages.append(kilometrage)
        effective_transects.append(effective_transect)
        effective_widths.append(effective_transect.length)
else:
    for transect in final_transects:
        effective_transect=LineString((transect))
        #center_intersect=shapely.intersection(effective_transect,linestring_center[0])
        #kilometrage=linestring_center[0].project(center_intersect)
        #if math.isnan(kilometrage):
            #continue
        #kilometrages.append(kilometrage)
        effective_transects.append(effective_transect)
        effective_widths.append(effective_transect.length)
        
kilometrages_arr=np.array(kilometrages)
effective_widths_arr=np.array(effective_widths)
"""
with fiona.open("results/final_transects.shp", 'w', 'ESRI Shapefile', schema_line) as c:
    for line,km in zip(effective_transects,kilometrages):
        c.write({
            'geometry': mapping(line),
            'properties': {'width': line.length,'kilometrag':km/1000}
        })"""


with fiona.open("important/miss/final_transects.shp", 'w', 'ESRI Shapefile', schema_line) as c:
    for line in effective_transects:
        c.write({
            'geometry': mapping(line),
            'properties': {'width': line.length}
        })

np.save(f'kilometrages.npy', kilometrages_arr)
np.save(f'widths.npy', effective_widths_arr)

for key,value in transect_dict.items():
    export_line(value,key)
