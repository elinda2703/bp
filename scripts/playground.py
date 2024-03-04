import fiona
import shapely
from math import sqrt,inf
import numpy as np
from shapely.geometry import shape, mapping, LineString

def pythagoras(x1, y1, x2, y2):
    # Computes distance between two points using Pythagoras theorem.
    distance = sqrt((x1-x2)**2+(y1-y2)**2)
    return distance

with fiona.open("data/vzorecek_leva.shp") as levy_breh:
    for feature in levy_breh:
        geometry = shape(feature['geometry'])
        leva_array = np.array(geometry.coords)

with fiona.open("data/vzorecek_prava.shp") as pravy_breh:
    for feature in pravy_breh:
        geometry = shape(feature['geometry'])
        prava_array=np.array(geometry.coords)

result_lp = []
result_pl=[]

n=0
for point1 in leva_array:
    shortest_dist=inf
    m=0
    for point2 in prava_array:
        distance = pythagoras(point1[0],point1[1],point2[0],point2[1])
        if distance<shortest_dist:
            shortest_dist=distance
            radek={"vertex":{"index":n,"souradnice":point1,"closest":{"1nn":{"index":m,"souradnice":point2,"vzdalenost":distance}}}}
        m+=1
    result_lp.append(radek)
    n+=1

n=0
for point1 in prava_array:
    shortest_dist=inf
    m=0
    for point2 in leva_array:
        distance = pythagoras(point1[0],point1[1],point2[0],point2[1])
        if distance<shortest_dist:
            shortest_dist=distance
            radek={"vertex":{"index":n,"souradnice":point1,"closest":{"1nn":{"index":m,"souradnice":point2,"vzdalenost":distance}}}}
        m+=1
    result_pl.append(radek)
    n+=1
#print(result_pl[231])
pricky=[]
for row_lp in result_lp:
    index_q_lp=row_lp["vertex"]["index"]
    index_p_lp=row_lp["vertex"]["closest"]["1nn"]["index"]
    for row_pl in result_pl:
        index_q_pl=row_pl["vertex"]["index"]
        index_p_pl=row_pl["vertex"]["closest"]["1nn"]["index"]
        if index_q_lp==index_p_pl and index_p_lp==index_q_pl:
            coords_leva=row_lp["vertex"]["souradnice"]
            coords_prava=row_pl["vertex"]["souradnice"]
            pricka=LineString([coords_leva,coords_prava])
            pricky.append(pricka)
print (len(pricky))

schema = {
    'geometry': 'LineString',
    'properties': {'id': 'int'},
}
# Write a new Shapefile
with fiona.open('data/vysledek.shp', 'w', 'ESRI Shapefile', schema) as c:
    ## If there are multiple geometries, put the "for" loop here
    x=0
    for item in pricky:
        c.write({
            'geometry': mapping(item),
            'properties': {'id': x},
        })
        x+=1

