# Mesher3DForSJTU

## 1. To generate 3d mesh from facets

```
>> mesher3d cube.poly -s 0.5 -o outmesh.mesh
```
-s sizing: to set the maximum scalar sizing of the mesh(equals to the maximum edge length)  
-o filename: to set the name of the output *.mesh file

## 2. To generate 3d mesh from points
```
>> mesher3d -s 5 sample.mesh -o outmesh.mesh
```
-s sizing: to set the maximum scalar sizing of the mesh(equals to the maximum edge length, regardless the points with marker 0)  
-o filename: to set the name of the output *.mesh file

## 3. To remesh a 3d mesh adaptively
You must keep 3 files in same path, *.mesh, *.remesh, *.value
```
>>mesher3d -r test3d
````
the chars after "-r" with no postfix
## 4. Quiet the mesh generation detials

```
>>mesher3d -r test3d -q
````
-q quiet mesh generation detials

## 5. To remesh a 3d mesh adaptively(with outer mesh regeneration) 
You must keep 3 files in same path, *.mesh, *.remesh, *.value
```
>>mesher3d -rr test3d -hmax 15 -hmin 3
````
the chars after "-r" with no postfix

-hmax: max size

-hmin: min size 