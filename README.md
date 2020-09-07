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
