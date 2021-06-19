# Mesher3DForSJTU

## 1 Usage
### 1.1 To generate 3d mesh from facets

```
>> mesher3d cube.poly -s 0.5 -o outmesh.mesh
```
-s sizing: to set the maximum scalar sizing of the mesh(equals to the maximum edge length)  
-o filename: to set the name of the output *.mesh file

### 1.2 To generate 3d mesh from points
```
>> mesher3d -s 5 sample.mesh -o outmesh.mesh
```
-s sizing: to set the maximum scalar sizing of the mesh(equals to the maximum edge length, regardless the points with marker 0)  
-o filename: to set the name of the output *.mesh file

### 1.3 To remesh a 3d mesh adaptively
You must keep 3 files in same path, *.mesh, *.remesh, *.value
```
>>mesher3d -r test3d
````
the chars after "-r" with no postfix
### 1.4 Quiet the mesh generation detials

```
>>mesher3d -r test3d -q
````
-q quiet mesh generation detials

### 1.5 To remesh a 3d mesh adaptively(with outer mesh regeneration) 
You must keep 3 files in same path, *.mesh, *.remesh, *.value
```
>>mesher3d -rr test3d -hmax 15 -hmin 3
````
the chars after "-rr" with no postfix

-hmax: max size

-hmin: min size 

### 1.5 To remesh a 3d mesh with a handle atomic area

```
>>mesher3d -hd test3d -s 15 -o out.3d
````
the chars after "-rr" with no postfix

-hmax: max size

-hmin: min size 

## 2 Remark
### 2.1 Labels for nodes

| Label | Significance |Tip|
|:------|:-------|:-----|
|0|Nodes inside the atomic area||
|1|Nodes on the border of the continuous area||
|2|Nodes on the border of the atomic area||
|3|Nodes between the border of the continuous area and the border of the atomic area||
### 2.2 Labels for tets

| Label | Significance |Tip|
|:------|:-------|:-----|
|0|Tets of the atomic area||
|1|Tets of the continuous area||
<!-- |4|Nodes on top border of the atomic area|with -hd|
|5|Nodes on top border of the continuous area|with -hd|
|6|Nodes on bottom border of the atomic area|with -hd|
|7|Nodes on bottom border of the continuous area|with -hd| -->

