#pragma once
#include <unordered_map>
#include <vector>
#include "mesh.h"


class HashFacetTable{
public:
    std::unordered_map<int, std::vector<TriangleFacet> > columns;

    HashFacetTable(){}
    void insert(TriangleFacet aFacet);
    int key(TriangleFacet aFacet);
    bool searchAnother(TriangleFacet keyFacet, TriangleFacet &goalFacet);
    bool search(TriangleFacet keyFacet, TriangleFacet &goalFacet);
    bool remove(TriangleFacet keyFacet);
};




