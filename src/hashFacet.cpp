#include "hashFacet.h"




void HashFacetTable::insert(TriangleFacet aFacet){
    columns[key(aFacet)].push_back(aFacet);
}

int HashFacetTable::key(TriangleFacet aFacet){
    return 3*aFacet.orderedNodeIndices[0] + 5 * aFacet.orderedNodeIndices[1] + 7 *aFacet.orderedNodeIndices[2];
}


bool HashFacetTable::searchAnother(TriangleFacet keyFacet, TriangleFacet &goalFacet){

    int theKey = key(keyFacet);
    for(auto &f : columns[theKey]){

        if (f.tet != keyFacet.tet
        && f.orderedNodeIndices[0] == keyFacet.orderedNodeIndices[0] 
        && f.orderedNodeIndices[1] == keyFacet.orderedNodeIndices[1]
        && f.orderedNodeIndices[2] == keyFacet.orderedNodeIndices[2]){
            goalFacet = f;
            return true;
        }
    }

    return false;

}

bool HashFacetTable::search(TriangleFacet keyFacet, TriangleFacet &goalFacet){
    int theKey = key(keyFacet);
    for(auto &f : columns[theKey]){

        if (f.orderedNodeIndices[0] == keyFacet.orderedNodeIndices[0] 
        && f.orderedNodeIndices[1] == keyFacet.orderedNodeIndices[1]
        && f.orderedNodeIndices[2] == keyFacet.orderedNodeIndices[2]){
            goalFacet = f;
            return true;
        }
    }
    return false;    
}


bool HashFacetTable::remove(TriangleFacet keyFacet){
    int theKey = key(keyFacet);
    for(int i=0; i<columns[theKey].size(); i++){
        TriangleFacet f = columns[theKey][i];
        if (f.tet == keyFacet.tet
            && f.orderedNodeIndices[0] == keyFacet.orderedNodeIndices[0] 
            && f.orderedNodeIndices[1] == keyFacet.orderedNodeIndices[1]
            && f.orderedNodeIndices[2] == keyFacet.orderedNodeIndices[2]){
            columns[theKey].erase(columns[theKey].begin()+i);
            return true;
        }
    }

    return false;    
}