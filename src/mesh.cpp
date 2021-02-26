#include "mesh.h"
#include "hashFacet.h"
#include "kdtree.h"
#include <unordered_map>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>

double SixTimesTetrahedronVolume(Vector3D v0, Vector3D v1, Vector3D v2, Vector3D v3){
    Vector3D a = v0 - v2;
    Vector3D b = v1 - v0;
    Vector3D c = v3 - v0;
    double rst = c[0]*(a[1]*b[2]-b[1]*a[2]) - c[1]*(a[0]*b[2]-b[0]*a[2]) + c[2]*(a[0]*b[1]-b[0]*a[1]);
    return rst;
}





void Mesh::rebuildIndices(){
    for(int i=0; i<nodes.size(); i++){
        nodes[i]->index = i;
    }

    for(int i=0; i<tetrahedrons.size(); i++){
        tetrahedrons[i]->index = i;
    }
}


void Mesh::rebuildTetrahedronsAdjacency(){
    HashFacetTable table;
    for(auto e: tetrahedrons){
        for(int i=0; i<4; i++){
            if(e->adjacentTetrahedrons[i] == nullptr){
                TriangleFacet keyTri = e->facet(i);
                TriangleFacet goalTri;
                if (table.searchAnother(keyTri, goalTri)){
                    Tetrahedron *goalTet = goalTri.tet;
                    e->adjacentTetrahedrons[i] = goalTet;
                    goalTet->adjacentTetrahedrons[goalTri.localIndex] = e;
                    table.remove(goalTri);
                }
                else{
                    table.insert(keyTri);
                }
            }
        }
    }
}

void Mesh::extractBorderNodeIndicesWithLabels(std::vector<int> labels, std::set<int> &nodeIndices){
    HashFacetTable table;
    auto satisfy=
    [&labels]
    (int label){
        bool rst=false;
        for(auto i: labels){
            if (label == i){
                rst = true;
                break;
            }
        }
        return rst;
    };

    for(auto &e: tetrahedrons){
        if ( satisfy(e->nodes[0]->label) 
        && satisfy(e->nodes[1]->label)
        && satisfy(e->nodes[2]->label)
        && satisfy(e->nodes[3]->label)){
            for(int i=0; i<4; i++){
                TriangleFacet keyTri = e->facet(i);
                TriangleFacet goalTri;
                if (table.search(keyTri, goalTri)){
                    table.remove(goalTri);    
                }
                else{
                    table.insert(keyTri);
                }    
            }
        }
    }
    for(auto kv: table.columns){
        for(auto f: kv.second){
            for(auto i: f.orderedNodeIndices){
                nodeIndices.insert(i);
            }
        }
    }
}

void Mesh::extractBorderNodes(std::vector<Node *> &sNodes){

    HashFacetTable table;
    for(auto &e: tetrahedrons){
        for(int i=0; i<4; i++){
            if (e->adjacentTetrahedrons[i] == nullptr){
                TriangleFacet keyTri = e->facet(i);
                TriangleFacet goalTri;
                if (table.search(keyTri, goalTri)){
                    table.remove(goalTri);    
                }
                else{
                    table.insert(keyTri);
                }
            }
        }
    }

    std::set<Node *> nodeSet;
    for(auto &kv: table.columns){
        for(auto &f: kv.second){
            for(auto i: f.orderedNodeIndices){
                nodeSet.insert(this->nodes[i]);
            }
        }
    }
    sNodes.insert(sNodes.end(), nodeSet.begin(), nodeSet.end());
}

void Mesh::extractBorder(std::vector<Node *> &sNodes, std::vector<TriangleFacet> &sFacets){
    HashFacetTable table;
    std::set<Node *> nodeSet;
    for(auto &e: tetrahedrons){
        for(int i=0; i<4; i++){
            if (e->adjacentTetrahedrons[i] == nullptr){
                TriangleFacet keyTri = e->facet(i, true);
                TriangleFacet goalTri;
                if (table.search(keyTri, goalTri)){
                    table.remove(goalTri);    
                }
                else{
                    table.insert(keyTri);
                }
            }
        }
    }

    for(auto &kv: table.columns){
        for(auto f: kv.second){
            sFacets.push_back(f);
            for(auto n: f.sNodes){
                nodeSet.insert(n);
            }
        }
    }
    sNodes.insert(sNodes.end(), nodeSet.begin(), nodeSet.end());
}

void Mesh::deleteLargeScaleTetrahedronsPermanently(std::vector<Tetrahedron *> &delElements){

    for(auto e: tetrahedrons){
        e->edit = 0;
    }

    for(auto n: nodes){
        n->edit = -1;
    }

    for(auto e: delElements){
        for(auto ee: e->adjacentTetrahedrons){
            if (ee){
                ee->removeAdjacentTetrahedron(e);
            }
        }
        e->edit = -1;
    }

    for(int i=0; i<tetrahedrons.size(); i++){
        if (tetrahedrons[i]->edit == -1){
            delete tetrahedrons[i];
            tetrahedrons.erase(tetrahedrons.begin()+i);
            i--;
        }
        else{
            for(auto n: tetrahedrons[i]->nodes){
                n->edit = 0;
            }
        }
    }

    for(int i=0; i<nodes.size(); i++){
        if (nodes[i]->edit == -1){
            delete nodes[i];
            nodes.erase(nodes.begin()+i);
            i--;
        }
    }

    rebuildIndices();

}


void Mesh::deleteSmallScaleTetrahedronsPermanently(std::vector<Tetrahedron *> &delElements){
    for(auto e: delElements){
        for(auto ee: e->adjacentTetrahedrons){
            if (ee){
                ee->removeAdjacentTetrahedron(e);
            }
        }
        tetrahedrons.erase(std::find(tetrahedrons.begin(), tetrahedrons.end(), e));
        
    }
    //TODO deleteSmallScaleTetrahedronsPermanently



}


void Mesh::mergeMesh(Mesh &anotherMesh, std::vector<Node *> &mergeNodes, bool deleteUselessNodes){
    std::vector<Node *> borderNodes;
    std::vector<Node *> anotherBorderNodes;
    for(auto n:anotherMesh.nodes){
        n->tempDate = nullptr;
    }

    extractBorderNodes(borderNodes);
    anotherMesh.extractBorderNodes(anotherBorderNodes);     
    struct kdtree *kd = kd_create(3);
    for(auto n: borderNodes){
        kd_insert(kd, n->pos.data(), static_cast<void*>(n));
    }

    for(auto n:anotherBorderNodes){
        struct kdres *search_rst =  kd_nearest_range(kd, n->pos.data(), 1e-13);
        if (kd_res_size(search_rst)){

            void *data = kd_res_item_data(search_rst);
            Node *goalNode = static_cast<Node*>(data);
            n->tempDate = data;
            mergeNodes.push_back(goalNode);
        }
        kd_res_free(search_rst);
    }

    for (auto e: anotherMesh.tetrahedrons){
        for(int i=0; i<4; i++){
            if (e->nodes[i]->tempDate){
                e->nodes[i] = static_cast<Node *>(e->nodes[i]->tempDate);    
            }
        }
        this->tetrahedrons.push_back(e);
    }

    for(auto n:anotherMesh.nodes){
        if (n->tempDate==nullptr){
            this->nodes.push_back(n);
        }
        else{
            if( deleteUselessNodes && n->tempDate!=n)
            delete n;
        }
    }
    kd_free(kd);
    
    this->rebuildIndices();
    this->rebuildTetrahedronsAdjacency();
    
}

void Mesh::splitElementFromCentre(Tetrahedron *tet){
    Node *insertNode = new Node();
    for(auto n:tet->nodes){
        insertNode->pos += n->pos;
    }
    insertNode->pos /= 4;
    Node *n0 = tet->nodes[0];
    Node *n1 = tet->nodes[1];
    Node *n2 = tet->nodes[2];
    Node *n3 = tet->nodes[3];
    Tetrahedron *t0 = new Tetrahedron(insertNode, n1, n2, n3);
    Tetrahedron *t1 = new Tetrahedron(n0, insertNode, n2, n3);
    Tetrahedron *t2 = new Tetrahedron(n0, n1, insertNode, n3);
    Tetrahedron *t3 = new Tetrahedron(n0, n1, n2, insertNode);

    t0->adjacentTetrahedrons[0] = tet->adjacentTetrahedrons[0];
    t0->adjacentTetrahedrons[1] = t1;
    t0->adjacentTetrahedrons[2] = t2;
    t0->adjacentTetrahedrons[3] = t3;
    if (t0->adjacentTetrahedrons[0]){
        t0->adjacentTetrahedrons[0]->replaceAdjacentTetrahedron(tet, t0);
    }

    t1->adjacentTetrahedrons[0] = t0;
    t1->adjacentTetrahedrons[1] = tet->adjacentTetrahedrons[1];
    t1->adjacentTetrahedrons[2] = t2;
    t1->adjacentTetrahedrons[3] = t3;
    if (t1->adjacentTetrahedrons[1]){
        t1->adjacentTetrahedrons[1]->replaceAdjacentTetrahedron(tet, t1);
    }

    t2->adjacentTetrahedrons[0] = t0;
    t2->adjacentTetrahedrons[1] = t1;
    t2->adjacentTetrahedrons[2] = tet->adjacentTetrahedrons[2];
    t2->adjacentTetrahedrons[3] = t3;
    if (t2->adjacentTetrahedrons[2]){
        t2->adjacentTetrahedrons[2]->replaceAdjacentTetrahedron(tet, t2);
    }

    t3->adjacentTetrahedrons[0] = t0;
    t3->adjacentTetrahedrons[1] = t1;
    t3->adjacentTetrahedrons[2] = t2;
    t3->adjacentTetrahedrons[3] = tet->adjacentTetrahedrons[3];
    if (t3->adjacentTetrahedrons[3]){
        t3->adjacentTetrahedrons[3]->replaceAdjacentTetrahedron(tet, t3);
    }

    tetrahedrons.erase(std::find(tetrahedrons.begin(), tetrahedrons.end(), tet));
    delete tet;
    
    this->nodes.push_back(insertNode);
    this->tetrahedrons.push_back(t0);
    this->tetrahedrons.push_back(t1);
    this->tetrahedrons.push_back(t2);
    this->tetrahedrons.push_back(t3);
}


void Mesh::readyForSpatialSearch(bool toBuildTetKDTree, bool toBuildNodeKDTree, bool toEstimateSizing){
    if (toBuildNodeKDTree){
        free(nodeKDTree);
        nodeKDTree = kd_create(3);
        for(auto n: nodes){
            kd_insert(nodeKDTree, n->pos.data(), static_cast<void*>(n));
        }
    }


    // estimateSizing();
    // tetRTree.RemoveAll();
    // for(auto e:tetrahedrons){
    //     e->generateBoundingBox();
    //     tetRTree.Insert(e->boundingBox.minimum.data(), e->boundingBox.minimum.data(), e);
    // }
    if(toEstimateSizing){
        estimateSizing();
    }


    if (toBuildTetKDTree){

        free(tetKDTree);
        tetKDTree = kd_create(3);
        for(auto e: tetrahedrons){
            e->generateBoundingBox();
            kd_insert(tetKDTree, e->center().data(), static_cast<void*>(e));

        }        
    }

}



bool Mesh::searchTetrahedronContain(Vector3D pos, Tetrahedron* &goalTet){
    bool rst = false;
    kdres *set = kd_nearest_range(tetKDTree, pos.data(), maxSizing);
    while (!kd_res_end(set)){
        Tetrahedron *tet = static_cast<Tetrahedron *>(kd_res_item_data(set));
        if (tet->boundingBox.contain(pos) && tet->contain(pos)){
            goalTet = tet;
            rst = true;
            break;
        }
        kd_res_next(set);
    }
    kd_res_free(set);
    return rst;

    // bool findGoal = false;
    // auto rtreeCallBack =
    // [&goalTet, &pos, &findGoal]
    // (Tetrahedron* const &tet)
    // {
    //     bool rst = true;

    //     if (tet->contain(pos)){
    //         goalTet = tet;
    //         findGoal = true;
    //         rst = false;
    //     }
    //     return rst;
    // };

    // Vector3D max = pos + maxSizing;
    // Vector3D min = pos - maxSizing;
    // tetRTree.Search(min.data(), max.data(), rtreeCallBack);
    // return findGoal;
}

bool Mesh::searchTetrahedronContain(Vector3D pos,  Tetrahedron* &goalTet, std::array<double, 4> &weights){

    bool rst = false;
    kdres *set = kd_nearest_range(tetKDTree, pos.data(), maxSizing);
    while (!kd_res_end(set)){
        Tetrahedron *tet = static_cast<Tetrahedron *>(kd_res_item_data(set));
        if (tet->boundingBox.contain(pos) && tet->contain(pos, weights)){
            goalTet = tet;
            rst = true;
            break;
        }
        kd_res_next(set);
    }
    kd_res_free(set);
    return rst;  

    // bool findGoal = false;
    // auto rtreeCallBack =
    // [&goalTet, &pos, &findGoal, &weights]
    // (Tetrahedron* const &tet)
    // {
    //     bool rst = true;

    //     if (tet->contain(pos, weights)){
    //         goalTet = tet;
    //         findGoal = true;
    //         rst = false;
    //     }
    //     return rst;
    // };

    // Vector3D max = pos + maxSizing;
    // Vector3D min = pos - maxSizing;
    // tetRTree.Search(min.data(), max.data(), rtreeCallBack);
    // return findGoal;    
}

void Mesh::clone(const Mesh &aMesh){
    std::unordered_map<Node*, Node *> originNodesMap;
    std::unordered_map<Tetrahedron*, Tetrahedron*> originElementsMap;
    for(auto n: aMesh.nodes){
        Node *aNode = new Node();
        *aNode = *n;
        originNodesMap[n] = aNode;
        this->nodes.push_back(aNode);
    }

    for(auto e: aMesh.tetrahedrons){
        Tetrahedron *aTet = new Tetrahedron();
        *aTet = *e;
        originElementsMap[e] = aTet;
        this->tetrahedrons.push_back(aTet);
        for(auto &n: aTet->nodes){
            n = originNodesMap[n];
        }
    }

    for(auto e: aMesh.tetrahedrons){
        for(auto &adj: e->adjacentTetrahedrons){
            if (adj){
                adj = originElementsMap[adj];
            }
        }
    }
}




void Mesh::checkBooleanRemove(Mesh &anotherMesh, int outerLayers){

    AABBox anotherBox;
    for(auto n: anotherMesh.nodes){
        anotherBox.insert(n->pos);
    }

    for(auto e: tetrahedrons){
        e->edit = 0;
        for(auto n: e->nodes){
            if(anotherBox.contain(n->pos)){
                Tetrahedron *goalTet;
                if (anotherMesh.searchTetrahedronContain(n->pos, goalTet)){
                    e->edit=-1;
                }
            }

        }
    }

    if (outerLayers>0){
        std::vector<Tetrahedron *> removeTets;
        for(auto n: nodes){
            n->edit = 0;
        }

        for(auto e: tetrahedrons){
            if(e->edit==-1){
                removeTets.push_back(e);
            }
        }
        for(auto e: removeTets){
                for(auto n: e->nodes){
                    n->edit = -1;
                }
        }

        std::vector<Tetrahedron *> freshTets;
        for(int i=0; i<outerLayers; i++){

            for(auto e:tetrahedrons){
                if (e->edit == 0){
                    for(auto n: e->nodes){
                        if (n->edit == -1){
                            e->edit = -1;
                            freshTets.push_back(e);
                            break;
                        }
                    }
                }
            }

            for(auto e: freshTets){
                for(auto n: e->nodes){
                    n->edit = -1;
                }
            }

            if(freshTets.empty()){
                break;
            }    
            freshTets.clear();
        }
    }
} 

void Mesh::loadNodeValues(const std::string &filePath){
    std::ifstream file(filePath);
    if (file.is_open()){
        while(file){
            std::string line;
            std::getline(file, line);
            std::stringstream lineStream(line);
            std::string keyString;
            lineStream >> keyString;
            if(keyString == "scalar"){
                std::string name;
                lineStream >> name;
                scalarValueNames.push_back(name);
                for(int i=0; i<nodes.size(); i++){
                    std::string line;
                    std::getline(file, line);
                    std::stringstream lineStream(line);
                    double value;
                    lineStream >> value;
                    nodes[i]->scalarValues.push_back(value); 
                }                
            }
            else if (keyString == "vector"){
                std::string name;
                lineStream >> name;
                vectorValueNames.push_back(name);
                for(int i=0; i<nodes.size(); i++){
                    std::string line;
                    std::getline(file, line);
                    std::stringstream lineStream(line);
                    double u, v, w;
                    lineStream >> u >> v>> w;
                    Vector3D vec(u, v, w);
                    nodes[i]->vectorValues.push_back(vec); 
                }    
            }
            
        }

        file.close();
    }
}



void Mesh::loadMESH(const std::string &filePath){
    std::ifstream file(filePath);
    if (file.is_open()){
        while (file)
        {
            std::string line;
            std::string keyString;
            std::getline(file, line);
            std::stringstream lineStream(line);
            lineStream >> keyString;

            if (keyString == "Vertices"){
                std::getline(file, line);
                std::stringstream lineStream(line);
                int nv;
                lineStream >> nv;
                for(int i=0; i<nv; i++){
                    std::getline(file, line);
                    std::stringstream vertexLineStream(line);
                    double x, y, z;
                    int label;
                    vertexLineStream >> x >> y >> z >> label;
                    Node *n = new Node(x, y, z);
                    n->label = label;
                    nodes.push_back(n);
                }
            }
            else if (keyString == "Tetrahedra"){
                std::getline(file, line);
                std::stringstream lineStream(line);
                int nt;
                lineStream >> nt;
                for(int i=0; i<nt;  i++){
                    std::getline(file, line);
                    std::stringstream tetLineStream(line);
                    int n0, n1, n2, n3, label;
                    tetLineStream >> n0 >> n1 >> n2 >> n3 >> label;
                    Tetrahedron *tet = new Tetrahedron(nodes[n0-1], nodes[n1-1], nodes[n2-1], nodes[n3-1]);
                    tet->label = label;
                    tetrahedrons.push_back(tet);
                }
            }
            else if (keyString == "End"){
                break;
            }

        }
        rebuildIndices();
        file.close();
    }
}

void Mesh::loadTETGENIO(tetgenio &in, bool withLabel){
    int numNodes = in.numberofpoints;
    int numTets = in.numberoftetrahedra;
    for(int i=0; i<numNodes; i++){
        Node *n = new Node(&(in.pointlist[i*3]));

        if (withLabel){
            n->label = in.pointmarkerlist[i];
        }
        
        nodes.push_back(n);
    }

    int first = in.firstnumber;
    for(int i=0; i<numTets; i++){
        int base = i*4;
        Tetrahedron *tet = new Tetrahedron(nodes[in.tetrahedronlist[base]-first]
        , nodes[in.tetrahedronlist[base+1] - first]
        , nodes[in.tetrahedronlist[base+2] - first]
        , nodes[in.tetrahedronlist[base+3] - first]);
        tetrahedrons.push_back(tet);
    }

    rebuildIndices();
}


void Mesh::exportMESH(const std::string &filePath){
    std::ofstream file(filePath);
    file << "MeshVersionFormatted 2\n";
    file << "Dimension\n         3\n";

    file << "Vertices\n";
    file << nodes.size() << "\n";
    for(auto n:nodes){
        file << n->pos[0] 
        << "  " << n->pos[1] 
        << "  " << n->pos[2]
        << "  " << n->label <<"\n";
    }

    file << "Tetrahedra\n";
    file << tetrahedrons.size() << "\n";
    for(auto e: tetrahedrons){
        file << e->nodes[0]->index + 1
        << "  " << e->nodes[1]->index + 1
        << "  " << e->nodes[2]->index + 1
        << "  " << e->nodes[3]->index + 1
        << "  " << e->label << "\n"; 
    }

    file << "End\n";
    file.close();
}


void Mesh::readyForCavityBasedOP(bool toRebuildAdjacency, bool toReadyForSpatialSearch, bool toEstimateSizing){
    if(toRebuildAdjacency){
        rebuildTetrahedronsAdjacency();
    }
    
    if(toEstimateSizing){
        estimateSizing();
    }

    if (toReadyForSpatialSearch){
        readyForSpatialSearch(true, true, !toEstimateSizing);
    }

    for(auto e: tetrahedrons){
        e->generateCircumSphere();
    }
}

void Mesh::CavityBasedInsertNode(Tetrahedron *tet, Node *insertNode){

    std::vector<TriangleFacet> sFacets;
    std::vector<Tetrahedron *> rmTets;
    rmTets.push_back(tet);
    for(int i=0; i<rmTets.size(); i++){
        Tetrahedron *e = rmTets[i];
        for(auto adje: e->adjacentTetrahedrons){
            if(e 
            && !e->fixed 
            && std::find(rmTets.begin(), rmTets.end(), e) == rmTets.end()
            && e->circumsphere.contain(insertNode->pos) ){
                rmTets.push_back(e);
            }
        }
    }
    extractBorderFacets(rmTets, sFacets);
    insertNode->index = nodes.size();
    nodes.push_back(insertNode);

    for(auto &f: sFacets){
        Tetrahedron *insertTet = new Tetrahedron(f.sNodes[0], f.sNodes[1], f.sNodes[2], insertNode);
        insertTet->generateBoundingBox();
        insertTet->generateCircumSphere();
        insertTet->index = tetrahedrons.size();
        tetrahedrons.push_back(insertTet);
    }
    
    deleteLargeScaleTetrahedronsPermanently(rmTets);
    rebuildTetrahedronsAdjacency();
}

void extractBorderNodes(std::vector<Tetrahedron *> &tets, std::vector<Node *> &borderNodes){
    std::set<Node *> nodeSet;
    HashFacetTable facetTable;
    for(auto e: tets){
        for(int i=0; i<4; i++){
            TriangleFacet keyFacet = e->facet(i, true);
            TriangleFacet goalFacet;
            if (facetTable.search(keyFacet, goalFacet)){
                facetTable.remove(goalFacet);
            }
            else{
                facetTable.insert(keyFacet);
            }
        }
    }

    for(auto kv:facetTable.columns){
        for(auto f: kv.second){
            for(auto n: f.sNodes){
                nodeSet.insert(n);
            }
        }
    }

    borderNodes.insert(borderNodes.end(), nodeSet.begin(), nodeSet.end());

    
}

void extractBorderFacets(std::vector<Tetrahedron *> &tets, std::vector<TriangleFacet> &borderFacets){
    HashFacetTable facetTable;
    for(auto e: tets){
        for(int i=0; i<4; i++){
            TriangleFacet keyFacet = e->facet(i, true);
            TriangleFacet goalFacet;
            if (facetTable.search(keyFacet, goalFacet)){
                facetTable.remove(goalFacet);
            }
            else{
                facetTable.insert(keyFacet);
            }
        }
    }

    for(auto kv:facetTable.columns){
        for(auto f: kv.second){
            borderFacets.push_back(f);
        }
    }    
}


void extratctBorder(std::vector<Tetrahedron *> &tets, std::vector<Node *> &borderNodes, std::vector<TriangleFacet> &borderFacets){
    std::set<Node *> nodeSet; 
    HashFacetTable facetTable;
    for(auto e: tets){
        for(int i=0; i<4; i++){
            TriangleFacet keyFacet = e->facet(i, true);
            TriangleFacet goalFacet;
            if (facetTable.search(keyFacet, goalFacet)){
                facetTable.remove(goalFacet);
            }
            else{
                facetTable.insert(keyFacet);
            }
        }
    }

    for(auto kv:facetTable.columns){
        for(auto f: kv.second){
            borderFacets.push_back(f);
            for(auto n: f.sNodes){
                nodeSet.insert(n);
            }
        }
    }
    borderNodes.insert(borderNodes.end(), nodeSet.begin(), nodeSet.end());
}

void transportNodesToTETGENIO(const std::vector<Node *> &sNodes, tetgenio &out){
    out.numberofpoints = sNodes.size();
    out.pointlist = new double[out.numberofpoints*3];
    out.pointmarkerlist = new int[out.numberofpoints];
    for(int i=0; i<out.numberofpoints; i++){
        for(int j=0; j<3; j++){
            out.pointlist[3*i+j] = sNodes[i]->pos[j]; 
        }
        out.pointmarkerlist[i] = sNodes[i]->label;
    }
}

void transportFacetsToTETGENIO(std::vector<Node *> &sNodes, std::vector<TriangleFacet> &facets, std::vector<Vector3D> &holes, tetgenio &out){
    
    std::unordered_map<Node *, int> nodeIndices;
    for(int i=0; i<sNodes.size(); i++){
        nodeIndices[sNodes[i]] = i;
    }
    
    out.firstnumber = 0;
    out.numberofpoints = sNodes.size();
    out.pointlist = new double[out.numberofpoints * 3];
    out.numberofpointmtrs = 1;
    out.pointmtrlist = new double[out.numberofpoints];
    out.pointmarkerlist = new int[out.numberofpoints];
    for(int i=0; i<out.numberofpoints; i++){
        int base =i*3;
        out.pointlist[base] = sNodes[i]->pos[0];
        out.pointlist[base+1] = sNodes[i]->pos[1];
        out.pointlist[base+2] = sNodes[i]->pos[2];
        out.pointmarkerlist[i] = sNodes[i]->label;
        out.pointmtrlist[i] = sNodes[i]->sizing;

    }

    out.numberoffacets = facets.size();
    out.facetlist = new tetgenio::facet[out.numberoffacets];
    out.facetmarkerlist = new int[out.numberoffacets];
    for(int i=0; i<out.numberoffacets; i++){
        tetgenio::facet &f = out.facetlist[i];
        f.numberofpolygons = 1;
        f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
        f.numberofholes = 0;
        f.holelist = NULL;
        
        tetgenio::polygon &p = f.polygonlist[0];
        p.numberofvertices = 3;
        p.vertexlist = new int[p.numberofvertices];
        p.vertexlist[0] = nodeIndices[facets[i].sNodes[0]];
        p.vertexlist[1] = nodeIndices[facets[i].sNodes[1]];
        p.vertexlist[2] = nodeIndices[facets[i].sNodes[2]];
    }

    out.numberofholes = holes.size();
    out.holelist = new double[out.numberofholes*3];
    for(int i=0; i<out.numberofholes; i++){
        int base = i*3;
        out.holelist[base] = holes[i][0];
        out.holelist[base+1] = holes[i][1];
        out.holelist[base+2] = holes[i][2];
    }
}


void instructTetrahedronConnectByTETGENIO(std::vector<Node *> &nodes, tetgenio &in, std::vector<Tetrahedron *> &tets){
    int numTets = in.numberoftetrahedra;
    for(int i=0; i<numTets; i++){
        int base = i*4;
        int firstNumber = in.firstnumber;
        int n0 = in.tetrahedronlist[base]-firstNumber;
        int n1 = in.tetrahedronlist[base+1]-firstNumber;
        int n2 = in.tetrahedronlist[base+2]-firstNumber;
        int n3 = in.tetrahedronlist[base+3]-firstNumber;
        Tetrahedron *tet = new Tetrahedron(nodes[n0], nodes[n1], nodes[n2], nodes[n3]);
        tets.push_back(tet);
    }
}



void Mesh::estimateSizing(){
    //TODO accelerate
    maxSizing = std::numeric_limits<double>::min();
    minSizing = std::numeric_limits<double>::max();
    for(auto n: nodes){
        n->sizing = std::numeric_limits<double>::max();
    }
    
    for(auto e: tetrahedrons){
        for(int i=0; i<4; i++){
            for(int j=i+1; j<4; j++){
                Node *n0 = e->nodes[i];
                Node *n1 = e->nodes[j];
                double d = distance(n0->pos, n1->pos);
                if(d<n0->sizing){
                    n0->sizing = d;
                }

                if (d<n1->sizing){
                    n1->sizing = d;
                }

                if (d>maxSizing){
                    maxSizing = d;
                }

                if (d<minSizing){
                    minSizing = d;
                }
            }
        }
    }

}

void Mesh::exportNodeValues(const std::string &filePath){
    std::ofstream file(filePath);
    for(int i=0; i<scalarValueNames.size(); i++){
        file << "scalar " << scalarValueNames[i] << std::endl;    
        for (auto n: nodes){
            file << n->scalarValues[i] << std::endl;
        }
    }

    for(int i=0; i<vectorValueNames.size(); i++){
        file << "vector " << vectorValueNames[i] << std::endl;    
        for (auto n: nodes){
        file << n->vectorValues[i][0] << "  "<<n->vectorValues[i][1] << "  " << n->vectorValues[i][2] << std::endl;
        }
    }
    file.close();
}


void Mesh::interpolateNodeValuesForAnotherMesh(Mesh &anotherMesh){
    estimateSizing();
    readyForSpatialSearch(true);
    anotherMesh.scalarValueNames = this->scalarValueNames;
    anotherMesh.vectorValueNames = this->vectorValueNames;
    for(auto n: anotherMesh.nodes){
        kdres *set = kd_nearest_range(nodeKDTree, n->pos.data(), 1e-10);
        if (kd_res_size(set)){
            double pos[3];
            Node *goalNode = static_cast<Node*>(kd_res_item_data(set));
            n->vectorValues = goalNode->vectorValues;
            n->scalarValues = goalNode->scalarValues;
            n->edit = 0;
            kd_res_free(set);
            continue;
        }
        kd_res_free(set);


        Tetrahedron *goalTet;
        std::array<double, 4> weights;
        if (searchTetrahedronContain(n->pos, goalTet, weights)){
            std::vector<double> scalars;
            std::vector<Vector3D> vecs;
            goalTet->interpolateNodeValue(weights, scalars, vecs);
            n->scalarValues = scalars;
            n->vectorValues = vecs;
            n->edit=0;
        }
        else{
            n->edit=-1;
        }
        
    }


}

//TODO CavityBasedInsert
void Mesh::CavityBasedInsert(std::vector<Vector3D> &positions, bool toRestartReady){


    if(toRestartReady){
        readyForCavityBasedOP();
    }

    for(auto pos: positions){
        for(auto e:tetrahedrons){
            if (e->boundingBox.contain(pos) && e->contain(pos)){
                Node *n = new Node(pos);
                n->label=-1;
                CavityBasedInsertNode(e, n);
                break;
            }
        }
    }

    rebuildIndices();

}


