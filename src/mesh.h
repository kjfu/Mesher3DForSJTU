#pragma once
#include <array>
#include <vector>
#include <algorithm>
#include <set>
#include <iostream>
#include "aabbox.h"
#include "vector3d.h"
#include "RTree.h"
#include "tetgen.h"
#include "sphere.h"
#include "edge.h"
double SixTimesTetrahedronVolume(Vector3D v0, Vector3D v1, Vector3D v2, Vector3D v3);

class Tetrahedron;
class TriangleElement;
struct TriangleFacet;




class TriangleElement{
public: 
    std::array<Node *, 3> nodes;
    AABBox boundingBox;
    int index;
    int label;
    bool fixed = false;
    int edit = 0;
    std::array<std::vector<TriangleElement* >, 3> adjacentTriangles; 

    TriangleElement(){

    }

    TriangleElement(Node *n0, Node *n1, Node *n2){
        nodes[0] = n0; nodes[1] = n1; nodes[2] = n2;
    }



    Edge edge(int index, bool withNodes=false){
        int i0 = (index+1)%3;
        int i1 = (index+2)%3;
        Edge e(nodes[i0]->index, nodes[i1]->index);
        if (withNodes){
            e.sNodes[0]=nodes[i0];
            e.sNodes[1]=nodes[i1];
        }
        e.tri = this;
        e.localIndexOfTriangle = index;
        return e;
    }    
};

struct TriangleFacet
{
    Tetrahedron *tet = nullptr;
    int localIndex;
    std::array<int, 3> orderedNodeIndices;
    std::vector<Node *> sNodes;
    TriangleFacet() = default;
    TriangleFacet(int i0, int i1, int i2){
        orderedNodeIndices[0] = i0;
        orderedNodeIndices[1] = i1;
        orderedNodeIndices[2] = i2;
        std::sort(orderedNodeIndices.begin(), orderedNodeIndices.end());
    }
};


class Tetrahedron 
{   

public:
    std::array<Node*, 4> nodes;
    std::array<Tetrahedron*, 4> adjacentTetrahedrons;
    AABBox boundingBox;
    Sphere circumsphere;
    int index;
    int label=0;
    bool fixed = false;
    int edit = 0;
    Tetrahedron(){
        adjacentTetrahedrons.fill(nullptr);
    }
    Tetrahedron(Node *n0, Node *n1, Node *n2, Node *n3){
        nodes[0]=n0; nodes[1]=n1; nodes[2]=n2; nodes[3]=n3;
        adjacentTetrahedrons.fill(nullptr);
    }

    Vector3D center(){
        Vector3D pos;
        for(auto n: nodes){
            pos += n->pos;
        }
        pos /= 4.0;
        return pos;
    }



    TriangleFacet facet(int index, bool withNodes = false){
        static int facetIndices[4][3]=
        { {1, 3, 2}
        , {0, 2, 3}
        , {0, 3, 1}
        , {0, 1, 2}};
        TriangleFacet rst(nodes[facetIndices[index][0]]->index
                        , nodes[facetIndices[index][1]]->index
                        , nodes[facetIndices[index][2]]->index);
        if (withNodes){
            rst.sNodes.push_back(nodes[facetIndices[index][0]]);
            rst.sNodes.push_back(nodes[facetIndices[index][1]]);
            rst.sNodes.push_back(nodes[facetIndices[index][2]]);
        }
        rst.localIndex = index;
        rst.tet = this;
        return rst;
    }

    Node *nodeOnFacet(int facetIndex, int nodeIndex){
        static int facetIndices[4][3]=
        { {1, 3, 2}
        , {0, 2, 3}
        , {0, 3, 1}
        , {0, 1, 2}};
        return nodes[facetIndices[facetIndex][nodeIndex]];
    }

    bool removeAdjacentTetrahedron(Tetrahedron *tet){
        bool rst = false;
        for(int i=0; i<4; i++){
            if(adjacentTetrahedrons[i] == tet){
                adjacentTetrahedrons[i] = nullptr;
                rst = true;
                break;
            }
        }
        return rst;
    }
    bool dualFacet(Node* n, TriangleFacet &goalFacet){
        bool rst = true;
        int goalIndex = -1;
        for(int i=0; i<4; i++){
            if(nodes[i]==n){
                goalIndex = i;
                break;
            }
        }
        if (goalIndex>=0){
            goalFacet = facet(goalIndex);
        }
        else{
            rst = false;
        }
        return rst;
    }

    bool replaceAdjacentTetrahedron(Tetrahedron *origin, Tetrahedron *another){
        bool rst = false;
        for(int i=0; i<4; i++){
            if(adjacentTetrahedrons[i] == origin){
                adjacentTetrahedrons[i] = another;
                rst = true;
                break;
            }
        }
        return rst;        
    }


    bool contain(Vector3D pos){
        static int nodeIndices[4][3]=
        { {1, 3, 2}
        , {0, 2, 3}
        , {0, 3, 1}
        , {0, 1, 2}};
        bool rst = true;
        for(int i=0; i<4; i++){
            if (SixTimesTetrahedronVolume(nodes[nodeIndices[i][0]]->pos, nodes[nodeIndices[i][1]]->pos, nodes[nodeIndices[i][2]]->pos, pos)<0){
                rst = false;
                break;
            }
        }
        return rst;        
    }

    bool contain(Vector3D pos, std::array<double, 4> &weights){
        static int nodeIndices[4][3]=
        { {1, 3, 2}
        , {0, 2, 3}
        , {0, 3, 1}
        , {0, 1, 2}};
        bool rst = true;
        for(int i=0; i<4; i++){
            double weight = SixTimesTetrahedronVolume(nodes[nodeIndices[i][0]]->pos, nodes[nodeIndices[i][1]]->pos, nodes[nodeIndices[i][2]]->pos, pos);
            if (weight<-std::numeric_limits<double>::epsilon()){
                rst = false;
                break;
            }
            else{
                weights[i] = weight;
            }
        }
        return rst;        
    }

    double volume(){
        return SixTimesTetrahedronVolume(nodes[0]->pos, nodes[1]->pos, nodes[2]->pos, nodes[3]->pos)/6.0;
    }

    bool isValid(){
        return (SixTimesTetrahedronVolume(nodes[0]->pos, nodes[1]->pos, nodes[2]->pos, nodes[3]->pos)>0);
    }

    void interpolateNodeValue(Vector3D pos, std::vector<double> &scalars, std::vector<Vector3D> &vectors){

        double v0 = SixTimesTetrahedronVolume(nodes[1]->pos, nodes[3]->pos, nodes[2]->pos, pos);
        double v1 = SixTimesTetrahedronVolume(nodes[0]->pos, nodes[2]->pos, nodes[3]->pos, pos);
        double v2 = SixTimesTetrahedronVolume(nodes[0]->pos, nodes[3]->pos, nodes[1]->pos, pos);
        double v3 = SixTimesTetrahedronVolume(nodes[0]->pos, nodes[1]->pos, nodes[2]->pos, pos);
        double v= v0 + v1 + v2 + v3;
        std::array<double, 4> weights{{v0/v, v1/v, v2/v, v3/v}};
        int numScalars = nodes[0]->scalarValues.size();
        int numVectors = nodes[0]->vectorValues.size();
        for(int i=0; i<numScalars; i++){
            double value =0;
            for(int j=0; j<4; j++){
                value += weights[j] * nodes[i]->scalarValues[j];
            }
            scalars.push_back(value);
        }

        for(int i=0; i<numVectors; i++){
            Vector3D vec(0.0);
            for(int j=0; j<4; j++){
                vec += weights[j] * nodes[i]->vectorValues[j];
            }
            vectors.push_back(vec);
        }

    }

    void interpolateNodeValue(std::array<double, 4> &sWeights, std::vector<double> &scalars, std::vector<Vector3D> &vectors){
        std::array<double, 4> weights = sWeights;
        double sumWeights = 0;
        for(auto w: weights){
            sumWeights += w;
        }
        for(auto &w: weights){
            w /= sumWeights;
        }

        int numScalars = nodes[0]->scalarValues.size();
        int numVectors = nodes[0]->vectorValues.size();
        for(int i=0; i<numScalars; i++){
            double value =0;
            for(int j=0; j<4; j++){
                value += weights[j] * nodes[i]->scalarValues[j];
            }
            scalars.push_back(value);
        }

        for(int i=0; i<numVectors; i++){
            Vector3D vec(0.0);
            for(int j=0; j<4; j++){
                vec += weights[j] * nodes[i]->vectorValues[j];
            }
            vectors.push_back(vec);
        }
    }

    void generateBoundingBox(){
        boundingBox.reset();
        for(auto n: nodes){
            boundingBox.insert(n->pos.XYZ());
        }
    }

    void generateCircumSphere(){
        circumsphere.initialize(nodes[0]->pos
        , nodes[1]->pos
        , nodes[2]->pos
        , nodes[3]->pos);
    }

};


class Mesh{
public:
    std::vector<std::string> scalarValueNames;
    std::vector<std::string> vectorValueNames;
    std::vector<Node *> nodes;
    std::vector<Tetrahedron *> tetrahedrons;
    // RTree<Tetrahedron*, double, 3, double, 10000> tetRTree;
    struct kdtree *nodeKDTree = nullptr;
    struct kdtree *tetKDTree = nullptr;
    AABBox aabbox;
    Mesh(){}
    ~ Mesh(){
        for(int i=0; i<nodes.size(); i++){
            delete nodes[i];
        }
        for (int i = 0; i < tetrahedrons.size(); i++){
            delete tetrahedrons[i];
        }
    }
    void clone(const Mesh &aMesh);
    void splitElementFromCentre(Tetrahedron *tet); 
    void rebuildIndices();
    void rebuildTetrahedronsAdjacency();
    void extractBorderNodeIndicesWithLabels(std::vector<int> labels, std::set<int> &nodeIndices);
    void extractBorderNodes(std::vector<Node *> &sNodes);
    void extractBorder(std::vector<Node *> &sNodes, std::vector<TriangleFacet> &sFacets);
    void deleteLargeScaleTetrahedronsPermanently(std::vector<Tetrahedron *> &delElements);
    void deleteSmallScaleTetrahedronsPermanently(std::vector<Tetrahedron *> &delElements);
    void mergeMesh(Mesh &another, std::vector<Node *> &mergeNodes);// useless
    void mergeMesh(Mesh &another, double eps=std::numeric_limits<double>::epsilon());
    void mergeMeshLegacy(Mesh &anotherMesh, std::vector<Node *> &mergeNodes, bool deleteUselessNodes = true);
    void checkBooleanRemove(Mesh &anotherMesh, int outerLayers=0); 
    
    //
    double maxSizing = std::numeric_limits<double>::min();
    double minSizing = std::numeric_limits<double>::max();
    void estimateSizing();
    void interpolateNodeValuesForAnotherMesh(Mesh &anotherMesh);

    //cavity-based op
    void readyForCavityBasedOP(bool toRebuildAdjacency = true, bool toReadyForSpatialSearch = true, bool toEstimateSizing=true);
    void CavityBasedInsertNode(Tetrahedron *tet, Node *insertNode);
    void CavityBasedInsert(std::vector<Vector3D> &positions, bool toRestartReady=true);

    //Spatial search
    double searchRangeSize;
    void readyForSpatialSearch(bool toBuildTetKDTree=true, bool toBuildNodeKDTree = true, bool toEstimateSizing = true);
    bool searchTetrahedronContain(Vector3D pos,  Tetrahedron* &goalTet);
    bool searchTetrahedronContain(Vector3D pos,  Tetrahedron* &goalTet, std::array<double, 4> &weights);

    //IO
    void loadMESH(const std::string &filePath);
    void exportVTK(const std::string &filePath);
    void loadNodeValues(const std::string &filePath);
    void exportNodeValues(const std::string &filePath);
    void loadTETGENIO(tetgenio &in, bool withLabel = false);
    void exportMESH(const std::string &filePath);


    void rebuildAABBox();
};






void extractBorderNodes(std::vector<Tetrahedron *> &tets, std::vector<Node *> &borderNodes);
void extractBorderFacets(std::vector<Tetrahedron *> &tets, std::vector<TriangleFacet> &borderFacets);
void extratctBorder(std::vector<Tetrahedron *> &tets, std::vector<Node *> &borderNodes, std::vector<TriangleFacet> &borderFacets);
void transportVector3dsToTETGENIO(const std::vector<Vector3D> &vec3ds, tetgenio &out);
void transportNodesToTETGENIO(const std::vector<Node *> &sNodes, tetgenio &out);
void transportFacetsToTETGENIO(std::vector<Node *> &sNodes, std::vector<TriangleFacet> &facets, std::vector<Vector3D> &holes, tetgenio &out);
void instructTetrahedronConnectByTETGENIO(std::vector<Node *> &nodes, tetgenio &in, std::vector<Tetrahedron *> &tets);