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
#include "SubEntities.h"
#include "GroupEntities.h"
#include "common.h"
#include "tetrahedron.h"
#include "SpatialSearcher.h"

class Tetrahedron;
class TriangleElement;
struct TriangleFacet;
class SurfaceMesh;



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

// struct TriangleFacet
// {
//     Tetrahedron *tet = nullptr;
//     int localIndex;
//     std::array<int, 3> orderedNodeIndices;
//     std::vector<Node *> sNodes;
//     TriangleFacet() = default;
//     TriangleFacet(int i0, int i1, int i2){
//         orderedNodeIndices[0] = i0;
//         orderedNodeIndices[1] = i1;
//         orderedNodeIndices[2] = i2;
//         std::sort(orderedNodeIndices.begin(), orderedNodeIndices.end());
//     }

//     bool testIntersection(TriangleFacet &another){
//         tetgenmesh tet;
//         return tet.tri_tri_inter(this->sNodes[0]->pos.data(), this->sNodes[1]->pos.data(), this->sNodes[2]->pos.data()
//         ,another.sNodes[0]->pos.data(), another.sNodes[1]->pos.data(), another.sNodes[2]->pos.data());
//     }
    
// };
    inline bool testIntersection(SubTriangle &my, SubTriangle &another){
        tetgenmesh tet;
        return tet.tri_tri_inter(my.forms[0]->pos.data(), my.forms[1]->pos.data(), my.forms[2]->pos.data()
        ,another.forms[0]->pos.data(), another.forms[1]->pos.data(), another.forms[2]->pos.data());
    }




class Mesh{
public:
    std::vector<std::string> scalarValueNames;
    std::vector<std::string> vectorValueNames;
    std::vector<Node *> nodes;
    std::vector<Tetrahedron *> tetrahedrons;

    // RTree<Tetrahedron*, double, 3, double, 10000> tetRTree;
    struct kdtree *nodeKDTree = nullptr;
    SpatialSearcher aSearcher;


    // struct kdtree *tetKDTree = nullptr;
    AABBox aabbox;
    Mesh(){
        aSearcher.mesh=this;
    }
    ~ Mesh(){
        for(int i=0; i<nodes.size(); i++){
            delete nodes[i];
        }
        for (int i = 0; i < tetrahedrons.size(); i++){
            delete tetrahedrons[i];
        }
    }


    Node* addNode(Vector3D vec){
        Node* n = new Node(vec);
        nodes.push_back(n);
        return n;
    }

    Tetrahedron* addTetrahedron(Node* n0, Node* n1, Node* n2, Node* n3){
        Tetrahedron *tet = new Tetrahedron(n0, n1, n2, n3);
        tetrahedrons.push_back(tet);
        return tet;
    }

    // VolumeBall getVolumeBall(Tetrahedron *tet, int iLocal){
    //     VolumeBall aBall;
    //     Node *keyNode = tet->nodes[iLocal];
    //     std::vector<Node *> nodes;
    //     return aBall;
    // }

    bool getVolumeShell(Tetrahedron *tet, int iLocal, VolumeShell &aShell){
        aShell.closed= true;
        aShell.tetrahedronEdges.clear();
        aShell.nodes.clear();
        SubEdge keyEdge = tet->getSubEdge(iLocal);
        aShell.tetrahedronEdges.push_back(keyEdge);

        int iFacet0 = TetrahedronFacetShareEdge[iLocal][0];
        int iFacet1 = TetrahedronFacetShareEdge[iLocal][1];
        Tetrahedron *prev = tet;
        Tetrahedron *nxt = tet->adjacentTetrahedrons[iFacet0];
        Node *ringVertex = tet->nodes[iFacet0];
        aShell.nodes.push_back(ringVertex);
        while(nxt!=nullptr && nxt!=tet){
            int subLocal = nxt->getLocalIndex(keyEdge);
            SubEdge sub = nxt->getSubEdge(subLocal);
            aShell.tetrahedronEdges.push_back(sub);
            int iiFacet0 = TetrahedronFacetShareEdge[sub.iLocal][0];
            int iiFacet1 = TetrahedronFacetShareEdge[sub.iLocal][1];
            Tetrahedron *tmp = nxt->adjacentTetrahedrons[iiFacet0];
            Node *tmpV = nxt->nodes[iiFacet0];
            if (tmp==prev){
                tmp = nxt->adjacentTetrahedrons[iiFacet1];
                tmpV = nxt->nodes[iiFacet1];
            }
            aShell.nodes.push_back(tmpV);
            prev = nxt;
            nxt = tmp;
        }

        if (nxt==nullptr){
            aShell.closed = false;
            nxt = tet->adjacentTetrahedrons[iFacet1];
            while (nxt!=nullptr && nxt!=tet){

                int subLocal = nxt->getLocalIndex(keyEdge);
                SubEdge sub = nxt->getSubEdge(subLocal);
                aShell.tetrahedronEdges.push_back(sub);
                int iiFacet0 = TetrahedronFacetShareEdge[sub.iLocal][0];
                int iiFacet1 = TetrahedronFacetShareEdge[sub.iLocal][1];
                Tetrahedron *tmp = nxt->adjacentTetrahedrons[iiFacet0];
                Node *tmpV = nxt->nodes[iiFacet0];
                if (tmp==prev){
                    tmp = nxt->adjacentTetrahedrons[iiFacet1];
                    tmpV = nxt->nodes[iiFacet1];
                }
                aShell.nodes.push_back(tmpV);
                prev = nxt;
                nxt = tmp;
                
            }
        }        
        
        return true;
    }




    void clone(const Mesh &aMesh);
    void splitElementFromCentre(Tetrahedron *tet); 
    void rebuildIndices();
    void rebuildTetrahedronsAdjacency();
    void extractBorderNodeIndicesWithLabels(std::vector<int> labels, std::set<int> &nodeIndices);
    void extractBorderNodes(std::vector<Node *> &sNodes);
    void extractBorder(std::vector<Node *> &sNodes, std::vector<SubTriangle> &sFacets);
    void extractBorder(SurfaceMesh &aSurface);    
    void deleteLargeScaleTetrahedronsPermanently(std::vector<Tetrahedron *> &delElements);
    void deleteSmallScaleTetrahedronsPermanently(std::vector<Tetrahedron *> &delElements);
    void mergeMesh(Mesh &another, std::vector<Node *> &mergeNodes);// useless
    void mergeMesh(Mesh &another, double eps=std::numeric_limits<double>::epsilon());
    void mergeMeshLegacy(Mesh &anotherMesh, std::vector<Node *> &mergeNodes, bool deleteUselessNodes = true);
    void checkBooleanRemove(Mesh &anotherMesh, int outerLayers=0); 
    void checkBooleanRemoveSpecial(Mesh &anotherMesh, int outerLayers=0); 
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
    bool searchTetrahedronIntersect(Tetrahedron *keyTet, Tetrahedron* &goalTet);

    //IO
    void loadMESH(const std::string &filePath);
    void exportVTK(const std::string &filePath);
    void loadNodeValues(const std::string &filePath);
    void exportNodeValues(const std::string &filePath);
    void loadTETGENIO(tetgenio &in, bool withLabel = false);
    void exportMESH(const std::string &filePath);

    //Perform staining
    void getSubRegionCenters(std::vector<Vector3D> &positions);

    void rebuildAABBox();
};






void extractBorderNodes(std::vector<Tetrahedron *> &tets, std::vector<Node *> &borderNodes);
void extractBorderFacets(std::vector<Tetrahedron *> &tets, std::vector<SubTriangle> &borderFacets);
void extratctBorder(std::vector<Tetrahedron *> &tets, std::vector<Node *> &borderNodes, std::vector<SubTriangle> &borderFacets);
void transportVector3dsToTETGENIO(const std::vector<Vector3D> &vec3ds, tetgenio &out);
void transportNodesToTETGENIO(const std::vector<Node *> &sNodes, tetgenio &out);
void transportFacetsToTETGENIO(std::vector<Node *> &sNodes, std::vector<SubTriangle> &facets, std::vector<Vector3D> &holes, tetgenio &out);
void instructTetrahedronConnectByTETGENIO(std::vector<Node *> &nodes, tetgenio &in, std::vector<Tetrahedron *> &tets);