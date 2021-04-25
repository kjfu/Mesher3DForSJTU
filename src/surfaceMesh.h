#pragma once


#include <vector>
#include <string>
#include <unordered_map>
#include "tetgen.h"
#include "mesh.h"
#include "node.h"
#include "triangle.h"


enum PROJECTION_TYPE{
    XY_PLANE = 0,
    YZ_PLANE = 1,
    ZX_PLANE = 2,
};
class SurfaceMesh{
public:
    std::vector<Node *> nodes;
    std::vector<TriangleElement *> triangles;
    std::unordered_map<int, std::string> zoneNames;
    SurfaceMesh(){}
    ~SurfaceMesh(){

        for(int i=0; i<nodes.size(); i++){
            delete nodes[i];
        }
        for(int i=0; i<triangles.size(); i++){
            delete triangles[i];
        }     
    }


    void cloneSurfaceMesh(SurfaceMesh &another);

    void addSubSurfaceMesh(SurfaceMesh &another);//Ensure no intersect elements with current mesh


    void mergeSurfaceMesh(SurfaceMesh &another, double tolerance = std::numeric_limits<double>::epsilon());

    void rebuildTriangleAdjacency();

    // void deleteTriangles(std::vector<TriangleElement*> delTriangles);


    void checkNonManifoldEdges(std::vector<Edge> &nonManifoldEdges);



    void deleteTriangles(std::vector<TriangleElement *> &toDelTriangles);

    
    void projectTRIANGULATEIO(triangulateio &in, PROJECTION_TYPE projectionType, double offset);

    //
    void loadMESH(const std::string &filePath);
    void loadVTK(const std::string &filePath);
    void loadPLY(const std::string &filePath);
    void loadTETGENIO(tetgenio &in);
    void exportPLY(const std::string &filePath);
    void exportVTK(const std::string &filePath);
    void exportMESH(const std::string &filePath);
    void exportSU2(const std::string &filePath);
    void exportTETGENIO(tetgenio &out);
    void rebuildIndices();

    double maxSizing = std::numeric_limits<double>::min();
    double minSizing = std::numeric_limits<double>::max();
    void estimateSizing();
};

