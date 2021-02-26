#include "mesher3d_io.h"
#include <fstream>
#include <sstream>
#include <iostream>

void loadMesh(tetgenio *in, std::string filePath){
	in->mesh_dim = 3;
	in->numberofpointattributes = 0;  // no point attribute.
	in->numberofpointmtrs = 0;
	in->firstnumber = 1;
	std::ifstream inFile(filePath);
	if (inFile.is_open()){

		while (inFile){		
			std::string line;
			std::string keystring;
			std::getline(inFile, line);
			std::stringstream lineStream(line);
			lineStream >> keystring;
			if (keystring == "Vertices"){
				std::getline(inFile, line);
				lineStream.clear();
				lineStream.str(line);				
				int nv;
				lineStream >> nv;
				in->numberofpoints = nv;
				in->pointlist = new REAL[nv*3];
				in->pointmarkerlist = new int[nv];
				for(int i=0; i<nv; i++){
					std::getline(inFile, line);
					lineStream.clear();
					lineStream.str(line);
					lineStream >> in->pointlist[3*i] >> in->pointlist[3*i+1] >> in->pointlist[3*i+2] >> in->pointmarkerlist[i];
				}
			}
            else if (keystring == "Tetrahedra"){
                std::getline(inFile, line);
                std::stringstream subLineStream(line);
                int nt;
                subLineStream >> nt;
                in->numberoftetrahedra = nt;
                in->numberofcorners = 4;
                // in->numberoftetrahedronattributes = 1;
                in->tetrahedronlist = new int[4*in->numberoftetrahedra];
                // in->tetrahedronattributelist = new REAL[in->numberoftetrahedra];
                for(int i=0; i<nt; i++){
                    std::getline(inFile, line);
                    std::stringstream subSubLineStream(line);
                    subSubLineStream >> in->tetrahedronlist[4*i] >> in->tetrahedronlist[4*i+1] >> in->tetrahedronlist[4*i+2] >> in->tetrahedronlist[4*i+3];
                    // in->tetrahedronattributelist[i] = 0;
                }
            }
		}
	}
}

void loadREMESH(std::vector<int> &elements, std::vector<std::array<double,3>> &points, std::string filePath){
    std::ifstream inFile(filePath);
    if (inFile.is_open()){
        while(inFile){
            std::string line;
            std::string keyString;
            std::getline(inFile, line);
            std::stringstream lineStream(line);
            lineStream >> keyString;
            if (keyString == "Append_points"){
                std::getline(inFile, line);
                lineStream.clear();
                lineStream.str(line);
                int numPoints;
                lineStream >> numPoints;
                for (int i=0; i<numPoints; i++){
                    std::getline(inFile, line);
                    std::stringstream subLineStream(line);
                    std::array<double,3> pos;
                    subLineStream >> pos[0] >> pos[1] >> pos[2];
                    points.push_back(pos);
                }
            }
            else if(keyString == "Refine_elements"){
                std::getline(inFile, line);
                lineStream.clear();
                lineStream.str(line);
                int numElements;
                lineStream >> numElements;
                for (int i=0; i<numElements; i++){
                    std::getline(inFile, line);
                    std::stringstream subLineStream(line);
                    int id;
                    subLineStream >> id;
                    elements.push_back(id);
                }
            }

        }
    }
    
}

void saveAsMESH(tetgenio *out, std::string filePath){
    std::ofstream outfile(filePath);
    outfile << "MeshVersionFormatted 2\n";
    outfile << "Dimension\n         3\n";
    outfile << "Vertices\n";
    outfile << out->numberofpoints << "\n";
    if(out->pointmarkerlist!=nullptr){
        for (int i=0; i<out->numberofpoints; i++){
            outfile << out->pointlist[i*3]
                    << "  " << out->pointlist[i*3+1]
                    << "  " << out->pointlist[i*3+2]
                    << "  " << out->pointmarkerlist[i] << "\n";
        }
    }
    else{
        for (int i=0; i<out->numberofpoints; i++){
            outfile << out->pointlist[i*3]
                    << "  " << out->pointlist[i*3+1]
                    << "  " << out->pointlist[i*3+2]
                    << "  " << 0 << "\n";
        }        
    }



    outfile << "Tetrahedra\n";
    outfile << out->numberoftetrahedra << "\n";
    for (int i=0; i<out->numberoftetrahedra; i++){
        outfile << out->tetrahedronlist[i*4]
                << "  " << out->tetrahedronlist[i*4+1]
                << "  " << out->tetrahedronlist[i*4+2]
                << "  " << out->tetrahedronlist[i*4+3]
                << "   0\n";
    }

    outfile << "End\n";
    outfile.close();

}

void saveAsMESH(tetgenio *out, std::string filePath, std::vector<int> tetMarkers){
    std::ofstream outfile(filePath);
    outfile << "MeshVersionFormatted 2\n";
    outfile << "Dimension\n         3\n";
    outfile << "Vertices\n";
    outfile << out->numberofpoints << "\n";
    for (int i=0; i<out->numberofpoints; i++){
        outfile << out->pointlist[i*3]
                << "  " << out->pointlist[i*3+1]
                << "  " << out->pointlist[i*3+2]
                << "  " << out->pointmarkerlist[i] << "\n";
    }

    outfile << "Tetrahedra\n";
    outfile << out->numberoftetrahedra << "\n";
    for (int i=0; i<out->numberoftetrahedra; i++){
        outfile << out->tetrahedronlist[i*4]
                << "  " << out->tetrahedronlist[i*4+1]
                << "  " << out->tetrahedronlist[i*4+2]
                << "  " << out->tetrahedronlist[i*4+3]
                << "  " << tetMarkers[i] << "\n";
    }

    outfile << "End\n";
    outfile.close();	
}

