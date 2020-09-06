#include "mesher3d_io.h"
#include <fstream>
#include <sstream>
#include <iostream>

void loadMESH(tetgenio *in, std::string filePath){
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
		}
	}
}


void saveAsMESH(tetgenio *out, std::string filePath){
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
                << "   0\n";
    }

    outfile << "End\n";
    outfile.close();

}

