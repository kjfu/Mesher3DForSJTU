#include "tetgen.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include <vector>
#include <array>
#include "mesher3d_io.h"
#include "mesher3d_core.h"
#include "vector3d.h"
#include <iostream>
#include <set>
#include "mesh.h"
#include <fstream>
#include <ctime>
#include <unordered_map>
int main(int argc, char *argv[]){
	// std::ofstream f("/home/kjfu/research/Mesher3DForSJTU/examples/refine_case/test3d.value");
	// f << "vector displacement\n";
	// for(int i=0 ;i<6117; i++){
	// 	f << rand()%100 << "  " << rand()%100 << "  " << rand()%10 << std::endl;
	// }
	// f << "scalar density\n";
	// for(int i=0 ;i<6117; i++){
	// 	f << rand()%100 << std::endl;
	// }
	// return 0;
    tetgenio in, out;
	Mesh goalMesh;
	Mesh backgroundMesh;

    char inFilePath[1024];
	int choice = -1;

    std::string str;
    REAL size = 0;

	std::string refineFileHeadIn;
    std::string refineFileHeadOut;

	std::string fileIn;
    std::string fileOut = "out3d.mesh";
	bool quiet = false;
	for(int i=1; i<argc; i++){
        str = std::string(argv[i]);

        if (str == "-s"){
            i++;
            size = atof(argv[i]);
        }
		else if(str=="-q"){
			quiet = true;
		}
        else if (str == "-o"){
            i++;
            fileOut = std::string(argv[i]);

        }
		else if (str == "-r"){
			i++;
			refineFileHeadIn = argv[i];
			refineFileHeadOut = std::string(argv[i])+"_out";
			choice = 3;
		}
		else if(str == "-p"){
			i++;
			fileIn = argv[i];
			choice = 4;
		}
        else if (str.length()>2){
            
			if (str.find(".poly") != std::string::npos){
				size_t cut = str.find(".poly");
				strcpy(inFilePath, str.substr(0, cut).c_str());
				in.load_poly(inFilePath);
				choice = choice>0?choice:1;
			}
			else if (str.find(".mesh") != std::string::npos){
				// size_t cut = str.find(".mesh");
				// strcpy(inFilePath, str.substr(0, cut).c_str());
				// in.load_medit(inFilePath, 1);
				if (choice==3||choice==4){
					continue;
				}				
				fileIn = argv[i];
				choice = choice>0?choice:2;

			}
        }
    }

	//TEST
	// choice = 2;
	// size = 15;
	// fileIn = "/home/kjfu/research/Mesher3DForSJTU/examples/bugCase/cp.mesh";
	// fileOut = "/home/kjfu/research/Mesher3DForSJTU/examples/bugCase/output.mesh";
	//TEST
    if (choice == -1){
		fprintf(stderr, "No input file!\n");
		exit(1);
    }
	if (choice == 1){
			constrainedTetrahedralization(&in, &out, size, quiet);    
			saveAsMESH(&out, fileOut);
	}
	else if (choice == 2){
		std::vector<int> tetMarkers;

		delaunayTetrahedralization(fileIn, fileOut, size, quiet);

	}
	else if (choice == 3){
		refineMesh(refineFileHeadIn, refineFileHeadOut, quiet);
	}
	else if(choice == 4){
		generatePeriodicBoundaryConditionMesh(fileIn, fileOut, size, quiet);
	}
	
    return 0;
}
