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
    tetgenio in, out;
	Mesh goalMesh;
	Mesh backgroundMesh;

    char inFilePath[1024];
	int choice = -1;
    std::string outFilePath = "out3d.mesh";
    std::string str;
    REAL size = 0;

	std::string refineFileHeadIn;
    std::string refineFileHeadOut;

	for(int i=1; i<argc; i++){
        str = std::string(argv[i]);

        if (str == "-s"){
            i++;
            size = atof(argv[i]);
        }
        else if (str == "-o"){
            i++;
            outFilePath = std::string(argv[i]);

        }
		else if (str == "-r"){
			i++;
			refineFileHeadIn = argv[i];
			refineFileHeadOut = std::string(argv[i])+"_out";
			choice = 3;
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
				if(choice==3){
					continue;
				}				
				loadMesh(&in, str);
				choice = choice>0?choice:2;

			}
        }
    }


	choice = 3;

    if (choice == -1){
		fprintf(stderr, "No input file!\n");
		exit(1);
    }
	if (choice == 1){
			constrainedTetrahedralization(&in, &out, size);    
			saveAsMESH(&out, outFilePath);
	}
	else if (choice == 2){
		std::vector<int> tetMarkers;
		delaunayTetrahedralization(&in, &out, size, tetMarkers);
		saveAsMESH(&out, outFilePath, tetMarkers);
	}
	else if (choice == 3){
		refineMesh(refineFileHeadIn, refineFileHeadOut);
	}
	
    return 0;
}
