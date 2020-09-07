#include "tetgen.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include "mesher3d_io.h"
#include "mesher3d_core.h"

int main(int argc, char *argv[]){
    tetgenio in, out;


    char inFilePath[1024];
	int choice = -1;
    std::string outFilePath = "out3d.mesh";
    std::string str;
    REAL size =0;
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
        else if (str.length()>2){
            
			if (str.find(".poly") != std::string::npos){
				size_t cut = str.find(".poly");
				strcpy(inFilePath, str.substr(0, cut).c_str());
				in.load_poly(inFilePath);
				choice = 1;
			}
			else if (str.find(".mesh") != std::string::npos){				
				loadMESH(&in, str);
				choice = 2;

			}
        }
    }




    if (choice == -1){
		fprintf(stderr, "No input file!\n");
		exit(1);
    }
	if (choice == 1){
			constrainedTetrahedralization(&in, &out, size);
	}
	else if (choice == 2){
		delaunayTetrahedralization(&in, &out, size);
	}
	

	
    saveAsMESH(&out, outFilePath);

    return 0;
}
