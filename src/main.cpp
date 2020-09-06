#include "tetgen.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include "mesher3d_io.h"
#include "mesher3d_core.h"

int main(int argc, char *argv[]){
    tetgenio in, out;


    char inFilePath[1024];
	
    bool getInFile = false;
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
				size_t cut = str.find_last_of(".poly");
				strcpy(inFilePath, str.substr(0, cut).c_str());
				in.load_poly(inFilePath);
				getInFile = true;
			}
			else if (str.find(".mesh") != std::string::npos){				
				loadMESH(&in, str);
				getInFile = true;

			}
        }
    }




    if (!getInFile){
		fprintf(stderr, "No input file!\n");
		exit(1);
    }

	if (size > 1e-10){
			constrainedTetrahedralization(&in, &out, size);
	}
	else{
		delaunayTetrahedralization(&in, &out);
	}
	

	
    saveAsMESH(&out, outFilePath);

    return 0;
}
