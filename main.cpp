#include "tetgen.h"
#include "stdio.h"
#include <string>
#include "mesher3d_io.h"

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
            int cut = str.find_first_of(".poly");
            strcpy(inFilePath, str.substr(0, cut).c_str());
            getInFile = true;
        }
    }




    if (getInFile){
        in.load_poly(inFilePath);
    }
    else{
        exit(1);
    }


    in.numberofpointmtrs = 1;
    in.pointmtrlist = new REAL[in.numberofpoints];
    for(int i=0; i<in.numberofpoints; i++){

        in.pointmtrlist[i] = size;
    }



    tetrahedralize("pqm", &in, &out);

    saveAsMESH(&out, outFilePath);

    return 0;
}
