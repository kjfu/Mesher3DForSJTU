#include "mesher3d_core.h"


void delaunayTetrahedralization(tetgenio *in, tetgenio *out){
	tetgenbehavior b;
	tetrahedralize(&b, in, out);
	
}

void constrainedTetrahedralization(tetgenio *in, tetgenio *out, REAL size){
	in->numberofpointmtrs = 1;
    in->pointmtrlist = new REAL[in->numberofpoints];
    for(int i=0; i<in->numberofpoints; i++){

        in->pointmtrlist[i] = size;
    }

	tetgenbehavior b;
	
	char command[] = "pqm";
	b.parse_commandline(command);
    tetrahedralize(&b, in, out);
}