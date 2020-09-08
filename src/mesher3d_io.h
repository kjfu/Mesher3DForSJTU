#ifndef _MESHER3D_IO_
#define _MESHER3D_IO_

#include "tetgen.h"
#include <string>
#include <vector>


void loadMESH(tetgenio *in, std::string filePath);

void saveAsMESH(tetgenio *out, std::string filePath);

void saveAsMESH(tetgenio *out, std::string filePath, std::vector<int> tetMarkers);



#endif
