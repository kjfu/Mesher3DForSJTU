#ifndef _MESHER3D_IO_
#define _MESHER3D_IO_

#include "tetgen.h"
#include <string>
#include <vector>
#include <array>
#include "vector3d.h"




void loadMesh(tetgenio *in, std::string filePath);


void loadREMESH(std::vector<int> &elements, std::vector<std::array<double,3>> &points, std::string filePath);

void saveAsMESH(tetgenio *out, std::string filePath);

void saveAsMESH(tetgenio *out, std::string filePath, std::vector<int> tetMarkers);


#endif
