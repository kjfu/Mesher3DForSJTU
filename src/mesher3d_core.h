#pragma once
#include "tetgen.h"
#include <vector>
#include <string>





void delaunayTetrahedralization(tetgenio *in, tetgenio *out, REAL size, std::vector<int> &tetMarkers);

void constrainedTetrahedralization(tetgenio *in, tetgenio *out, REAL size);


void refineMesh(const std::string &fileInHead, const std::string &fileOutHead);


