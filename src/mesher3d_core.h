#pragma once
#include "tetgen.h"
#include <vector>
#include <string>





void delaunayTetrahedralization(tetgenio *in, tetgenio *out, REAL size, std::vector<int> &tetMarkers, bool beQuiet=false);

void constrainedTetrahedralization(tetgenio *in, tetgenio *out, REAL size, bool beQuiet=false);


void refineMesh(const std::string &fileInHead, const std::string &fileOutHead, bool beQuiet=false);


