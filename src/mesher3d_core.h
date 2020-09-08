#pragma once
#include "tetgen.h"
#include <vector>
void delaunayTetrahedralization(tetgenio *in, tetgenio *out, REAL size, std::vector<int> &tetMarkers);

void constrainedTetrahedralization(tetgenio *in, tetgenio *out, REAL size);