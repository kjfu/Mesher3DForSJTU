#pragma once
#include "tetgen.h"

void delaunayTetrahedralization(tetgenio *in, tetgenio *out);

void constrainedTetrahedralization(tetgenio *in, tetgenio *out, REAL size);