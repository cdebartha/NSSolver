#pragma once

#include <vector>
#include "CooMatrix.h"
#include "sputils.h"

CooMatrix computeA(int N);

std::vector<CooMatrix> genAdata(std::vector<std::vector<int>>& leveldata);

CooMatrix computeP(int nf, int nc);

CooMatrix computeR(int nc, int nf);

std::vector<CooMatrix> genProlongData(std::vector<std::vector<int>>& leveldata);

std::vector<CooMatrix> genRestrictData(std::vector<std::vector<int>>& leveldata);

