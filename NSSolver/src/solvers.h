#pragma once

#include "CooMatrix.h"
#include <vector>

double computeRMS(std::vector<double>& residual);

std::vector<double> computeResidual(CooMatrix& A, std::vector<double>& b, std::vector<double>& x);

void gauss_seidel_update(CooMatrix& A, std::vector<double>& b, std::vector<double>& x);

void gseidel(CooMatrix& A, std::vector<double>& b, std::vector<double>& x, int maxiters = 10000, double tol = 1e-6);

void mg_update(std::vector<CooMatrix>& A, std::vector<double>& b, std::vector<double>& x, std::vector<CooMatrix>& P, std::vector<CooMatrix>& R, size_t level, size_t nLevels);

void multigrid(std::vector<CooMatrix>& A, std::vector<double>& b, std::vector<double>& x, std::vector<CooMatrix>& P, std::vector<CooMatrix>& R, int maxiters = 10000, double tol = 1e-6);
