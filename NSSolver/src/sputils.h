#pragma once

#include "CooMatrix.h"

void printmat(CooMatrix& mat);

void print1D(const std::vector<double>& vec);

std::vector<double> linspace(double start, double end, int steps);

std::vector<double> zeros(int N);

std::vector<double> addvec(std::vector<double>& a, std::vector<double>& b);

void updatevec(std::vector<double>& x, const std::vector<double>& b);

CooMatrix eye(int N);

CooMatrix kron(CooMatrix& A, CooMatrix& B);

CooMatrix addmat(CooMatrix& A, CooMatrix& B);

CooMatrix kronsum(CooMatrix& A, CooMatrix& B);