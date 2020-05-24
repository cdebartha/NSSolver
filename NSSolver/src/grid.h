#pragma once

#include <vector>
#include "CooMatrix.h"

struct grid {
	double Lx;
	double Ly;
	double Lz;
	int Nx;
	int Ny;
	int Nz;
	double hx;
	double hy;
	double hz;
	std::vector<double> xc;
	std::vector<double> yc;
	std::vector<double> zc;
	std::vector<double> xe;
	std::vector<double> ye;
	std::vector<double> ze;
	std::vector<std::vector<int>> leveldata;
	std::vector<CooMatrix> A;
	std::vector<CooMatrix> prolongdata;
	std::vector<CooMatrix> restrictdata;

	grid(double Lx, double Ly, double Lz, int Nx, int Ny, int Nz);

	~grid() {};

};