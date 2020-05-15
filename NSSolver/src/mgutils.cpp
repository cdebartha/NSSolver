#include <iostream>
#include <vector>
#include "CooMatrix.h"
#include "utils.h"

CooMatrix computeA(int N)
{
	int nnz = N + 2 * (N - 1);
	std::vector<int> rowind;
	std::vector<int> colind;
	std::vector<double> val;
	rowind.reserve(nnz);
	colind.reserve(nnz);
	val.reserve(nnz);
	for (int i = 0; i < N; i++) {
		if (i == 0) {
			rowind.push_back(i);
			rowind.push_back(i);
			colind.push_back(i);
			colind.push_back(i + 1);
			val.push_back(-2.0);
			val.push_back(1.0);
		}
		else if (i == N - 1) {
			rowind.push_back(i);
			rowind.push_back(i);
			colind.push_back(i - 1);
			colind.push_back(i);
			val.push_back(1.0);
			val.push_back(-2.0);
		}
		else {
			rowind.push_back(i);
			rowind.push_back(i);
			rowind.push_back(i);
			colind.push_back(i - 1);
			colind.push_back(i);
			colind.push_back(i + 1);
			val.push_back(1.0);
			val.push_back(-2.0);
			val.push_back(1.0);
		}
	}
	return CooMatrix(N, N, nnz, rowind, colind, val);
}

std::vector<CooMatrix> genAdata(std::vector<std::vector<int>>& leveldata) 
{
	size_t nLevels = leveldata[0].size();
	std::vector<CooMatrix> A_mg;
	A_mg.reserve(nLevels);
	for (int i = 0; i < nLevels; ++i) {
		int Nx = leveldata[0][i] - 2;
		int Ny = leveldata[1][i] - 2;
		int Nz = leveldata[2][i] - 2;
		CooMatrix Ax = computeA(Nx);
		CooMatrix Ay = computeA(Ny);
		CooMatrix Az = computeA(Nz);
		CooMatrix Axy = kronsum(Ax, Ay);
		CooMatrix Axyz = kronsum(Axy, Az);
		A_mg.push_back(Axyz);
	}
	return A_mg;
}

CooMatrix computeP(int nf, int nc) 
{
	int nnz = 3 * nc;
	std::vector<int> rowind;
	std::vector<int> colind;
	std::vector<double> val;
	rowind.reserve(nnz);
	colind.reserve(nnz);
	val.reserve(nnz);

	rowind.push_back(0);
	for (int i = 1; i < nf - 2; i = i + 2) {
		rowind.push_back(i);
		rowind.push_back(i + 1);
		rowind.push_back(i + 1);
	}
	rowind.push_back(nf - 2);
	rowind.push_back(nf - 1);

	for (int i = 0; i < nc; ++i) {
		colind.push_back(i);
		colind.push_back(i);
		colind.push_back(i);
		val.push_back(0.5);
		val.push_back(1.0);
		val.push_back(0.5);
	}
	return CooMatrix(nf, nc, (int)val.size(), rowind, colind, val);
}

CooMatrix computeR(int nc, int nf) 
{
	int nnz = 3 * nc;
	std::vector<int> rowind;
	std::vector<int> colind;
	std::vector<double> val;
	rowind.reserve(nnz);
	colind.reserve(nnz);
	val.reserve(nnz);

	colind.push_back(0);
	for (int i = 1; i < nf - 2; i = i + 2) {
		colind.push_back(i);
		colind.push_back(i + 1);
		colind.push_back(i + 1);
	}
	colind.push_back(nf - 2);
	colind.push_back(nf - 1);

	for (int i = 0; i < nc; ++i) {
		rowind.push_back(i);
		rowind.push_back(i);
		rowind.push_back(i);
		val.push_back(0.25);
		val.push_back(0.5);
		val.push_back(0.25);
	}	
	return CooMatrix(nc, nf, (int)val.size(), rowind, colind, val);
}

std::vector<CooMatrix> genProlongData(std::vector<std::vector<int>>& leveldata) 
{
	size_t nLevels = leveldata[0].size();
	std::vector<CooMatrix> P_mg;
	P_mg.reserve(nLevels);
	for (size_t i = 1; i < nLevels; ++i) {
		int ncx = leveldata[0][i] - 2;
		int nfx = leveldata[0][i - 1] - 2;
		int ncy = leveldata[1][i] - 2;
		int nfy = leveldata[1][i - 1] - 2;
		int ncz = leveldata[2][i] - 2;
		int nfz = leveldata[2][i - 1] - 2;
		CooMatrix Px = computeP(nfx, ncx);
		CooMatrix Py = computeP(nfy, ncy);
		CooMatrix Pz = computeP(nfz, ncz);
		CooMatrix Pyx = kron(Py, Px);
		CooMatrix Pzyx = kron(Pz, Pyx);
		P_mg.push_back(Pzyx);
	}
	return P_mg;
}

std::vector<CooMatrix> genRestrictData(std::vector<std::vector<int>>& leveldata) 
{
	size_t nLevels = leveldata[0].size();
	std::vector<CooMatrix> R_mg;
	R_mg.reserve(nLevels);
	for (size_t i = 0; i < nLevels - 1; ++i) {
		int ncx = leveldata[0][i + 1] - 2;
		int nfx = leveldata[0][i] - 2;
		int ncy = leveldata[1][i + 1] - 2;
		int nfy = leveldata[1][i] - 2;
		int ncz = leveldata[2][i + 1] - 2;
		int nfz = leveldata[2][i] - 2;
		CooMatrix Rx = computeR(ncx, nfx);
		CooMatrix Ry = computeR(ncy, nfy);
		CooMatrix Rz = computeR(ncz, nfz);
		CooMatrix Ryx = kron(Ry, Rx);
		CooMatrix Rzyx = kron(Rz, Ryx);
		R_mg.push_back(Rzyx);
	}
	return R_mg;
}