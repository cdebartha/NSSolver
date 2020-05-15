#include <iostream>
#include <vector>
#include <chrono>
#include "CooMatrix.h"
#include "sputils.h"
#include "solvers.h"
#include "mgutils.h"

int main()
{
	const int Nx = 129;
	const int Ny = 129;
	const int Nz = 129;
	int nLevels = 7;
	int nx = Nx - 2;
	int ny = Ny - 2;
	int nz = Nz - 2;
	int neq = nx * ny * nz;
	std::vector<std::vector<int>> leveldata{ {129, 65, 33, 17, 9, 5, 3},
								             {129, 65, 33, 17, 9, 5, 3},
								             {129, 65, 33, 17, 9, 5, 3} };
	std::vector<CooMatrix> A = genAdata(leveldata);
//	VECTOR_double b = computeb(neq);
	std::vector<double> b = linspace(0.0, 1.0, neq);
	std::vector<double> x = zeros(neq);
	std::vector<CooMatrix> prolongdata = genProlongData(leveldata);
	std::vector<CooMatrix> restrictdata = genRestrictData(leveldata);
//	int N = 50;
//	CooMatrix A = computeA(N);
//	CooMatrix A2 = kronsum(A, A);
//	CooMatrix A3 = kronsum(A2, A);
//	int N3 = N * N * N;
//	std::vector<double> b = linspace(0.0, 1.0, N3);
//	std::vector<double> x = zeros(N3);
	auto time1 = std::chrono::high_resolution_clock::now();
//	gseidel(A3, b, x);
	multigrid(A, b, x, prolongdata, restrictdata);
	auto time2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> solvetime = time2 - time1;
	std::cout << "Total time(s): " << solvetime.count() << "\n";
//	std::vector<double> x = { 1.5,2.0,3.2,4.5,5.0 };
//	std::vector<double> b = A * x;
//	print1D(b);
//	std::cout << A2 << std::endl;
//	printmat(A);
}