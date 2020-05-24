#include <iostream>
#include <vector>
#include <chrono>
#include "CooMatrix.h"
#include "sputils.h"
#include "solvers.h"
#include "mgutils.h"
#include "Array3D.h"
#include "grid.h"
#include "rkutils.h"

int main()
{
//	const int Nx = 129;
//	const int Ny = 129;
//	const int Nz = 129;
//	int nLevels = 6;
//	int nx = Nx - 2;
//	int ny = Ny - 2;
//	int nz = Nz - 2;
//	int neq = nx * ny * nz;
//	std::vector<std::vector<int>> leveldata{ {129, 65, 33, 17, 9, 5},
//								             {129, 65, 33, 17, 9, 5},
//								             {129, 65, 33, 17, 9, 5} };
//	std::vector<CooMatrix> A = genAdata(leveldata);
//	VECTOR_double b = computeb(neq);
//	std::vector<double> b = linspace(0.0, 1.0, neq);
//	std::vector<double> x = zeros(neq);
//	std::vector<CooMatrix> prolongdata = genProlongData(leveldata);
//	std::vector<CooMatrix> restrictdata = genRestrictData(leveldata);
//	int N = 50;
//	CooMatrix A = computeA(N);
//	CooMatrix A2 = kronsum(A, A);
//	CooMatrix A3 = kronsum(A2, A);
//	int N3 = N * N * N;
//	std::vector<double> b = linspace(0.0, 1.0, N3);
//	std::vector<double> x = zeros(N3);
//	auto time1 = std::chrono::high_resolution_clock::now();
//	gseidel(A3, b, x);
//	multigrid(A, b, x, prolongdata, restrictdata);
//	auto time2 = std::chrono::high_resolution_clock::now();
//	std::chrono::duration<double> solvetime = time2 - time1;
//	std::cout << "Total time(s): " << solvetime.count() << "\n";
//	std::vector<double> x = { 1.5,2.0,3.2,4.5,5.0 };
//	std::vector<double> b = A * x;
//	print1D(b);
//	std::cout << A2 << std::endl;
//	printmat(A);
	/*
	int Nx = 5;
	int Ny = 5;
	int Nz = 5;
	std::vector<Array3D> x = { Array3D(Nx, Ny, Nz) , Array3D(Nx, Ny, Nz) , Array3D(Nx, Ny, Nz) };
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			for (int k = 0; k < Nz; ++k) {
				x[0](i, j, k) = 100 * i + 10 * j + k;
			}
		}
	}
	for (int i = 0; i < x[0].data.size(); ++i) {
		x[0].data[i] = 2.0 * x[0].data[i];
	}
	*/
//	std::vector<double> x(Nx * Ny * Nz, 5);
//	print1D(x[0].data);
//	print1D(x[1].data);
//	print1D(x[2].data);
//	print1D(std::vector<double>(x.data.begin(), x.data.end()));
	auto time1 = std::chrono::high_resolution_clock::now();
	double Re = 5000;
	int nt = 13000;
	int Nx = 128;
	int Ny = 128;
	int Nz = 128;
	double Lx = 1.0;
	double Ly = 1.0;
	double Lz = 1.0;
	grid Grid(Lx, Ly, Lz, Nx, Ny, Nz);
	//Initialize variables
	std::vector<Array3D> var = { Array3D(Grid.Nx, Grid.Ny + 1, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny + 1, Grid.Nz) ,
							   Array3D(Grid.Nx + 1, Grid.Ny + 1, Grid.Nz + 1) };
	//Apply boundary conditions
	equatebc(var, var, Grid);
	//Initiate iteration and time
	int iteration = 0;
	double time = 0.0;

	//begin time-loop
	while (iteration < nt) {
		double dt = timestep(var, Grid);
//		std::cout << dt << std::endl;
		time += dt;
		iteration += 1;
		std::vector<Array3D> varOld = { var[0], var[1], var[2] };
		std::vector<Array3D> E = { Array3D(Grid.Nx, Grid.Ny + 1, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny + 1, Grid.Nz) };
		std::vector<Array3D> q = { Array3D(Grid.Nx, Grid.Ny + 1, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny + 1, Grid.Nz) };
		//RK substep 1
		explicitTerm(E, var, Grid, Re);
		rkss(0.0, 1.0 / 3.0, 1.0 / 3.0, dt, var, q, E, Grid, Re);
		equatebc(var, var, Grid);
		//RK substep 2
		explicitTerm(E, var, Grid, Re);
		rkss(-5.0 / 9.0, 15.0 / 16.0, 5.0 / 12.0, dt, var, q, E, Grid, Re);
		equatebc(var, var, Grid);
		//RK substep 3
		explicitTerm(E, var, Grid, Re);
		rkss(-153.0 / 128.0, 8.0 / 15.0, 1.0 / 4.0, dt, var, q, E, Grid, Re);
		equatebc(var, var, Grid);
		double res = residual(var, varOld, E, Grid, dt, Re);
		auto time2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> solvetime = time2 - time1;
		std::cout << "Iteration = " << iteration << ", time = " << time << ", CPU Time(s): " << solvetime.count() << ", |Residual| = " << res << "\n";
	}
	
}