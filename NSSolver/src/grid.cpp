#include "grid.h"
#include "sputils.h"
#include "mgutils.h"

grid::grid(double Lx, double Ly, double Lz, int Nx, int Ny, int Nz)
{
	this->Lx = Lx;
	this->Ly = Ly;
	this->Lz = Lz;
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
	this->hx = Lx / ((double)Nx - 1.0);
	this->hy = Ly / ((double)Ny - 1.0);
	this->hz = Lz / ((double)Nz - 1.0);
	this->xc = linspace(-hx / 2, Lx + hx / 2, Nx + 1);
	this->yc = linspace(-hy / 2, Ly + hy / 2, Ny + 1);
	this->zc = linspace(-hz / 2, Lz + hz / 2, Nz + 1);
	this->xe = linspace(0, Lx, Nx);
	this->ye = linspace(0, Ly, Ny);
	this->ze = linspace(0, Lz, Nz);
	int nx = Nx + 1;
	int ny = Ny + 1;
	int nz = Nz + 1;
	leveldata.resize(3);
	while (nx != 3 && ny != 3 && nz != 3) {
		this->leveldata[0].push_back(nx);
		this->leveldata[1].push_back(ny);
		this->leveldata[2].push_back(nz);
		nx = (int)((nx - 1) / 2) + 1;
		ny = (int)((ny - 1) / 2) + 1;
		nz = (int)((nz - 1) / 2) + 1;
	}
	this->A = genAdata(this->leveldata);
	this->prolongdata = genProlongData(this->leveldata);
	this->restrictdata = genRestrictData(this->leveldata);
}