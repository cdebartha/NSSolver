#include "rkutils.h"
#include "mgutils.h"
#include "solvers.h"
#include <cmath>

void thomas(std::vector<double>& a, std::vector<double> b, std::vector<double>& c, std::vector<double> d, std::vector<double>& x)
{
	int n = d.size();
	for (int i = 1; i < n; ++i) {
		double w = a[i - 1] / b[i - 1];
		b[i] -= w * c[i - 1];
		d[i] -= w * d[i - 1];
	}
	x[n - 1] = d[n - 1] / b[n - 1];
	for (int i = n - 2; i > -1; --i) {
		x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
	}
}

double timestep(std::vector<Array3D>& var, grid& Grid)
{
	double varMax[3];
	double dt = DBL_MAX;
	double c = 1.2;
	for (int ivar = 0; ivar < 3; ++ivar) {
		varMax[ivar] = var[ivar].data[0];
		for (int idata = 1; idata < var[ivar].data.size(); ++idata) {
			if (var[ivar].data[idata] > varMax[ivar]) {
				varMax[ivar] = var[ivar].data[idata];
			}
		}
	}
	if (varMax[0] > 1e-8) {
		double dt_max = c * Grid.hx / varMax[0];
		if (dt_max < dt) {
			dt = dt_max;
		}
	}
	if (varMax[1] > 1e-8) {
		double dt_max = c * Grid.hy / varMax[1];
		if (dt_max < dt) {
			dt = dt_max;
		}
	}
	if (varMax[2] > 1e-8) {
		double dt_max = c * Grid.hz / varMax[2];
		if (dt_max < dt) {
			dt = dt_max;
		}
	}
	return dt;
}

void equatebc(std::vector<Array3D>& varStar, std::vector<Array3D>& var, grid& Grid)
{
	//No-slip at x-min and x-max
	for (int k = 0; k < Grid.Nz + 1; ++k) {
		for (int j = 0; j < Grid.Ny + 1; ++j) {
			varStar[0](0, j, k) = 0;
			varStar[0](Grid.Nx - 1, j, k) = 0;
			varStar[3](0, j, k) = var[3](1, j, k);
			varStar[3](Grid.Nx, j, k) = var[3](Grid.Nx - 1, j, k);
			if (j != Grid.Ny) {
				varStar[1](0, j, k) = -var[1](1, j, k);
				varStar[1](Grid.Nx, j, k) = -var[1](Grid.Nx - 1, j, k);
			}
			if (k != Grid.Nz) {
				varStar[2](0, j, k) = -var[2](1, j, k);
				varStar[2](Grid.Nx, j, k) = -var[2](Grid.Nx - 1, j, k);
			}
		}
	}

	//Periodic in y
	for (int k = 0; k < Grid.Nz + 1; ++k) {
		for (int i = 0; i < Grid.Nx + 1; ++i) {
			varStar[1](i, 0, k) = var[1](i, Grid.Ny - 2, k);
			varStar[1](i, Grid.Ny - 1, k) = var[1](i, 1, k);
			varStar[3](i, 0, k) = var[3](i, Grid.Ny - 1, k);
			varStar[3](i, Grid.Ny, k) = var[3](i, 1, k);
			if (i != Grid.Nx) {
				varStar[0](i, 0, k) = var[0](i, Grid.Ny - 1, k);
				varStar[0](i, Grid.Ny, k) = var[0](i, 1, k);
			}
			if (k != Grid.Nz) {
				varStar[2](i, 0, k) = var[2](i, Grid.Ny - 1, k);
				varStar[2](i, Grid.Ny, k) = var[2](i, 1, k);
			}
		}
	}

	//No-slip at z-min and all-fixed at z-max
	for (int j = 0; j < Grid.Ny + 1; ++j) {
		for (int i = 0; i < Grid.Nx + 1; ++i) {
			varStar[2](i, j, 0) = 0.0;
			varStar[2](i, j, Grid.Nz - 1) = 0.0;
			varStar[3](i, j, 0) = var[3](i, j, 1);
			varStar[3](i, j, Grid.Nz) = var[3](i, j, Grid.Nz - 1);
			if (i != Grid.Nx) {
				varStar[0](i, j, 0) = -var[0](i, j, 1);
				varStar[0](i, j, Grid.Nz) = 2 - var[0](i, j, Grid.Nz - 1);
			}
			if (j != Grid.Ny) {
				varStar[1](i, j, 0) = -var[1](i, j, 1);
				varStar[1](i, j, Grid.Nz) = -var[1](i, j, Grid.Nz - 1);
			}
		}
	}
}

void rhsPPE(Array3D& rhs, Array3D& u, Array3D& v, Array3D& w, double dt, grid& Grid)
{
	for (int k = 1; k < Grid.Nz; ++k) {
		for (int j = 1; j < Grid.Ny; ++j) {
			for (int i = 1; i < Grid.Nx; ++i) {
				rhs(i, j, k) = (1 / dt) * ((u(i, j, k) - u(i - 1, j, k)) / (Grid.xe[i] - Grid.xe[i - 1]) + 
					(v(i, j, k) - v(i, j - 1, k)) / (Grid.ye[j] - Grid.ye[j - 1]) + 
					(w(i, j, k) - w(i, j, k - 1)) / (Grid.ze[k] - Grid.ze[k - 1]));
			}
		}
	}
}

void computeDelP(std::vector<Array3D>& delp, Array3D& p, grid& Grid)
{
	for (int k = 1; k < Grid.Nz; ++k) {
		for (int j = 1; j < Grid.Ny; ++j) {
			for (int i = 1; i < Grid.Nx; ++i) {
				if (i != Grid.Nx - 1) {
					delp[0](i, j, k) = (p(i + 1, j, k) - p(i, j, k)) / (Grid.xc[i + 1] - Grid.xc[i]);
				}
				if (j != Grid.Ny - 1) {
					delp[1](i, j, k) = (p(i, j + 1, k) - p(i, j, k)) / (Grid.yc[j + 1] - Grid.yc[j]);
				}
				if (k != Grid.Nz - 1) {
					delp[2](i, j, k) = (p(i, j, k + 1) - p(i, j, k)) / (Grid.zc[k + 1] - Grid.zc[k]);
				}
			}
		}
	}
}
void generateb(std::vector<double>& b, Array3D& f, Array3D& phi, grid& Grid)
{
	int nx = Grid.Nx - 1;
	int ny = Grid.Ny - 1;
	int nz = Grid.Nz - 1;
	double h = Grid.hx; //Assuming uniform mesh
	for (int k = 1; k < Grid.Nz; ++k) {
		for (int j = 1; j < Grid.Ny; ++j) {
			for (int i = 1; i < Grid.Nx; ++i) {
				b[(i - 1) + nx * (j - 1) + nx * ny * (k - 1)] = h * h * f(i, j, k) - phi(i - 1, j, k) - phi(i + 1, j, k) - phi(i, j - 1, k) - phi(i, j + 1, k) - phi(i, j, k - 1) - phi(i, j, k + 1);
			}
		}
	}
}

void initguess(std::vector<double>& x, std::vector<double>& d)
{
	int n = d.size();
	std::vector<double> lower(n-1, 1.0);
	std::vector<double> upper(n-1, 1.0);
	std::vector<double> diag(n, -6.0);
	thomas(lower, diag, upper, d, x);
}

void solvePPE(Array3D& rhsp, Array3D& p, grid& Grid)
{
	int nx = Grid.Nx - 1;
	int ny = Grid.Ny - 1;
	int nz = Grid.Nz - 1;
	int neq = nx * ny * nz;
	std::vector<double> b(neq);
	std::vector<double> p1d(neq);
	generateb(b, rhsp, p, Grid);
	initguess(p1d, b);
	multigrid(Grid.A, b, p1d, Grid.prolongdata, Grid.restrictdata);
	for (int k = 1; k < Grid.Nz; ++k) {
		for (int j = 1; j < Grid.Ny; ++j) {
			for (int i = 1; i < Grid.Nx; ++i) {
				p(i, j, k) = p1d[(i - 1) + nx * (j - 1) + nx * ny * (k - 1)];
			}
		}
	}
}

void rkss(double c1, double c2, double c3, double delt, std::vector<Array3D>& var, std::vector<Array3D>& q, std::vector<Array3D>& E, grid& Grid, double Re)
{
	int nq = q.size();
	int nvar = var.size();
	std::vector<Array3D> rhs = { Array3D(Grid.Nx, Grid.Ny + 1, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny + 1, Grid.Nz) };
	std::vector<Array3D> varStar = { Array3D(Grid.Nx, Grid.Ny + 1, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny + 1, Grid.Nz) ,
							   Array3D(Grid.Nx + 1, Grid.Ny + 1, Grid.Nz + 1) };
	double deltrk = c3 * delt;
	for (int i = 0; i < nq; ++i) {
		for (int idata = 0; idata < q[i].data.size(); ++idata) {
			q[i].data[idata] = c1 * q[i].data[idata] + delt * E[i].data[idata];
			rhs[i].data[idata] = var[i].data[idata] - c2 * q[i].data[idata];
		}
	}
	equatebc(varStar, var, Grid);
	for (int iq = 0; iq < nq; ++iq) {
		int nx = q[iq].shape[0];
		int ny = q[iq].shape[1];
		int nz = q[iq].shape[2];
		double cBnd = deltrk / (Re * Grid.hz * Grid.hz);
		double cLower = -cBnd;
		double cUpper = -cBnd;
		double cDiag = 1 + 2 * cBnd;
		std::vector<double> lower(nz - 3, cLower);
		std::vector<double> upper(nz - 3, cUpper);
		std::vector<double> diag(nz - 2, cDiag);
		for (int i = 1; i < nx - 1; ++i) {
			for (int j = 1; j < ny - 1; ++j) {
				std::vector<double> rhsTemp(nz - 2);
				for (int k = 1; k < nz - 1; ++k) {
					rhsTemp[k - 1] = rhs[iq](i, j, k);
				}
				rhsTemp[0] += cBnd * varStar[iq](i, j, 0);
				rhsTemp[nz - 3] += cBnd * varStar[iq](i, j, nz - 1);
				std::vector<double> varTemp(nz - 2);
				thomas(lower, diag, upper, rhsTemp, varTemp);
				for (int k = 1; k < nz - 1; ++k) {
					varStar[iq](i, j, k) = varTemp[k - 1];
				}
			}
		}
	}
	Array3D rhsP(Grid.Nx + 1, Grid.Ny + 1, Grid.Nz + 1);
	rhsPPE(rhsP, varStar[0], varStar[1], varStar[2], deltrk, Grid);
	solvePPE(rhsP, varStar[nvar - 1], Grid);

	std::vector<Array3D> delP = { Array3D(Grid.Nx, Grid.Ny + 1, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny, Grid.Nz + 1) ,
							   Array3D(Grid.Nx + 1, Grid.Ny + 1, Grid.Nz) };
	computeDelP(delP, varStar[nvar - 1], Grid);

	for (int ivar = 0; ivar < 3; ++ivar) {
		for (int idata = 0; idata < var[ivar].data.size(); ++idata) {
			var[ivar].data[idata] = varStar[ivar].data[idata] - deltrk * delP[ivar].data[idata];
		}
	}
	for (int ivar = 3; ivar < nvar; ++ivar) {
		for (int idata = 0; idata < var[ivar].data.size(); ++idata) {
			var[ivar].data[idata] = varStar[ivar].data[idata];
		}
	}
}

void explicitTerm(std::vector<Array3D>& E, std::vector<Array3D>& var, grid& Grid, double Re)
{
	for (int k = 1; k < Grid.Nz; ++k) {
		for (int j = 1; j < Grid.Ny; ++j) {
			for (int i = 1; i < Grid.Nx - 1; ++i) {
				//U2x
				double uavg2 = 0.5 * (var[0](i + 1, j, k) + var[0](i, j, k));
				double uavg1 = 0.5 * (var[0](i, j, k) + var[0](i - 1, j, k));

				double U2x = (uavg2 * uavg2 - uavg1 * uavg1) / (Grid.xc[i + 1] - Grid.xc[i]);

				//UVy
				double vavg2 = 0.5 * (var[1](i, j, k) + var[1](i + 1, j, k));
				double vavg1 = 0.5 * (var[1](i, j - 1, k) + var[1](i + 1, j - 1, k));
				uavg2 = 0.5 * (var[0](i, j, k) + var[0](i, j + 1, k));
				uavg1 = 0.5 * (var[0](i, j, k) + var[0](i, j - 1, k));

				double UVy = (uavg2 * vavg2 - uavg1 * vavg1) / (Grid.ye[j] - Grid.ye[j - 1]);

				//UWz
				double wavg2 = 0.5 * (var[2](i, j, k) + var[2](i + 1, j, k));
				double wavg1 = 0.5 * (var[2](i, j, k - 1) + var[2](i + 1, j, k - 1));
				uavg2 = 0.5 * (var[0](i, j, k) + var[0](i, j, k + 1));
				uavg1 = 0.5 * (var[0](i, j, k) + var[0](i, j, k - 1));

				double UWz = (uavg2 * wavg2 - uavg1 * wavg1) / (Grid.ze[k] - Grid.ze[k - 1]);

				//Ux2
				double Ux2 = ((var[0](i + 1, j, k) - var[0](i, j, k)) / (Grid.xe[i + 1] - Grid.xe[i]) - 
					(var[0](i, j, k) - var[0](i - 1, j, k)) / (Grid.xe[i] - Grid.xe[i - 1])) /
					(Grid.xc[i + 1] - Grid.xc[i]);

				//Uy2
				double Uy2 = ((var[0](i, j + 1, k) - var[0](i, j, k)) / (Grid.yc[j + 1] - Grid.yc[j]) -
					(var[0](i, j, k) - var[0](i, j - 1, k)) / (Grid.yc[j] - Grid.yc[j - 1])) /
					(Grid.ye[j] - Grid.ye[j - 1]);

				E[0](i, j, k) = U2x + UVy + UWz - (1.0 / Re) * (Ux2 + Uy2);
			}
		}
	}
	for (int k = 1; k < Grid.Nz; ++k) {
		for (int j = 1; j < Grid.Ny - 1; ++j) {
			for (int i = 1; i < Grid.Nx; ++i) {
				//V2y
				double vavg2 = 0.5 * (var[1](i, j + 1, k) + var[1](i, j, k));
				double vavg1 = 0.5 * (var[1](i, j, k) + var[1](i, j - 1, k));

				double V2y = (vavg2 * vavg2 - vavg1 * vavg1) / (Grid.yc[j + 1] - Grid.yc[j]);

				//VUx
				double uavg2 = 0.5 * (var[0](i, j + 1, k) + var[0](i, j, k));
				double uavg1 = 0.5 * (var[0](i - 1, j, k) + var[0](i - 1, j + 1, k));
				vavg2 = 0.5 * (var[1](i, j, k) + var[1](i + 1, j, k));
				vavg1 = 0.5 * (var[1](i, j, k) + var[1](i - 1, j, k));

				double VUx = (uavg2 * vavg2 - uavg1 * vavg1) / (Grid.xe[i] - Grid.xe[i - 1]);

				//VWz
				double wavg2 = 0.5 * (var[2](i, j, k) + var[2](i, j + 1, k));
				double wavg1 = 0.5 * (var[2](i, j, k - 1) + var[2](i, j + 1, k - 1));
				vavg2 = 0.5 * (var[1](i, j, k) + var[1](i, j, k + 1));
				vavg1 = 0.5 * (var[1](i, j, k) + var[1](i, j, k - 1));

				double VWz = (vavg2 * wavg2 - vavg1 * wavg1) / (Grid.ze[k] - Grid.ze[k - 1]);

				//Vx2
				double Vx2 = ((var[1](i + 1, j, k) - var[1](i, j, k)) / (Grid.xc[i + 1] - Grid.xc[i]) -
					(var[1](i, j, k) - var[1](i - 1, j, k)) / (Grid.xc[i] - Grid.xc[i - 1])) /
					(Grid.xe[i] - Grid.xe[i - 1]);

				//Vy2
				double Vy2 = ((var[1](i, j + 1, k) - var[1](i, j, k)) / (Grid.ye[j + 1] - Grid.ye[j]) -
					(var[1](i, j, k) - var[1](i, j - 1, k)) / (Grid.ye[j] - Grid.ye[j - 1])) /
					(Grid.yc[j + 1] - Grid.yc[j]);

				E[1](i, j, k) = V2y + VUx + VWz - (1.0 / Re) * (Vx2 + Vy2);
			}
		}
	}

	for (int k = 1; k < Grid.Nz - 1; ++k) {
		for (int j = 1; j < Grid.Ny; ++j) {
			for (int i = 1; i < Grid.Nx; ++i) {
				//W2z
				double wavg2 = 0.5 * (var[2](i, j, k + 1) + var[2](i, j, k));
				double wavg1 = 0.5 * (var[2](i, j, k) + var[2](i, j, k - 1));

				double W2z = (wavg2 * wavg2 - wavg1 * wavg1) / (Grid.zc[k + 1] - Grid.zc[k]);

				//WVy
				double vavg2 = 0.5 * (var[1](i, j, k) + var[1](i, j, k + 1));
				double vavg1 = 0.5 * (var[1](i, j - 1, k) + var[1](i, j - 1, k + 1));
				wavg2 = 0.5 * (var[2](i, j, k) + var[2](i, j + 1, k));
				wavg1 = 0.5 * (var[2](i, j, k) + var[2](i, j - 1, k));

				double WVy = (wavg2 * vavg2 - wavg1 * vavg1) / (Grid.ye[j] - Grid.ye[j - 1]);

				//WUx
				wavg2 = 0.5 * (var[2](i + 1, j, k) + var[2](i, j, k));
				wavg1 = 0.5 * (var[2](i, j, k) + var[2](i - 1, j, k));
				double uavg2 = 0.5 * (var[0](i, j, k) + var[0](i, j, k + 1));
				double uavg1 = 0.5 * (var[0](i - 1, j, k) + var[0](i - 1, j, k + 1));

				double WUx = (uavg2 * wavg2 - uavg1 * wavg1) / (Grid.xe[i] - Grid.xe[i - 1]);

				//Wx2
				double Wx2 = ((var[2](i + 1, j, k) - var[2](i, j, k)) / (Grid.xc[i + 1] - Grid.xc[i]) -
					(var[2](i, j, k) - var[2](i - 1, j, k)) / (Grid.xc[i] - Grid.xc[i - 1])) /
					(Grid.xe[i] - Grid.xe[i - 1]);

				//Wy2
				double Wy2 = ((var[2](i, j + 1, k) - var[2](i, j, k)) / (Grid.yc[j + 1] - Grid.yc[j]) -
					(var[2](i, j, k) - var[2](i, j - 1, k)) / (Grid.yc[j] - Grid.yc[j - 1])) /
					(Grid.ye[j] - Grid.ye[j - 1]);

				E[2](i, j, k) = W2z + WVy + WUx - (1.0 / Re) * (Wx2 + Wy2);
			}
		}
	}
}

double residual(std::vector<Array3D>& var, std::vector<Array3D>& varOld, std::vector<Array3D>& E, grid& Grid, double delt, double Re)
{
	double normU = 0;
	double normV = 0;
	double normW = 0;
	double normdiv = 0;
	for (int k = 1; k < Grid.Nz; ++k) {
		for (int j = 1; j < Grid.Ny; ++j) {
			for (int i = 1; i < Grid.Nx - 1; ++i) {
				double val = (var[0](i, j, k) - varOld[0](i, j, k)) / delt + E[0](i, j, k) +
					(var[3](i + 1, j, k) - var[3](i, j, k)) / (Grid.xc[i + 1] - Grid.xc[i]) -
					(1.0 / Re) * ((var[0](i, j, k + 1) - var[0](i, j, k)) / (Grid.zc[k + 1] - Grid.zc[k]) -
						(var[0](i, j, k) - var[0](i, j, k - 1)) / (Grid.zc[k] - Grid.zc[k - 1])) /
					(Grid.ze[k] - Grid.ze[k - 1]);
				normU += val * val;
			}
		}
	}
	for (int k = 1; k < Grid.Nz; ++k) {
		for (int j = 1; j < Grid.Ny - 1; ++j) {
			for (int i = 1; i < Grid.Nx; ++i) {
				double val = (var[1](i, j, k) - varOld[1](i, j, k)) / delt + E[1](i, j, k) +
					(var[3](i, j + 1, k) - var[3](i, j, k)) / (Grid.yc[j + 1] - Grid.yc[j]) -
					(1.0 / Re) * ((var[1](i, j, k + 1) - var[1](i, j, k)) / (Grid.zc[k + 1] - Grid.zc[k]) -
						(var[1](i, j, k) - var[1](i, j, k - 1)) / (Grid.zc[k] - Grid.zc[k - 1])) /
					(Grid.ze[k] - Grid.ze[k - 1]);
				normV += val * val;
			}
		}
	}
	for (int k = 1; k < Grid.Nz - 1; ++k) {
		for (int j = 1; j < Grid.Ny; ++j) {
			for (int i = 1; i < Grid.Nx; ++i) {
				double val = (var[2](i, j, k) - varOld[2](i, j, k)) / delt + E[2](i, j, k) +
					(var[3](i, j, k + 1) - var[3](i, j, k)) / (Grid.zc[k + 1] - Grid.zc[k]) -
					(1.0 / Re) * ((var[2](i, j, k + 1) - var[2](i, j, k)) / (Grid.ze[k + 1] - Grid.ze[k]) -
						(var[2](i, j, k) - var[2](i, j, k - 1)) / (Grid.ze[k] - Grid.ze[k - 1])) /
					(Grid.zc[k + 1] - Grid.zc[k]);
				normW += val * val;
			}
		}
	}
	for (int k = 1; k < Grid.Nz; ++k) {
		for (int j = 1; j < Grid.Ny; ++j) {
			for (int i = 1; i < Grid.Nx; ++i) {
				double val = (var[0](i, j, k) - var[0](i - 1, j, k)) / (Grid.xe[i] - Grid.xe[i - 1]) +
					(var[1](i, j, k) - var[1](i, j - 1, k)) / (Grid.ye[j] - Grid.ye[j - 1]) +
					(var[2](i, j, k) - var[2](i, j, k - 1)) / (Grid.ze[k] - Grid.ze[k - 1]);
				normdiv += val * val;
			}
		}
	}
	double res = pow((normU + normV + normW + normdiv) /
		(var[0].data.size() + var[1].data.size() + var[2].data.size() + var[3].data.size()), 0.5);
	return res;
}