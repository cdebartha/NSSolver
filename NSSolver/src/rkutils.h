#pragma once

#include <vector>
#include "Array3D.h"
#include "grid.h"

void thomas(std::vector<double>& a, std::vector<double> b, std::vector<double>& c, std::vector<double> d, std::vector<double>& x);

double timestep(std::vector<Array3D>& var, grid& Grid);

void equatebc(std::vector<Array3D>& varStar, std::vector<Array3D>& var, grid& Grid);

void explicitTerm(std::vector<Array3D>& E, std::vector<Array3D>& var, grid& Grid, double Re);

void rhsPPE(Array3D& rhs, Array3D& u, Array3D& v, Array3D& w, double dt, grid& Grid);

void computeDelP(std::vector<Array3D>& delp, Array3D& p, grid& Grid);

void generateb(std::vector<double>& b, Array3D& f, Array3D& phi, grid& Grid);

void initguess(std::vector<double>& x, std::vector<double>& d);

void solvePPE(Array3D& rhsp, Array3D& p, grid& Grid);

void rkss(double c1, double c2, double c3, double delt, std::vector<Array3D>& var, std::vector<Array3D>& q, std::vector<Array3D>& E, grid& Grid, double Re);

double residual(std::vector<Array3D>& var, std::vector<Array3D>& varOld, std::vector<Array3D>& E, grid& Grid, double delt, double Re);