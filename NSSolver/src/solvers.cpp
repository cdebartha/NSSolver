#include <iostream>
#include <cmath>
#include <vector>
#include "CooMatrix.h"
#include "utils.h"

double computeRMS(std::vector<double>& residual)
{
	double rmsValue = 0.0;
	for (int i = 0; i < residual.size(); ++i) {
		rmsValue += residual[i] * residual[i];
	}
	rmsValue = sqrt(rmsValue / ((double)residual.size()));
	return rmsValue;
}

std::vector<double> computeResidual(CooMatrix& A, std::vector<double>& b, std::vector<double>& x)
{
	std::vector<double> residual;
	residual.reserve(A.dim(0));
	int inz = 0;
	int Annz = A.nnz();
	for (int i = 0; i < A.dim(0); ++i) {
		double temp = 0;
		while (A.row(inz) == i) {
			temp += A.data(inz) * x[A.col(inz)];
			inz++;
			if (inz == Annz) {
				break;
			}
		}
		residual.push_back(b[i] - temp);
	}
	return residual;
}

void gauss_seidel_update(CooMatrix& A, std::vector<double>& b, std::vector<double>& x) {
	int inz = 0;
	int Annz = A.nnz();
	for (int i = 0; i < A.dim(0); i++) {
		double temp = 0;
		int idiag = 0;
		while (A.row(inz) == i) {
			if (A.col(inz) == A.row(inz)) {
				temp += b[i];
				idiag = inz;
				inz++;
				if (inz == Annz) {
					break;
				}				
			}
			else {
				temp -= A.data(inz) * x[A.col(inz)];
				inz++;
				if (inz == Annz) {
					break;
				}
			}
		}
		x[i] = (1 / A.data(idiag)) * temp;
	}
}

void gseidel(CooMatrix& A, std::vector<double>& b, std::vector<double>& x, int maxiters = 10000, double tol = 1e-6) {
	int iteration = 0;
	double error = tol + 1;
	size_t neq = x.size();
	std::vector<double> residual;
	residual.reserve(neq);
	while (iteration < maxiters && error > tol) {
		gauss_seidel_update(A, b, x);
		residual = computeResidual(A, b, x);
		error = computeRMS(residual);
		iteration = iteration + 1;
//		std::cout << "Iteration = " << iteration << "; |Residual| = " << error << std::endl;
	}
}

void mg_update(std::vector<CooMatrix>& A, std::vector<double>& b, std::vector<double>& x, std::vector<CooMatrix>& P, std::vector<CooMatrix>& R, size_t level, size_t nLevels) {
	if (level == nLevels - 2) {
		gseidel(A[level], b, x, 10);
	}
	else if (level == nLevels - 3) {
		gseidel(A[level], b, x, 5);
	}
	else if (level == nLevels - 4) {
		gseidel(A[level], b, x, 2);
	}
	else {
		gseidel(A[level], b, x, 1);
	}

	//	cout << "Pre-smoothing done for level " << level << endl;

	std::vector<double> resc = R[level] * computeResidual(A[level], b, x);

	//	cout << "Restriction 1 done for level " << level << endl;

	std::vector<double> eps = zeros(resc.size());

	if (level == nLevels - 2) {
		gseidel(A[level + 1], resc, eps);
	}
	else {
		mg_update(A, resc, eps, P, R, level + 1, nLevels);
	}

	updatevec(x, (std::vector<double>&)(P[level]*eps));

	//	cout << "Prolongation 1 done for level " << level << endl;

	std::vector<double> resc2 = R[level] * computeResidual(A[level], b, x);

	//	cout << "Restriction 2 done for level " << level << endl;

	if (level == nLevels - 2) {
		gseidel(A[level + 1], resc2, eps);
	}
	else {
		mg_update(A, resc2, eps, P, R, level + 1, nLevels);
	}

	updatevec(x, (std::vector<double>&)(P[level] * eps));

	//	cout << "Prolongation 2 done for level " << level << endl;

	if (level == nLevels - 2) {
		gseidel(A[level], b, x, 10);
	}
	else if (level == nLevels - 3) {
		gseidel(A[level], b, x, 5);
	}
	else if (level == nLevels - 4) {
		gseidel(A[level], b, x, 2);
	}
	else {
		gseidel(A[level], b, x, 1);
	}
	//	cout << "Post-smoothing done for level " << level << endl;
}

void multigrid(std::vector<CooMatrix>& A, std::vector<double>& b, std::vector<double>& x, std::vector<CooMatrix>& P, std::vector<CooMatrix>& R, int maxiters = 10000, double tol = 1e-6) {
	int iteration = 0;
	double error = tol + 1;
	size_t neq = x.size();
	size_t level0 = 0;
	std::vector<double> residual;
	size_t nLevels = P.size() + 1;
	std::cout << "Initiating multigrid solver for " << nLevels << " levels..." << std::endl;
	while (iteration < maxiters && error > tol) {
		mg_update(A, b, x, P, R, level0, nLevels);
		residual = computeResidual(A[0], b, x);
		error = computeRMS(residual);
		iteration = iteration + 1;
		std::cout << "Iteration = " << iteration << "; |Residual| = " << error << std::endl;
	}
}