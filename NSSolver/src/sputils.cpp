#include <iostream>
#include <vector>
#include "CooMatrix.h"

void print1D(const std::vector<double>& vec)
{
	std::cout << "[ " << std::flush;
	for (auto& element : vec) {
		std::cout << element << ", " <<std::flush;
	}
	std::cout << "]" << std::endl;
}

void printmat(CooMatrix& mat) {
	double zero = 0.0;
	int nz = 0;
	for (int i = 0; i < mat.dim(0); ++i) {
		int j = 0;
		std::cout << "[ " << std::flush;
		while (mat.row(nz) == i && j < mat.dim(1)) {
			int Aj = mat.col(nz);
			if (Aj == j) {
				std::cout << mat.data(nz) << ", " << std::flush;
				++j;
				++nz;
				if (nz == mat.nnz()) {
					break;
				}
			}
			else {
				std::cout << zero << ", " << std::flush;
				++j;
			}
		}
		while (j < mat.dim(1)) {
			std::cout << zero << ", " << std::flush;
			++j;
		}
		std::cout << "]," << std::endl;
	}
}

std::vector<double> linspace(double start, double end, int steps)
{
	std::vector<double> linearspace;
	linearspace.reserve(steps);
	double delta = (end - start) / ((double)steps - 1.0);
	for (int i = 0; i < steps - 1; ++i) {
		linearspace.push_back(start + delta * (double)i);
	}
	linearspace.push_back(end);
	return linearspace;
}

std::vector<double> zeros(int N) {
	std::vector<double> zeroArray;
	zeroArray.reserve(N);
	for (int i = 0; i < N; ++i) {
		zeroArray.push_back(0.0);
	}
	return zeroArray;
}

std::vector<double> addvec(std::vector<double>& a, std::vector<double>& b) 
{
	if (a.size() != b.size()) {
		throw "For addition, both arrays must have same dimensions.";
	}
	std::vector<double> sumArray;
	sumArray.reserve(a.size());
	for (int i = 0; i < a.size(); ++i) {
		sumArray.push_back(a[i] + b[i]);
	}
	return sumArray;
}

void updatevec(std::vector<double>& x, const std::vector<double>& b)
{
	if (x.size() != b.size()) {
		throw "For update, both arrays must have same dimensions.";
	}

	for (int i = 0; i < x.size(); ++i) {
		x[i] += b[i];
	}
}

CooMatrix eye(int N)
{
	std::vector<int> rowind;
	std::vector<int> colind;
	std::vector<double> val;
	for (int i = 0; i < N; i++) {
		rowind.push_back(i);
		colind.push_back(i);
		val.push_back(1.0);
	}
	return CooMatrix(N, N, N, rowind, colind, val);
}

CooMatrix kron(CooMatrix& A, CooMatrix& B) {
	unsigned int M = A.dim(0) * B.dim(0);
	unsigned int N = A.dim(1) * B.dim(1);
	unsigned int Annz = A.nnz();
	unsigned int Bnnz = B.nnz();
	unsigned int nnz = Annz * Bnnz;
	std::vector<int> rowind;
	std::vector<int> colind;
	std::vector<double> val;
	rowind.reserve(nnz);
	colind.reserve(nnz);
	val.reserve(nnz);
	int inz = 0;
	int jnz = 0;
	int counteri = 0;
	int counterj = 0;
	for (int i = 0; i < A.dim(0); ++i) {
		counterj = 0;
		for (int j = 0; j < B.dim(0); ++j) {
			inz = counteri;
			while (A.row(inz) == i) {
				jnz = counterj;
				while (B.row(jnz) == j) {
					rowind.push_back(A.row(inz) * B.dim(0) + B.row(jnz));
					colind.push_back(A.col(inz) * B.dim(1) + B.col(jnz));
					val.push_back(A.data(inz) * B.data(jnz));
					++jnz;
					if (jnz == Bnnz) {
						break;
					}
				}
				inz++;
				if (inz == Annz) {
					break;
				}
			}
			counterj = jnz;
		}
		counteri = inz;
	}
	return CooMatrix(M, N, nnz, rowind, colind, val);
}

CooMatrix addmat(CooMatrix& A, CooMatrix& B) {
	if (A.dim(0) != B.dim(0) || A.dim(1) != B.dim(1)) {
		throw "For addition, both arrays must have same dimensions.";
	}
	unsigned int M = A.dim(0);
	unsigned int N = A.dim(1);
	unsigned int Annz = A.nnz();
	unsigned int Bnnz = B.nnz();
	std::vector<int> rowind;
	std::vector<int> colind;
	std::vector<double> val;
	rowind.reserve((size_t)Annz + (size_t)Bnnz);
	colind.reserve((size_t)Annz + (size_t)Bnnz);
	val.reserve((size_t)Annz + (size_t)Bnnz);
	int nz = 0;
	int inz = 0;
	int jnz = 0;
	for (int i = 0; i < A.dim(0); ++i) {
		while (A.row(inz) == i && B.row(jnz) == i) {
			int Aj = A.col(inz);
			int Bj = B.col(jnz);

			if (Aj == Bj) {
				double result = A.data(inz) + B.data(jnz);
				if (result != 0) {
					rowind.push_back(i);
					colind.push_back(Aj);
					val.push_back(result);
					++nz;
				}
				++inz;
				++jnz;
				if (inz == Annz) {
					break;
				}
				if (jnz == Bnnz) {
					break;
				}
			}
			else if (Aj < Bj) {
				double result = A.data(inz);
				if (result != 0) {
					rowind.push_back(i);
					colind.push_back(Aj);
					val.push_back(result);
					++nz;
				}
				++inz;
				if (inz == Annz) {
					break;
				}
			}
			else {
				double result = B.data(jnz);
				if (result != 0) {
					rowind.push_back(i);
					colind.push_back(Bj);
					val.push_back(result);
					++nz;
				}
				++jnz;
				if (jnz == Bnnz) {
					break;
				}
			}
		}
		//tail
		if (inz != Annz) {
			while (A.row(inz) == i) {
				double result = A.data(inz);
				if (result != 0) {
					rowind.push_back(i);
					colind.push_back(A.col(inz));
					val.push_back(result);
					++nz;
				}
				++inz;
				if (inz == Annz) {
					break;
				}
			}
		}

		if (jnz != Bnnz) {
			while (B.row(jnz) == i) {
				double result = B.data(jnz);
				if (result != 0) {
					rowind.push_back(i);
					colind.push_back(B.col(jnz));
					val.push_back(result);
					++nz;
				}
				++jnz;
				if (jnz == Bnnz) {
					break;
				}
			}
		}
	}
	return CooMatrix(M, N, (unsigned int)val.size(), rowind, colind, val);
}

CooMatrix kronsum(CooMatrix& A, CooMatrix& B) {
	int M = A.dim(0);
	int N = B.dim(0);
	CooMatrix Im = eye(M);
	CooMatrix In = eye(N);
	CooMatrix kronIA = kron(In, A);
	CooMatrix kronBI = kron(B, Im);
	return addmat(kronIA, kronBI);
}