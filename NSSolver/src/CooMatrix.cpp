#include "CooMatrix.h"

CooMatrix::CooMatrix(unsigned int nrows, unsigned int ncols, unsigned int nnz, std::vector<int>& rowind, std::vector<int>& colind, std::vector<double>& val):
	nrows_(nrows), ncols_(ncols), nnz_(nnz), rowind_(rowind), colind_(colind), val_(val)
{
//	std::cout << "CooMatrix created." << std::endl;
}

CooMatrix::CooMatrix(const CooMatrix& other):
	nrows_(other.nrows_), ncols_(other.ncols_), nnz_(other.nnz_), rowind_(other.rowind_), colind_(other.colind_), val_(other.val_)
{
//	std::cout << "CooMatrix copied." << std::endl;
}

CooMatrix::~CooMatrix()
{
//	std::cout << "CooMatrix destructed." << std::endl;
}

CooMatrix& CooMatrix::operator = (const CooMatrix& other)
{
	this->nrows_ = other.nrows_;
	this->ncols_ = other.ncols_;
	this->nnz_ = other.nnz_;
	this->rowind_ = other.rowind_;
	this->colind_ = other.colind_;
	this->val_ = other.val_;
//	std::cout << "CooMatrix assigned." << std::endl;
	return *this;
}

int CooMatrix::dim(int i)
{
	if (i == 0) {
		return this->nrows_;
	}
	else if (i == 1) {
		return this->ncols_;
	}
	else {
		throw "Invalid dimension!";
	}
}

std::vector<double> CooMatrix::operator * (const std::vector<double>& other)
{
	std::vector<double> product(this->nrows_);
//	product.reserve(this->nrows_);
	unsigned int inz = 0;
	for (unsigned int i = 0; i < this->nrows_; ++i) {
		double temp = 0.0;
		while (this->rowind_[inz] == i) {
			temp += this->val_[inz] * other[this->colind_[inz]];
			++inz;
			if (inz == this->nnz_) {
				break;
			}
		}
		product[i] = temp;
	}
	return product;
}


std::ostream& operator << (std::ostream& os, const CooMatrix& mat)
{
	for (unsigned int inz = 0; inz < mat.nnz_; ++inz) {
		os << mat.rowind_[inz] << '\t' << mat.colind_[inz] << '\t' << mat.val_[inz] << std::endl;
	}
	return os;
}