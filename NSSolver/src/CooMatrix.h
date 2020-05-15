#pragma once

#include <vector>
#include <iostream>

class CooMatrix
{
private:
	std::vector<int> rowind_;
	std::vector<int> colind_;
	std::vector<double> val_;
	unsigned int nrows_;
	unsigned int ncols_;
	unsigned int nnz_;

public:
	//Constructors
	CooMatrix(unsigned int nrows, unsigned int ncols, unsigned int nnz, std::vector<int>& rowind, std::vector<int>& colind, std::vector<double>& val);
	//Copy Constructor
	CooMatrix(const CooMatrix& other);
	//Destructor
	~CooMatrix();

	//Assignment Operator
	CooMatrix& operator = (const CooMatrix& other);

	//Access info functions
	double& data(unsigned int i) { return val_[i]; };
	int& row(unsigned int i) { return rowind_[i]; };
	int& col(unsigned int i) { return colind_[i]; };
	int dim(int i);
	int nnz() { return nnz_; };

	//Vector Multiplication
	std::vector<double> operator * (const std::vector<double>& other);

	//print
	friend std::ostream& operator << (std::ostream& os, const CooMatrix& mat);
};

