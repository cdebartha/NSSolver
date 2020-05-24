#pragma once

#include <iostream>
#include <vector>

struct Array3D
{
	unsigned int shape[3];
	std::vector<double> data;

	//
	Array3D() {};
	//Constructor
	Array3D(int i, int j, int k);
	//Copy Constructor
	Array3D(const Array3D& other);
	//Destructor
	~Array3D();

	//Assignment Operator
	Array3D& operator = (const Array3D& other);

	//Access info functions
	double& operator()(int i, int j, int k);
	double& operator()(int i, int j, int k, int istride, int jstride);
};
