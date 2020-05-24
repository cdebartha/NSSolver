#include "Array3D.h"

Array3D::Array3D(int i, int j, int k)
{
	this->shape[0] = i;
	this->shape[1] = j;
	this->shape[2] = k;
	this->data.resize(i * j * k);
}

Array3D::Array3D(const Array3D& other) 
{
	this->shape[0] = other.shape[0];
	this->shape[1] = other.shape[1];
	this->shape[2] = other.shape[2];
	this->data = other.data;
}

Array3D::~Array3D()
{

}

Array3D& Array3D::operator = (const Array3D& other) 
{
	this->shape[0] = other.shape[0];
	this->shape[1] = other.shape[1];
	this->shape[2] = other.shape[2];
	this->data = other.data;
	return *this;
}

double& Array3D::operator()(int i, int j, int k)
{
	return this->data[i + this->shape[0] * j + this->shape[0] * this->shape[1] * k];
}

double& Array3D::operator()(int i, int j, int k, int istride, int jstride)
{
	return this->data[i + istride * j + istride * jstride * k];
}