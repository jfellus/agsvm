/*
 * matrix.cpp
 *
 *  Created on: 22 nov. 2013
 *      Author: jfellus
 */

#include "matrix.h"
#include "math.h"
#include "../vector_sparse/VectorSparse.h"

static MatrixFloat* cur_cov_mt_C;
static MatrixFloat* cur_cov_mt_X;
static __multithread__(compute_cov_mt) (int i) {
	MatrixFloat& X = *cur_cov_mt_X;
	MatrixFloat& C = *cur_cov_mt_C;

	for(int j=i; j<X.width; j++) {
		C[i+j*X.width] = 0;
		for(int k=0; k<X.height; k++) C[i+j*X.width] += X[i + k*X.width] * X[j + k*X.width];
		C[i+j*X.width] /= X.height;

		C[j+i*X.width] = C[i+j*X.width];
	}
}


void MatrixFloat::covariance_mt(MatrixFloat& X) {
	cur_cov_mt_C = this;
	cur_cov_mt_X = &X;
	compute_cov_mt(X.width);
}

void MatrixFloat::rand(float min, float max) {
	for(int i=0; i<width*height; i++)
		data[i] = frand(min,max);
}

MatrixFloat& MatrixFloat::operator=(const VectorSparse& v) {
	memset(data, 0, width*height*sizeof(float));
	for(auto i=v.entries.begin(); i!=v.entries.end(); i++) {
		data[(*i).i] = (*i).val;
	}
	return *this;
}

MatrixFloat& MatrixFloat::operator+=(const VectorSparse& v) {
	for(auto i = v.entries.begin(); i!=v.entries.end(); i++) {
		data[(*i).i] += (*i).val;
	}
	return *this;
}


MatrixFloat& MatrixFloat::operator-=(const VectorSparse& v) {
	for(auto i = v.entries.begin(); i!=v.entries.end(); i++) {
		data[(*i).i] -= (*i).val;
	}
	return *this;
}
