/*
 * VectorSparse.cpp
 *
 *  Created on: 17 nov. 2014
 *      Author: jfellus
 */

#include "VectorSparse.h"

double vector_ps_float(const VectorSparse& v1, const VectorSparse& v2) {
	auto i1 = v1.entries.begin();
	auto i2 = v2.entries.begin();
	double s = 0;
	while(i1!=v1.entries.end() && i2!=v2.entries.end()) {
		if((*i1).i == (*i2).i) {
			s += (*i1).val * (*i2).val;
			i1++; i2++;
		}
		else if((*i1).i < (*i2).i) i1++; else i2++;
	}
	return s;
}

double vector_ps_float(const VectorSparse& v1, float* v2, size_t dim) {
	double s = 0;
	for(auto i1 = v1.entries.begin(); i1!=v1.entries.end(); i1++) {
//		DBG((*i1).val << "*" << v2[(*i1).i] << "at" << (*i1).i);
//		sleep(1);
		s += (*i1).val * v2[(*i1).i];
//		DBGV(s);
	}
	return s;
}

double vector_ps_float(float* v1, const VectorSparse& v2, size_t dim) {
	return vector_ps_float((const VectorSparse&)v2,(float*)v1, dim);
}
