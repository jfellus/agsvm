/*
 * VectorSparse.h
 *
 *  Created on: 17 nov. 2014
 *      Author: jfellus
 */

#ifndef VECTORSPARSE_H_
#define VECTORSPARSE_H_

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <fstream>
#include "../common/utils.h"
#include "../common/matrix.h"




class VectorSparse;


double vector_ps_float(const VectorSparse& v1, const VectorSparse& v2) ;
inline double vector_ps_float(const VectorSparse& v1, const VectorSparse& v2, size_t dim) { return vector_ps_float(v1,v2); }

double vector_ps_float(const VectorSparse& v1, float* v2, size_t dim);

double vector_ps_float(float* v1, const VectorSparse& v2, size_t dim);



typedef struct {
	uint i;
	double val;
} entry;

class VectorSparse {
public:
	std::list<entry> entries;
	uint width = 0;

	double* dense = NULL;
	uint dense_d = 0;
	bool bUpdateDense = true;
public:
	VectorSparse() {}
	virtual ~VectorSparse() {}
	VectorSparse(const VectorSparse& v) {(*this)=v;	}
	VectorSparse(float* v, size_t dim) {set(v,dim);}
	VectorSparse(const MatrixFloat& m) {set(m.data,m.width*m.height);}


	double get(size_t i) {
		for(auto j = entries.begin(); j!=entries.end(); j++) {
			if((*j).i == i) return (*j).val;
			if((*j).i > i) return 0;
		}
		return 0;
	}

	VectorSparse& operator=(const VectorSparse& v) {
		entries.clear();
		entries = v.entries;
		width = v.width;
		bUpdateDense = true;
		return *this;
	}

	VectorSparse& set(float* v, size_t dim) {
		bUpdateDense = true;
		entries.clear();
		for(uint i=0; i<dim; i++) {
			if(v[i]==0) continue;
			entry e = {i,v[i]};
			entries.push_back(e);
		}
		width = dim;
		return *this;
	}

	VectorSparse& set(uint i, double val) {
		bUpdateDense=true;
		entry e = {i,val};
		for(auto j = entries.begin(); j!=entries.end(); j++) {
			if((*j).i == i) {(*j).val = val; return *this;}
			else if((*j).i > i) {entries.insert(j, e); return *this;}
		}
		entries.push_back(e);
		width = e.i+1;
		return *this;
	}

	VectorSparse& operator*=(const VectorSparse& v) {
		auto i1 = entries.begin();
		auto i2 = v.entries.begin();
		while(i1!=entries.end() && i2!=v.entries.end()) {
			if((*i1).i == (*i2).i) {
				(*i1).val *= (*i2).val;
				i1++; i2++;
			}
			else if((*i1).i < (*i2).i) { entries.erase(i1++); }
			else i2++;
		}
		while(i1!=entries.end()) { entries.erase(i1++);}
		bUpdateDense = true;
		return *this;
	}

	VectorSparse& operator+=(const VectorSparse& v) {
		auto i1 = entries.begin();
		auto i2 = v.entries.begin();
		while(i1!=entries.end() && i2!=v.entries.end()) {
		//	DBG((*i1).i << " <-> " << (*i2).i);
			if((*i1).i == (*i2).i) {
				(*i1).val += (*i2).val;
				i1++; i2++;
			}
			else if((*i1).i > (*i2).i) {
				entries.insert(i1, (*i2));
				i2++;
			} else i1++;
		}
	//	DBG("---------------\n\n");
		while(i2!=v.entries.end()) {	entries.push_back(*i2); width = width>((*i2).i+1) ? width : ((*i2).i+1); i2++; }
		bUpdateDense = true;
		return *this;
	}

	VectorSparse& operator-=(const VectorSparse& v) {
		auto i1 = entries.begin();
		auto i2 = v.entries.begin();
		while(i1!=entries.end() && i2!=v.entries.end()) {
			if((*i1).i == (*i2).i) {
				(*i1).val -= (*i2).val;
				i1++; i2++;
			}
			else if((*i1).i > (*i2).i) {
				entry e = *i2; e.val = -e.val;
				entries.insert(i1, e);
				i2++;
			} else i1++;
		}
		while(i2!=v.entries.end()) {
			entry e = *i2; e.val = -e.val;
			entries.push_back(e); width = width>((*i2).i+1) ? width : ((*i2).i+1); i2++;
		}
		bUpdateDense = true;
		return *this;
	}

	VectorSparse& operator*=(float f) {
		for(auto i = entries.begin(); i!=entries.end(); i++) {
			(*i).val *= f;
		}
		bUpdateDense = true;
		return *this;
	}

	VectorSparse& operator/=(float f) {
		for(auto i = entries.begin(); i!=entries.end(); i++) {
			(*i).val /= f;
		}
		bUpdateDense = true;
		return *this;
	}


	operator double*() {
		compute_dense();
		return dense;
	}

	void set_width(uint d) {
		if(width < d) {	set(d-1, 0); bUpdateDense = true;}
	}

	void compute_dense() {
		if(dense!=NULL && dense_d != width) {delete dense; dense=0; dense_d=0;}
		if(dense==NULL) {dense_d = width; dense = new double[dense_d];}
		if(bUpdateDense) {
			memset(dense, 0, sizeof(double)*dense_d);
			for(auto i = entries.begin(); i!=entries.end(); i++) {
				dense[(*i).i] = (*i).val;
			}
			bUpdateDense = false;
		}
	}

	void dump() {
		std::cout << "[ ";
		auto i1 = entries.begin();
		for(uint i=0; i<width; i++) {
			if((*i1).i == i) { std::cout << (*i1).val << " "; i1++;}
			else std::cout << "0 ";
		}
		std::cout << "]\n";
	}

	void dump_internal() {
		std::cout << "[ ";

		for(auto i1 = entries.begin(); i1!=entries.end(); i1++) {
			std::cout << (*i1).i << ":" << (*i1).val << " ";
		}
		std::cout << "]\n";
	}

	double norm_l2() const {
		return vector_ps_float(*this, *this);
	}

	void clear() {
		entries.clear(); width = 0;
		bUpdateDense = true;
	}

	void read(std::ifstream& f) {
		int c = 0;
		entry e; e.i = -1;
		while(f.good() && c!='\n') {
			f >> e.i;
			f.get();
			f >> e.val;
			entries.push_back(e);
			c = f.get();
		}
		width = e.i+1;
	}
};


inline float l2(const VectorSparse& v) {return v.norm_l2();}


inline VectorSparse operator *(float f, const VectorSparse& v) {
	VectorSparse v2 = v;
	v2 *= f;
	return v2;
}

inline VectorSparse operator *(const VectorSparse& v, float f) {
	VectorSparse v2 = v;
	v2 *= f;
	return v2;
}



inline VectorSparse operator +(const VectorSparse& v1, const VectorSparse& v2) {
	VectorSparse v3 = v1;
	v3 += v2;
	return v3;
}

inline VectorSparse operator -(const VectorSparse& v1, const VectorSparse& v2) {
	VectorSparse v3 = v1;
	v3 -= v2;
	return v3;
}


inline MatrixFloat operator -(const MatrixFloat& v1, const VectorSparse& v2) {
	MatrixFloat v = v1;
	for(auto i = v2.entries.begin(); i != v2.entries.end(); i++) {
		v.data[(*i).i] -= (*i).val;
	}
	return v;
}

inline MatrixFloat operator +(const MatrixFloat& v1, const VectorSparse& v2) {
	MatrixFloat v = v1;
	for(auto i = v2.entries.begin(); i != v2.entries.end(); i++) {
		v.data[(*i).i] += (*i).val;
	}
	return v;
}




#endif /* VECTORSPARSE_H_ */
