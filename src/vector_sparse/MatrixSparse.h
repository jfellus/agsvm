/*
 * MatrixSparse.h
 *
 *  Created on: 17 nov. 2014
 *      Author: jfellus
 */

#ifndef MATRIXSPARSE_H_
#define MATRIXSPARSE_H_

#include "../common/matrix.h"
#include "VectorSparse.h"

class MatrixSparse {
public:
	VectorSparse* vectors;
	size_t width,height;
public:
	MatrixSparse() { width = height = 0; vectors = NULL;}
	MatrixSparse(Matrix& m) {
		width = m.width;
		height = m.height;
		vectors = new VectorSparse[height];
		for(size_t i=0; i<height; i++) vectors[i].set(m.get_row(i), width);
	}

	MatrixSparse(size_t w, size_t h) : width(w), height(h) {
		vectors = new VectorSparse[h];
	}

	virtual ~MatrixSparse() {
		delete vectors;
	}

	inline VectorSparse& get_row(size_t i) { return vectors[i]; }

	inline double get(size_t i, size_t j) {return vectors[i].get(j);}

	inline VectorSparse& operator[](int i) {return get_row(i);}

	inline double set(size_t i, size_t j, double x) {
		get_row(i).set(j, x);
		if(width<get_row(i).width) width = get_row(i).width;
		return x;
	}


	void create(size_t w, size_t h) {
		width = w;
		height = h;
		vectors = new VectorSparse[h];
	}

	void create_ref(const MatrixSparse& m, size_t first_row, size_t nrows) {
		width = m.width;
		height = nrows;
		vectors = &m.vectors[first_row];
	}


	void load(const char* filename) {
		DBG("Load " << filename);
		std::ifstream f(filename);
		height = 0;
		while(f.good()) { int c = f.get(); if(c=='\n') height++; }
		f.close();

		vectors = new VectorSparse[height];

		width = 0;
		std::ifstream ff(filename);
		for(int i=0; ff.good(); i++) {
			int label; ff >> label;
			vectors[i].read(ff);
			if(width<vectors[i].width) width =vectors[i].width;
		}
		ff.close();
		DBG("ok");
	}

	operator bool() {return vectors!=NULL;}

	void clear() {
		for(uint i=0; i<height; i++) get_row(i).clear();
	}

};

#endif /* MATRIXSPARSE_H_ */
