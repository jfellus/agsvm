/*
 * test.cpp
 *
 *  Created on: 18 nov. 2014
 *      Author: jfellus
 */


#include "VectorSparse.h"
#include "../common/utils.h"

int main(int argc, char **argv) {
	float* v1 = new float[10];
		for(uint i=0; i<10; i++) v1[i] = rand()%2 ? 3 : 0;
	float* v3 = new float[10];
		for(uint i=0; i<10; i++) v3[i] = rand()%2 ? 4 : 0;

	VectorSparse v2(v1, 10);
	v2.dump();

	VectorSparse v4(v3, 10);
	v4.dump();

	VectorSparse v5 = v2;
	v5 *= v4;

	v5.dump();

	VectorSparse v6;
	v6.set(2, 0.5);
	v6.set(4, 0.5);
	v6.set(6, 0.5);
	v6.dump();

	v5 += v6;

	v5.dump();

	v5 *= 0.5;
	v5.dump();

	v5.set(2,1);
	v5.set(3,1);
	v5.set(4,1);
	v5.set(6,2);
	v5.dump();

	DBG("Norm L2 : " << v5.norm_l2());

}
