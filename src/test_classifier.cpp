/*
 * test_classifier.cpp
 *
 *  Created on: 27 août 2015
 *      Author: jfellus
 */


#define MONOTHREAD

#include "common/math.h"
#include "common/multithread.h"
#include "common/utils.h"
#include "common/plot.h"
#include "common/gossip.h"
#include "retin/toolbox/imcodec/ImageCodec.h"
#include <vector>
#include <string>
#include <libgen.h>
#include <unistd.h>

using namespace std;
int N = 1;

Matrix X;
Matrix w;
Matrix res;

int main(int argc, char **argv) {
	if(argc < 4) {
		DBG("USAGE : " << argv[0] << " <classifier.fvec> <data.fvec> <outfile.txt>");
		exit(1);
	}
	string classifier = argv[1];
	string dataset = argv[2];
	string out = argv[3];


	// Init
	Matrix X_nobias;
	X_nobias.load(dataset.c_str());
	X.create(X_nobias.width+1,X_nobias.height);
	for(int i=0; i<X.height; i++) {
		memcpy(X.get_row(i),X_nobias.get_row(i),X_nobias.width*sizeof(float));
		X[i*X.width + X_nobias.width] = 1;
	}

	w.load(classifier.c_str());


	// Test
	res.create(1, X.height);
	for(int i=0; i<X.height; i++) {
		res[i] = vector_ps_float(w, X.get_row(i), X.width);
	}

	// Write
	FILE* f = fopen(out.c_str(), "w");
	for(int i=0; i<X.height; i++) {
		fprintf(f, "%f\n", res[i]);
	}
	fclose(f);
}

