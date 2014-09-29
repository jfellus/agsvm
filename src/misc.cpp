/*
 * misc.cpp
 *
 *  Created on: 13 nov. 2013
 *      Author: jfellus
 */
#include "retin/toolbox/imcodec/ImageCodec.h"

using namespace retin;


///////////
// CONNECTIVITY




void dump_degrees() {
	for(int i=0; i<N; i++) fappend("data/stats/degrees.txt", fmt("%u\n", node[i].nbNeighbors));
}


void create_network() {
	node = new Node[N];
	if(NB_NEIGHBORS==-1) {
		DBG("Fully connected network");
		for(int i=0; i<N; i++) {	node[i].nbNeighbors = N-1; node[i].neighbors = NULL;	}
		NB_EDGES = N*(N-1);
	} else {
		DBG("BA-model network");
		DBG("  max degree = " << NB_NEIGHBORS);
		// BA-model (preferential attachment)
		for(int i=0; i<N; i++) { node[i].neighbors = new int[NB_NEIGHBORS]; node[i].nbNeighbors = 0; }
		node[0].connect(1);

		for(int i=2; i<N; i++) {
			for(int j=0; j<i-1; j++) {
				if(node[j].nbNeighbors >= NB_NEIGHBORS) continue;
				if(node[j].nbNeighbors == 0) continue;
				int e = rand()%NB_EDGES;
				if(e<=node[j].nbNeighbors) node[i].connect(j);
			}
			node[i].connect(i-1);
		}

		dump_degrees();
	}
}


bool str_has_extension(const char* s, const char* e) {
	if(strlen(s)-1-strlen(e) < 0 || s[strlen(s)-1-strlen(e)]!='.') return false;
	return !strcmp(&s[strlen(s)-strlen(e)], e);
}



void deinit() {
	delete[] node;
}



void dumpXY(const char* file) {
	FILE* f = fopen(file, "w");
	for(int i=0; i<n; i++) {
		fprintf(f, "%f %f %i\n", X[i*X.width], X[i*X.width +1], y[i]==-1 ? 0 : 1);
	}
	fclose(f);
}


void loadClassesFromJpg(const char* datafile) {
	DBG("Load from JPEG : " << datafile);
	DBGV(n);
	string s = datafile;
	unsigned char* pxl;
	unsigned char* palette;
	size_t w,h,c;

	X.create(2, n);
	y.create(1, n);

	loadImage(pxl, palette, w,h,c, s);


	for(int i=0; i<n; i++) {
		size_t _x = rand()%w;
		size_t _y = rand()%h;
		X.get_row(i)[0] = _x;
		X.get_row(i)[1] = _y;
		unsigned char* p = &pxl[((h-_y-1)*w + _x)*3];
		y[i] =  (p[0] > 128) ? 1 : -1;
	}

	dumpXY("data.txt");
}

void loadTestClassesFromJpg(const char* datafile) {
	DBG("Load test data from JPEG : " << datafile);
	DBGV(n);
	string s = datafile;
	unsigned char* pxl;
	unsigned char* palette;
	size_t w,h,c;

	X_test.create(2, n);
	y_test.create(1, n);

	loadImage(pxl, palette, w,h,c, s);


	for(int i=0; i<n; i++) {
		size_t _x = rand()%w;
		size_t _y = rand()%h;
		X.get_row(i)[0] = _x;
		X.get_row(i)[1] = _y;
		unsigned char* p = &pxl[((h-_y-1)*w + _x)*3];
		y[i] =  (p[0] > 128) ? 1 : -1;
	}

	dumpXY("data_test.txt");
}

