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
		fprintf(f, "%f %f %i\n", X.get(i,0), X.get(i,1), y[i]==-1 ? 0 : 1);
	}
	fclose(f);
}



void dumpXY_sparse(const char* file) {
	FILE* f = fopen(file, "w");
	for(int i=0; i<n; i++) {
		fprintf(f, "%f %f %i\n", X.get(i,0), X.get(i,1), y[i]==-1 ? 0 : 1);
	}
	fclose(f);
}




void load_labels(Matrix& y, const char* filename) {
	std::ifstream f(filename);
	int n = 0;
	while(f.good()) { int c = f.get(); if(c=='\n') n++; }
	f.close();
	y.create(1,n);

	std::ifstream ff(filename);
	for(int i=0; ff.good(); i++) {
		ff >> y[i];
		int c;
		while(ff.good() && (c=ff.get())!='\n') {}
	}
}
