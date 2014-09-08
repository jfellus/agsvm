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
	if(strlen(s)-2-strlen(e) < 0 || s[strlen(s)-2-strlen(e)]!='.') return false;
	return !strcmp(&s[strlen(s)-1-strlen(e)], e);
}



void deinit() {
	delete[] node;
}





void loadClassesFromJpg(const char* datafile) {
	string s = datafile;
	unsigned char* pxl;
	unsigned char* palette;
	size_t w,h,c;
	loadImage(pxl, palette, w,h,c, s);
	X.create(2, n);
	y = new bool[n];
	for(int i=0; i<n; i++) {
		size_t _x = rand()%w;
		size_t _y = rand()%h;
		X.get_row(i)[0] = _x;
		X.get_row(i)[1] = _y;
		unsigned char* p = &pxl[(_y*w + _x)*3];
		y[i] =  (p[0] > 128);
	}
}
