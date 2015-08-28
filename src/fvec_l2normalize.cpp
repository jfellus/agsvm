/*
 * fvec_dump.cpp
 *
 *  Created on: 27 févr. 2015
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


std::string file;
std::string format = "txt";
int N = 1;



int main(int argc, char **argv) {
	if(argc<3) {
		DBG("Usage : " << argv[0] << " <in.fvec> <out.fvec>");
	}
	file = argv[1];

	Matrix X;
	X.load(file.c_str());
	for(uint i=0; i<X.height; i++) {
		vector_sdiv_float(&X[i*X.width], vector_n2_float(&X[i*X.width], X.width), X.width);
	}
	X.write(argv[2]);
}
