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




//////////
// DUMP //
//////////

string fE, fEstddev;
std::ofstream ffE, ffEstddev;

string fmt_padd_float(float f) {
	return fmt("%012.3f",f);
}

void dump_classifier() {
	shell(TOSTRING("mkdir -p data" << PREFIX << "w/" << ALGO << "/N" << N << "/M" << NB_MESSAGES << "/l" << LEARNING_RATE << "/"));
	node[0].w.write(TOSTRING("data" << PREFIX << "w/" << ALGO << "/N" << N << "/M" << NB_MESSAGES << "/l" << LEARNING_RATE << "/" << fmt_padd_float((((float)nbgradients_evaluated/::N))) << ".fvec").c_str());
}



bool file_exists(const char* s) {
	return access(s, F_OK) != -1;
}

string newfilename(const char* s) {
	size_t ii = 0;
	string fE = fmt(s, ii++);
	while(file_exists(fE.c_str())) fE = fmt(s, ii++);
	return fE;
}




void init() {
	DBG("INIT");
	DBGV(NBTHREADS);

	shell("mkdir -p data/w");
	shell(TOSTRING("mkdir -p data" << PREFIX));
	shell(TOSTRING("mkdir -p data" << PREFIX << "w"));

	//	system("rm -rf data/*");
	shell("rm -rf plots/*");
	if(ALGO=="STAG" || ALGO=="STAGR") {
		if(N==1) {
		fE = newfilename(fmt("data%sE_%s_N%u_%u_%f_%%u.txt", PREFIX.c_str(), ALGO.c_str(),N, STAG_BUFFER_SIZE, LEARNING_RATE));
		fEstddev = newfilename(fmt("data%sDEV_%s_%u_%f_%%u.txt", PREFIX.c_str(), ALGO.c_str(), STAG_BUFFER_SIZE, LEARNING_RATE));
		} else {
			shell(fmt("mkdir -p data%sN%u", PREFIX.c_str(), N));
			shell(fmt("mkdir -p data%sN%u/stddev", PREFIX.c_str(), N));
			fE = newfilename(fmt("data%sN%u/E_%s_N%u_%u_%f_%%u.txt", PREFIX.c_str(), N, ALGO.c_str(), N, STAG_BUFFER_SIZE, LEARNING_RATE));
			fEstddev = newfilename(fmt("data%sN%u/stddev/DEV_%s_N%u_%u_%f_%%u.txt", PREFIX.c_str(), N, ALGO.c_str(), N, STAG_BUFFER_SIZE, LEARNING_RATE));
		}
		ffE.open(fE, ios_base::app);
	}
	else {
		fE = newfilename(fmt("data%sE_%s_N%u_%f_M%f_%%u.txt", PREFIX.c_str(), ALGO.c_str(), N, LEARNING_RATE, NB_MESSAGES));
//		fEstddev = newfilename(fmt("data/DEV_%s_N%u_%f_%%u.txt", ALGO.c_str(), N, LEARNING_RATE));
		ffE.open(fE, ios_base::app);
//		ffEstddev.open(fEstddev, ios_base::app);
	}


	if(str_has_extension(dataset.c_str(), "jpg")) {
		loadClassesFromJpg(dataset.c_str());
		loadTestClassesFromJpg(dataset_test.c_str());
	}
	else {
		DBGV(dataset);
		DBGV(labels);
		DBGV(dataset_test);
		DBGV(labels_test);

		if(ADD_BIAS) {
				Matrix X_nobias;
				X_nobias.load(dataset.c_str());
				X.create(X_nobias.width+1,X_nobias.height);
				for(int i=0; i<X.height; i++) {
					memcpy(X.get_row(i),X_nobias.get_row(i),X_nobias.width*sizeof(float));
					X[i*X.width + X_nobias.width] = 1;
				}
		}
		else
			X.load(dataset.c_str());
		y.load(labels.c_str());

		if(ADD_BIAS) {
			Matrix X_test_nobias;
			X_test_nobias.load(dataset_test.c_str());
			X_test.create(X_test_nobias.width+1,X_test_nobias.height);
			for(int i=0; i<X_test.height; i++) {
				memcpy(X_test.get_row(i),X_test_nobias.get_row(i),X_test_nobias.width*sizeof(float));
				X_test[i*X_test.width + X_test_nobias.width] = 1;
			}
		}
		else X_test.load(dataset_test.c_str());
		y_test.load(labels_test.c_str());

		if(CATEGORY != -1) {
			for(int i=0; i<y.height; i++) y[i] = y[i]==CATEGORY ? 1 : -1;
			for(int i=0; i<y_test.height; i++) y_test[i] = y_test[i]==CATEGORY ? 1 : -1;
		}
	}
//	if(LIMIT_NDATA!=-1 && X.height > LIMIT_NDATA) X.height = LIMIT_NDATA;
	n = X.height;
	D = X.width;

	create_network();

	int ndo = 0;
	for(int i=0; i<N-1; i++) {
		node[i].init(X, y, ndo, n/N);
		ndo += n/N;
	}
	node[N-1].init(X, y, ndo, n-ndo);

	DBGV(N);
	DBGV(D);
	DBGV(n);
	DBGV(X_test.width);
	DBGV(X_test.height);
	sleep(1);
}
