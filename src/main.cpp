/*
 * main.cpp
 *
 *  Created on: 8 nov. 2013
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

using namespace std;



////////////
// PARAMS //
////////////

int N = get_config("N", 1);
int n = get_config("n", 100);
int D = get_config("D", 2);
string dataset = get_config_str("dataset", "data.jpg");

float LAMBDA = get_config("LAMBDA", 0.1);

int T_MAX = get_config("T_MAX", 1000);
int NB_NEIGHBORS = get_config("NB_NEIGHBORS", -1);
int NB_EDGES = 0;


//////////
// DATA //
//////////

Matrix X; // Sample
int* y; // Labels

float learningRate = 1;


int t = 0;

//////////////////////


#include "algebra.cpp"
int __node_last_id = 0;

//////////
// NODE //
//////////

class Node;
Node* node;
int last_sender, last_receiver;

class Node {
public:
	int id;

	Matrix X;
	int* y;
	Matrix w;
	float b;
	int n;


	Node() {this->id = __node_last_id++;}
	void init(Matrix& X, int* y, int first, int n) {
		this->y = new int[n];
		memcpy(this->y, &y[first], sizeof(int)*n);
		this->X.create_ref(&X[first*X.width],X.width,n);

		w.create(D, 1);
		b = 0;

		this->n = n;
	}

	~Node() {}


	// 2) GOSSIP

	void init_gossip() {
	}

	void iteration() {
		optimize();

		// Send
		last_receiver = gossip_choose_receiver();
		send(node[last_receiver]);
	}

	void optimize() {
		int i = rand()%n;
		float* sample = X.get_row(i);
		if( y[i] * vector_ps_float(w, sample, D) < 1) {
			for(int d=0; d<D; d++) w[d] = (1-learningRate*LAMBDA)*w[d] + learningRate*y[i]*sample[d];
		} else {
			for(int d=0; d<D; d++) w[d] = (1-learningRate*LAMBDA)*w[d];
		}
	}

	int send(Node& node) {


		node.receive(0);
		return 0;
	}

	void receive(int arg) {
		// TODO
	}


	void compute_estimate() {
	}


	//////////////////
	// CONNECTIVITY //
	//////////////////

	int* neighbors;
	int nbNeighbors;

	bool is_network_complete() {return neighbors==NULL;}
	int get_neighbors_count() {return nbNeighbors;}

	void connect(int neighbor) {
		node[neighbor].neighbors[node[neighbor].nbNeighbors++] = id;
		neighbors[nbNeighbors++] = neighbor;
		NB_EDGES += 2;
	}

	int gossip_choose_receiver() {
		if(N==1) return 0;
		if(is_network_complete()) {
			int r = rand()%(N-1);
			if(r>=id) r++;
			return r;
		} else {
			int r = rand()%get_neighbors_count();
			return neighbors[r];
		}
	}

};


#include "misc.cpp"



//////////
// DUMP //
//////////


void dump_classifier() {
	FILE* f = fopen(fmt("data/classifier_%u.txt", t), "w");
	fprintf(f, "%f %f %f", node[0].w[0], node[0].w[1], node[0].b);
	fclose(f);
}

void compute_errors() {
	if(last_sender!=-1) node[last_sender].compute_estimate();
	if(last_receiver!=-1) node[last_receiver].compute_estimate();
//
//	// Relative Error to Consensus
//	float FCnorm = FULL_COVARIANCE.n2();
//	float REC = 0;
//	for(int i=0; i<N; i++) REC += FULL_COVARIANCE.l2(node[i].C_rec);
//	REC /= (N*FCnorm);

//	fappend("data/E.txt", fmt("%f\n", REC));

	dump_classifier();
}



void init(const char* datafile) {
	DBG("INIT");
	DBGV(NBTHREADS);
	system("rm -rf data/*");
	system("rm -rf plots/*");

	if(str_has_extension(datafile, "jpg")) loadClassesFromJpg(datafile);
	else X.load(datafile);
//	if(LIMIT_NDATA!=-1 && X.height > LIMIT_NDATA) X.height = LIMIT_NDATA;
	n = X.height;
	D = X.width;

	node = new Node[N];
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
}


int main(int argc, char **argv) {
	try {
	DBG_START("Init ");
	if(argc<=1) {
		init(dataset.c_str());
	} else {
		DBG("with : " << argv[1]);
		init(argv[1]);
		//chdir(dirname(argv[1]));
	}

	if(system("mkdir -p data")) {}

	DBG_END();

	////////////////////////

	for(int i=0; i<N; i++) node[i].init_gossip();

	for(t=1; t<T_MAX; t++) {
		DBGV(t);

		// PEGASOS GLOBAL LEARNING RATE = 1/(lambda*t)
		learningRate = 1/(LAMBDA*t);

		last_sender = gossip_choose_sender();
		node[last_sender].iteration();
		compute_errors();
	}

	//////////////////////

	DBG("finished");
	}catch(const char* c) {DBG("ERROR : " << c );}
}
