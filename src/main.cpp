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

int T_MAX = get_config("T_MAX", 1000);
int NB_NEIGHBORS = get_config("NB_NEIGHBORS", -1);
int NB_EDGES = 0;


//////////
// DATA //
//////////

Matrix X;
bool* y;

Matrix w;
float b;

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

	Node() {this->id = __node_last_id++;}
	void init(Matrix& X, int first, int n) {
	}

	~Node() {}


	// 2) GOSSIP

	void init_gossip() {
	}

//	int send(Node& node) {
//		node.receive(U,L,mu,w);
//		return 0;
//	}
//
//	void receive(Matrix& U2, Matrix& L2, Matrix& mu2, float w2) {
//		Matrix cr(D,D);
//		cr.reconstruct(U2,L2);
//
//		w += w2;	mu += mu2;		C += cr;
//	}
//

	void compute_estimate() {
	}


	// CONNECTIVITY
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

/////////////////////

void gossip() {
	last_sender = gossip_choose_sender();
	last_receiver = node[last_sender].gossip_choose_receiver();

//	last_nb_projs_sent = node[last_sender].send(node[last_receiver]);

//	last_msg_size = (last_nb_projs_sent+1)*(D+1);
//	total_msg_size += last_msg_size;

//	fappend("data/stats/total_msg_size.txt",fmt("%d\n",total_msg_size));
//	fappend("data/stats/nb_projs_sent.txt",fmt("%d\n",last_nb_projs_sent));
//	fappend("data/stats/total_msg_size_N.txt",fmt("%f %f\n",(float)t/N, (float)total_msg_size/N));
}


//////////
// DUMP //
//////////

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
		node[i].init(X, ndo, n/N);
		ndo += n/N;
	}
	node[N-1].init(X,ndo, n-ndo);

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

//	if(system("mkdir -p data")) {}

	DBG_END();

	////////////////////////

	for(int i=0; i<N; i++) node[i].init_gossip();

	for(t=0; t<T_MAX; t++) {
		DBGV(t);
		gossip();
		compute_errors();
	}

	//////////////////////

	DBG("finished");
	}catch(const char* c) {DBG("ERROR : " << c );}
}
