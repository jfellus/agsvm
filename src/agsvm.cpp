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
#include <unistd.h>

using namespace std;



////////////
// PARAMS //
////////////

int N = get_config("N", 1);
int n = get_config("n", 100);
int D = get_config("D", 2);
string dataset = get_config_str("dataset", "data.jpg");
string labels = get_config_str("labels", "");
string dataset_test = get_config_str("dataset_test", "data.jpg");
string labels_test = get_config_str("labels_test", "");

int CATEGORY = get_config("CATEGORY", 1);

float LAMBDA = get_config("LAMBDA", 0.1);

int T_MAX = get_config("T_MAX", 2000000*N);
bool ADD_BIAS = get_config("ADD_BIAS", true);
int NB_NEIGHBORS = get_config("NB_NEIGHBORS", -1);
int NB_EDGES = 0;

bool TRICK_LAMBDA_MUL = get_config("TRICK_LAMBDA_MUL", true);

int STAG_BUFFER_SIZE = get_config("STAG_BUFFER_SIZE", 20);

string ALGO = get_config_str("ALGO", "STAG");

float LEARNING_RATE = get_config("LEARNING_RATE", 0.02);

bool SHUFFLE_DATASET = get_config("SHUFFLE_DATASET", false);

bool EXACT_REGUL = get_config("EXACT_REGUL", false);

float NB_MESSAGES = get_config("MESSAGES", 10.0);
int ACCURACY = get_config("ACCURACY", 100);

int E_START = get_config("E_START", 1);
int E_END = get_config("E_END", T_MAX*NB_MESSAGES*N);


string PREFIX = get_config_str("PREFIX", "/");

bool B_DBG_W = get_config("DBG_W", false);
bool B_DBG_TEST_ERROR = get_config("DBG_TEST_ERROR", true);

bool B_UPDATE_ON_RECV = get_config("UPDATE_ON_RECV", false);



//////////
// DATA //
//////////

Matrix X; // Sample
Matrix y; // Labels

Matrix X_test; // Test sample
Matrix y_test; // Test labels


int t = 0;
int nbgradients_evaluated = 0;

//////////////////////



#include <time.h>
#include <unistd.h>
#include <sys/time.h>
static struct timeval ts;
bool tic(int ms) {
//	if(ms>0) {
//		struct timeval tv;
//		gettimeofday(&tv, 0);
//		float dt = (tv.tv_sec-ts.tv_sec)*1000 + 0.001*(tv.tv_usec-ts.tv_usec);
//		if(dt > ms) {
//			ts = tv;
//			return true;
//		}
//		return false;
//	} else {
		return (t)%((int)(ms))==0;
//	}
}
int TICTIC = 1;


///////////////////////

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
	Matrix y;
	Matrix w;
	Matrix w_avg;
	float b;
	int n;

	Matrix s;

	float cost;
	double weight;

	int nbmessages_sent;


	Node() {
		nbmessages_sent = NB_MESSAGES;
		shuffled_indices = 0;nbNeighbors = 0;neighbors = 0;isample = 0;
		this->id = __node_last_id++; b=cost=n=0;
		weight = 0;
	}
#include "funcs.cpp"


	void iteration() {
		if(N==1) {
			optimize();
		} else {
			if(nbmessages_sent>=NB_MESSAGES) {
				optimize_gossip();
				nbmessages_sent = 0;
			} else {
				last_receiver = gossip_choose_receiver();
				if(send(node[last_receiver])) nbmessages_sent++;
			}
		}
	}


	////////////////
	// ALGORITHMS //
	////////////////

	void SGD(float learningRate) {
		int i = draw_sample();
		float* sample = X.get_row(i);

		// Gradient for l2-regularized hinge loss

		if( hinge_loss(sample, y[i],w) > 0) {
			for(int d=0; d<D; d++) w[d] -= learningRate * (LAMBDA*w[d] - y[i]*sample[d]);
		} else {
			for(int d=0; d<D; d++) w[d] -= learningRate * LAMBDA * w[d];
		}

	}

	Matrix averagedGradient;
	Matrix gradientsMemory;
	int curbufsize;
	void SAG(float learningRate) {
		if(!gradientsMemory) {
			gradientsMemory.create(D, n); gradientsMemory.clear();
			averagedGradient.create(D, 1); averagedGradient.clear();
			curbufsize = 0;
		}

		int i = draw_sample();
		float* sample = X.get_row(i);

		// Update averaged gradient
		vector_sub_float(averagedGradient,gradientsMemory.get_row(i), D);
		if( hinge_loss(sample, y[i],w) > 0) {
			for(int d=0; d<D; d++) gradientsMemory[i*D + d] = - y[i]*sample[d];
		} else {
			for(int d=0; d<D; d++) gradientsMemory[i*D + d] = 0;
		}
		vector_add_float(averagedGradient,gradientsMemory.get_row(i), D);

		if(curbufsize < n) curbufsize++;

		// Learn
		for(int d=0; d<D; d++) w[d] -= learningRate*(LAMBDA*w[d] + averagedGradient[d]/n);
	}


	void SAG_gossip(float learningRate) {
		if(!gradientsMemory) {
			gradientsMemory.create(D, n); gradientsMemory.clear();
			averagedGradient.create(D, 1); averagedGradient.clear();
			curbufsize = 0;
		}

		int i = draw_sample();
		float* sample = X.get_row(i);

		// Update averaged gradient
		vector_sub_float(averagedGradient,gradientsMemory.get_row(i), D);
		if( hinge_loss(sample, y[i],w) > 0) {
			for(int d=0; d<D; d++) gradientsMemory[i*D + d] = - y[i]*sample[d];
		} else {
			for(int d=0; d<D; d++) gradientsMemory[i*D + d] = 0;
		}
		vector_add_float(averagedGradient,gradientsMemory.get_row(i), D);

		if(curbufsize < n) curbufsize++;

		// Learn
	//	for(int d=0; d<D; d++) w[d] *= (1 - learningRate * LAMBDA);
		for(int d=0; d<D; d++) s[d] = (1-learningRate*LAMBDA)*w[d] - learningRate/curbufsize*averagedGradient[d];
		for(int d=0; d<D; d++) w[d] = s[d]/weight;

	}

	void SGD_gossip(float learningRate) {
		int i = draw_sample();
		float* sample = X.get_row(i);

		// Gradient for l2-regularized hinge loss
		if( hinge_loss(sample, y[i],w) > 0) {
			for(int d=0; d<D; d++) s[d] -= learningRate * (LAMBDA*w[d] - y[i]*sample[d]);
		} else {
			for(int d=0; d<D; d++) s[d] -= learningRate * LAMBDA * w[d];
		}

		for(int d=0; d<D; d++) w[d] = s[d]/weight;
	}

	//////////
	// FLOW //
	//////////

	int iterations = 1;
	void optimize() {
		//SAG(LEARNING_RATE);
		SGD(0.1/pow(iterations,1/4.0));
		nbgradients_evaluated++;
		iterations++;
	}

	void optimize_gossip() {
		//SAG_gossip(LEARNING_RATE);
		nbgradients_evaluated++;
	}

	////////////
	// GOSSIP //
	////////////

	void init_gossip() {
		w = 0.0;
		s.create(D,1);
		s = 0.0;
		weight = 1.0;
	}

	int send(Node& node) {
		if(weight<10e-7) return 0;
		weight *= 0.5;
		for(int d=0; d<D; d++) s[d] *= 0.5;
		node.receive(id);
		return 1;
	}

	void receive(int sender) {
		weight += node[sender].weight;
		for(int d=0; d<D; d++) s[d] += node[sender].s[d];

		//if(B_UPDATE_ON_RECV && weight>1e-5)	for(int d=0; d<D; d++) w[d] = (1-LEARNING_RATE*LAMBDA)*w[d] - LEARNING_RATE*averagedGradient[d];
		for(int d=0; d<D; d++) w[d] = s[d]/weight;
	}

};


#include "misc.cpp"


void compute_errors() {
	if(B_DBG_TEST_ERROR) {
		double avgcost = 0;
		double cost2 = 0;
		int N = ::N > 10 ? 10 : ::N;
		for(int i=0; i<N; i++) node[i].compute_estimate();
		for(int i=0; i<N; i++) {
			avgcost += node[i].cost;
			cost2 += node[i].cost * node[i].cost;
		}

		avgcost /= N;
		cost2 /= N;

		//		int i= rand()%N;
		//		node[i].compute_estimate();
		//		double avgcost = node[i].cost;

		ffE << ((float)nbgradients_evaluated/::N) << " " << avgcost << "\n";
		//	ffEstddev << ((float)nbgradients_evaluated/::N) << " " << sqrt(cost2 - avgcost*avgcost) << "\n";
		ffE.flush();
		//	ffEstddev.flush();
		//	setenv("GSVM_E_", avgcost);
	}
	if(B_DBG_W) dump_classifier();
}


int main(int argc, char **argv) {
	srand(clock());

	try {
		DBG_START("Init ");
		init();

		sys("mkdir -p data");

		DBG_END();

		////////////////////////

		if(N!=1) {for(int i=0; i<N; i++) node[i].init_gossip();}

		t = 0;
		compute_errors();

		DBG("start");
		for(t=1; t/N<T_MAX; t++) {
			last_sender = gossip_choose_sender();
			node[last_sender].iteration();
			if(tic(TICTIC)) {
				FILE* f = fopen("TICTIC", "r"); if(f) {int r = fscanf(f, "%d", &TICTIC); (void)r; fclose(f);}
				DBGV(TICTIC);
				DBG("t=" << ((float)nbgradients_evaluated/N));
				compute_errors();
			}
		}

		//////////////////////

		DBG("finished");
	}catch(const char* c) {DBG("ERROR : " << c );}
}
