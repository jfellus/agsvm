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

int T_MAX = get_config("T_MAX", 100000000);
bool ADD_BIAS = get_config("ADD_BIAS", true);
int NB_NEIGHBORS = get_config("NB_NEIGHBORS", -1);
int NB_EDGES = 0;

bool TRICK_LAMBDA_MUL = get_config("TRICK_LAMBDA_MUL", true);

int STAG_BUFFER_SIZE = get_config("STAG_BUFFER_SIZE", 20);

//////////
// DATA //
//////////

Matrix X; // Sample
Matrix y; // Labels

Matrix X_test; // Test sample
Matrix y_test; // Test labels

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
	Matrix y;
	Matrix w;
	Matrix w_avg;
	float b;
	int n;

	float cost;

	int iterations;


	Node() {this->id = __node_last_id++; iterations = 1; }
	void init(Matrix& X, Matrix& y, int first, int n) {
		this->y.create_ref(&y[first], 1, n);
		this->X.create_ref(&X[first*X.width],X.width,n);

		w.create(D, 1); w.clear();
		w_avg.create(D, 1); w.clear();
		b = 0;

		this->n = n;
	}

	~Node() {}


	void iteration() {
		optimize();

		// Send
		last_receiver = gossip_choose_receiver();
		send(node[last_receiver]);
	}


	// Single hinge loss = max(0,1-y<w,x>)
	float hinge_loss(float* x, float y, float* w) {
		float r = 1-y*vector_ps_float(w,x,D);
		return MAX(r,0);
	}

	// Overall Hinge loss = sum_{x,y} max(0,1-y<w,x>)
	float hinge_loss() {
		float E = 0;
		for(int i=0; i< ::n; i++) E += hinge_loss(::X.get_row(i), ::y[i], w);
		return E / ::n;
	}

	// SVM cost function (l2 regularization + hinge loss)
	float compute_cost() {
		return LAMBDA*0.5*vector_n2p2_float(w, D) +  hinge_loss();
	}

	float compute_cost_avg_w() {
		return  LAMBDA*0.5*vector_n2p2_float(w_avg, D)  +  hinge_loss();
	}

	////////////////
	// ALGORITHMS //
	////////////////

	void SGD(float learningRate) {
		int i = rand()%n;
		float* sample = X.get_row(i);

		// Gradient for l2-regularized hinge loss
		if(TRICK_LAMBDA_MUL) trick_lambda_mul(sample, i);
		else {
			if( y[i] * vector_ps_float(w, sample, D) < 1) {
				for(int d=0; d<D; d++) w[d] -= learningRate * (LAMBDA*w[d] - y[i]*sample[d]);
			} else {
				for(int d=0; d<D; d++) w[d] -= learningRate * LAMBDA * w[d];
			}
		}
	}

	float mul = 1;
	void trick_lambda_mul(float* sample, int i) {
		if( y[i] * vector_ps_float(w, sample, D) * mul < 1) {
			for(int d=0; d<D; d++) w[d] = w[d]*mul - learningRate * (LAMBDA*w[d] - y[i]*sample[d]);
			mul = 1;
		} else {
			mul *= (1 - learningRate * LAMBDA);
		}
	}

	void project_on_L2_ball() {
		float n2p2 = vector_n2p2_float(w, D);
		if(n2p2 > 1.0/LAMBDA) {
			vector_smul_float(w, 1.0/sqrtf(n2p2) * 1.0/sqrtf(LAMBDA), D );
		}
	}

	void pegasos() {
		SGD(1.0/(LAMBDA*iterations)); // Except learning rate, its a classical SGD (projection ?)
		//project_on_L2_ball();
	}


	Matrix gradientsMemory;
	Matrix averagedGradient;
	void SAG(float learningRate) {
		if(!gradientsMemory) {
			gradientsMemory.create(n, 1); gradientsMemory.clear();
			averagedGradient.create(D, 1); averagedGradient.clear();
		}

		int i = rand()%n;
		float* sample = X.get_row(i);

		// Update averaged gradient
		for(int d=0; d<D; d++) averagedGradient[d] -= gradientsMemory[i] * y[i] * sample[d];
		gradientsMemory[i] = (y[i]*vector_ps_float(w, sample, D)) < 1 ? -1 : 0;
		for(int d=0; d<D; d++) averagedGradient[d] += gradientsMemory[i] * y[i] * sample[d];

		// Learn
		for(int d=0; d<D; d++) w[d] *= (1 - learningRate * LAMBDA);
		for(int d=0; d<D; d++) w[d] -= learningRate / n * averagedGradient[d];
	}

	int curi;
	void STAG(float learningRate) {
		if(!gradientsMemory) {
			curi = 0;
			gradientsMemory.create(D, STAG_BUFFER_SIZE); gradientsMemory.clear();
			averagedGradient.create(D, 1); averagedGradient.clear();
		}

		int i = rand()%n;
		float* sample = X.get_row(i);

		int lasti = (curi+STAG_BUFFER_SIZE)%STAG_BUFFER_SIZE;

		// Update averaged gradient
		for(int d=0; d<D; d++) averagedGradient[d] -= gradientsMemory[lasti*D+d];
		if( y[i] * vector_ps_float(w, sample, D) < 1) {
			for(int d=0; d<D; d++) gradientsMemory[curi*D + d] = (LAMBDA*w[d] - y[i]*sample[d]);
		} else {
			for(int d=0; d<D; d++) gradientsMemory[curi*D + d] = LAMBDA * w[d];
		}
		for(int d=0; d<D; d++) averagedGradient[d] += gradientsMemory[curi*D+d];

		curi = (curi+1)%STAG_BUFFER_SIZE;

		// Learn
		for(int d=0; d<D; d++) w[d] -= learningRate / n * averagedGradient[d];
	}

	void DA(float learningRate) {
		if(!averagedGradient) {averagedGradient.create(D, 1); averagedGradient.clear();}

		int i = rand()%n;
		float* sample = X.get_row(i);

		if( y[i] * vector_ps_float(w, sample, D) < 1) {
			for(int d=0; d<D; d++) averagedGradient[d] += (LAMBDA*w[d] - y[i]*sample[d]);
		} else {
			for(int d=0; d<D; d++) averagedGradient[d] += LAMBDA * w[d];
		}

		//dump(averagedGradient, D);

		// Learn
		for(int d=0; d<D; d++) w[d] -= learningRate * averagedGradient[d];
	}

	//////////
	// FLOW //
	//////////

	void optimize() {

		// Choose algorithm here
			// pegasos();
		 	 //SGD(0.00001);
			// SAG(0.1);
			 STAG(0.1);
			// DA(0.1/(LAMBDA*iterations));
		iterations++;

		// Compute average w (for Dual Averaging etc...)
		for(int d = 0; d < D; d ++) w_avg[d] = w_avg[d] * (iterations-1.0)/iterations + 1.0/iterations * w[d];
	}


	////////////
	// GOSSIP //
	////////////

	void init_gossip() {
	}

	int send(Node& node) {
		node.receive(0);
		return 0;
	}

	void receive(int arg) {
		// TODO
	}


	void compute_estimate() {
		this->cost = compute_cost();
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

	void dump(float* v, int D) {
		printf("[");
		for(int i=0; i<D; i++) printf(" %f", v[i]);
		printf(" ]\n");
		sleep(1);
	}
};


#include "misc.cpp"



//////////
// DUMP //
//////////

string fE;

void dump_classifier() {
	FILE* f = fopen(fmt("data/classifier_%u.txt", t), "w");
	fprintf(f, "%f %f %f", node[0].w[0], node[0].w[1], node[0].b);
	fclose(f);
}

void compute_errors() {
	if(last_sender!=-1) node[last_sender].compute_estimate();
	if(last_receiver!=-1) node[last_receiver].compute_estimate();

	DBG("Error = " << node[last_sender].cost);
	fappend(fE, fmt("%u %f\n", t, node[last_sender].cost));

	//dump_classifier();
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

	//	system("rm -rf data/*");
	system("rm -rf plots/*");
	fE = newfilename("data/E_%u.txt");


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
			X_test_nobias.load(dataset.c_str());
			for(int i=0; i<X_test.height; i++) {
				memcpy(X_test.get_row(i),X_test_nobias.get_row(i),X_test_nobias.width*sizeof(float));
				X_test[i*X_test.width + X_test_nobias.width] = 1;
			}
		}
		else X_test.load(dataset_test.c_str());
		y_test.load(labels_test.c_str());
		if(CATEGORY != -1) {
			for(int i=0; i<y.height; i++) y[i] = y[i]==CATEGORY ? 1 : -1;
		}
	}
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
		init();
	} else {
		DBG("with : " << argv[1]);
		dataset = argv[1];
		init();
		//chdir(dirname(argv[1]));
	}

	if(system("mkdir -p data")) {}

	DBG_END();

	////////////////////////

	for(int i=0; i<N; i++) node[i].init_gossip();

	for(t=1; t<T_MAX; t++) {
		DBGV(t);

		last_sender = gossip_choose_sender();
		node[last_sender].iteration();
		if(t % 100 == 0)
			compute_errors();
	}

	//////////////////////

	DBG("finished");
	}catch(const char* c) {DBG("ERROR : " << c );}
}
