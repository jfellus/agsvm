/*
 * main.cpp
 *
 *  Created on: 8 nov. 2013
 *      Author: jfellus
 */
#define MONOTHREAD

#include <stdio.h>
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
#include "vector_sparse/MatrixSparse.h"

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

string ALGO = get_config_str("ALGO", "STAG");

float LEARNING_RATE = get_config("LEARNING_RATE", 0.02);

bool SHUFFLE_DATASET = get_config("SHUFFLE_DATASET", false);




//////////
// DATA //
//////////

MatrixSparse X; // Sample
Matrix y; // Labels

MatrixSparse X_test; // Test sample
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

	MatrixSparse X;
	Matrix y;
	Matrix w;
	Matrix w_avg;
	float b;
	int n;

	float cost;

	int iterations;


	Node() {
		curi = 0; shuffled_indices = 0;nbNeighbors = 0;neighbors = 0;isample = 0;
		this->id = __node_last_id++; iterations = 1; b=cost=n=0;
	}

	void init(MatrixSparse& X, Matrix& y, int first, int n) {
		this->y.create_ref(&y[first], 1, n);
		this->X.create_ref(X, first,n);

		w.create(D, 1);
		w_avg.create(D, 1);
		b = 0;

		this->n = n;
		if(SHUFFLE_DATASET) {
			shuffled_indices = new int[n];
			for(int i=0; i<n; i++) { shuffled_indices[i] = i;}
			for (int i = n - 1; i>=1; i--) {
				int j = rand()%i;
				int tmp = shuffled_indices[j];
				shuffled_indices[j]  = shuffled_indices[i];
				shuffled_indices[i]  = tmp;
			}
		}
	}

	~Node() {}


	void iteration() {
		DBG("optim");
		optimize();

		// Send
		last_receiver = gossip_choose_receiver();
		send(node[last_receiver]);
		DBG("ok");
	}


	inline float eval_classifier(const Matrix& w, const VectorSparse& x) {
		return vector_ps_float(w.data,x,D);
	}

	// Single hinge loss = max(0,1-yf(x))
	inline float hinge_loss(const VectorSparse& x, float y, const Matrix& w) {
		float r = 1-y*eval_classifier(w,x);
		return MAX(r,0);
	}

	// Overall Hinge loss = sum_{x,y} max(0,1-y<w,x>)
	inline float hinge_loss(const Matrix& w) {
		float E = 0;
		for(int i=0; i< ::n; i++) E += hinge_loss(::X.get_row(i), ::y[i], w);
		DBGV(E);
		return E / ::n;
	}

	inline float l2(float* w) { return vector_n2p2_float(w, D); }

	// SVM cost function (l2 regularization + hinge loss)
	inline float compute_cost() {
		return LAMBDA*0.5*l2(w) +  hinge_loss(w);
	}

	inline float compute_cost_avg_w() {
		return  LAMBDA*0.5*l2(w_avg)  +  hinge_loss(w_avg);
	}

	////////////////
	// ALGORITHMS //
	////////////////

	void SGD(float learningRate) {
		int i = draw_sample();
		VectorSparse& sample = X.get_row(i);

		// Gradient for l2-regularized hinge loss
		if( hinge_loss(sample, y[i],w) > 0) {
			for(uint d=0;d<D; d++)	w[d] -= (learningRate * LAMBDA)*w[d];
			for(auto j=sample.entries.begin(); j!=sample.entries.end(); j++) {
				w[(*j).i] += (learningRate * y[i]) * (*j).val ;
			}
		} else for(uint d=0;d<D; d++)	w[d] -= (learningRate * LAMBDA)*w[d];

	}


	void pegasos() {
		SGD(LEARNING_RATE/(LAMBDA*iterations)); // Except learning rate, its a classical SGD (projection ?)
	}

	int* shuffled_indices;
	int isample;
	int draw_sample() {
		if(SHUFFLE_DATASET) {
			int i=shuffled_indices[isample++];
			if(isample>=n) isample = 0;
			return i;
		}
		return rand()%n;
	}


	MatrixSparse gradientsMemory;
	VectorSparse averagedGradient;
	void SAG(float learningRate) {
		if(!gradientsMemory) {
			//gradientsMemory.create(n, 1); gradientsMemory.clear();
			gradientsMemory.create(D, n); gradientsMemory.clear();
		}

		int i = draw_sample();
		VectorSparse& sample = X.get_row(i);

		// Update averaged gradient
		//		for(int d=0; d<D; d++) averagedGradient[d] -= gradientsMemory[i] * y[i] * sample[d];
		//		gradientsMemory[i] = (y[i]*vector_ps_float(w, sample, D)) < 1 ? -1 : 0;
		//		for(int d=0; d<D; d++) averagedGradient[d] += gradientsMemory[i] * y[i] * sample[d];

		// Update averaged gradient
		averagedGradient -= gradientsMemory[i];
		if( hinge_loss(sample, y[i],w) > 0) {
			gradientsMemory[i] = w;
			gradientsMemory[i] *= LAMBDA;
			gradientsMemory[i] -= (y[i]*sample);
		}
		else {
			gradientsMemory[i] = w;
			gradientsMemory[i] *= LAMBDA;
		};
		averagedGradient += gradientsMemory[i];

		// Learn
		w *= (1 - learningRate * LAMBDA);
		w -= (learningRate / n) * averagedGradient;
	}

	int curi;
	void STAG(float learningRate) {
		if(!gradientsMemory) {
			curi = 0;
			gradientsMemory.create(D, STAG_BUFFER_SIZE); gradientsMemory.clear();
		}

		int i = draw_sample();
		VectorSparse& sample = X.get_row(i);

		int lasti = (curi+STAG_BUFFER_SIZE-1)%STAG_BUFFER_SIZE;

		// Update averaged gradient
		averagedGradient -= gradientsMemory[lasti];
		if( hinge_loss(sample, y[i],w) > 0)  {
			gradientsMemory[curi] = w;
			gradientsMemory[curi] *= LAMBDA;
			gradientsMemory[curi] -= y[i]*sample;
		}
		else {
			gradientsMemory[curi] = w;
			gradientsMemory[curi] *= LAMBDA;
		}

		averagedGradient += gradientsMemory[curi];

		curi = (curi+1)%STAG_BUFFER_SIZE;

		// Learn
		w -= (learningRate / n) * averagedGradient;
	}

	void DA(float learningRate) {
		int i = draw_sample();
		VectorSparse& sample = X.get_row(i);

		if( hinge_loss(sample, y[i],w) > 0) {
			averagedGradient += (LAMBDA*w);
			averagedGradient -= y[i]*sample;
		}
		else averagedGradient += (LAMBDA * w);

		//dump(averagedGradient, D);

		// Learn
		w = averagedGradient; w *= learningRate;
	}

	//////////
	// FLOW //
	//////////

	void optimize() {

		DBG(ALGO << "(lr=" << LEARNING_RATE << ",bufsize=" << STAG_BUFFER_SIZE << ",n=" << n << ",D=" << D << ")");
		// Choose algorithm here
		if(ALGO=="STAG") STAG(LEARNING_RATE);
		else if(ALGO=="SAG") SAG(LEARNING_RATE);
		else if(ALGO=="DA") DA( LEARNING_RATE /(LAMBDA*iterations) );
		else if(ALGO=="pegasos") pegasos();
		else if(ALGO=="SGD") SGD(LEARNING_RATE);
		iterations++;

		// Compute average w (for Dual Averaging etc...)
		w_avg *= (iterations-1.0)/iterations;
		w_avg += (float*)((1.0f/iterations) * w);

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
		this->cost = ALGO=="DA" ? compute_cost_avg_w() : compute_cost();
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


#include "misc_sparse.cpp"



//////////
// DUMP //
//////////

string fE;

void dump_classifier() {
	FILE* f = fopen(fmt("data/classifier_%u.txt", t), "w");
	//fprintf(f, "%f %f %f", node[0].w[0], node[0].w[1], node[0].b);
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
	setlocale(LC_ALL, "C");
	DBG("INIT");
	DBGV(NBTHREADS);

	//	system("rm -rf data/*");
	shell("rm -rf plots/*");
	if(ALGO=="STAG") fE = newfilename(fmt("data/E_%s_%u_%f_%%u.txt", ALGO.c_str(), STAG_BUFFER_SIZE, LEARNING_RATE));
	else fE = newfilename(fmt("data/E_%s_%f_%%u.txt", ALGO.c_str(), LEARNING_RATE));


	DBGV(dataset);
	DBGV(labels);
	DBGV(dataset_test);
	DBGV(labels_test);

	if(ADD_BIAS) {
		X.load(dataset.c_str());
		for(int i=0; i<X.height; i++) X.set(i,X.width,1);
	}
	else
		X.load(dataset.c_str());
	load_labels(y,labels.c_str());

	if(ADD_BIAS) {
		X_test.load(dataset.c_str());
		for(int i=0; i<X_test.height; i++) X_test.set(i, X_test.width,1);
	}
	else X_test.load(dataset_test.c_str());
	load_labels(y_test, labels_test.c_str());

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
	sleep(1);
}


int main(int argc, char **argv) {
	//	srand(time(0));

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

		t = 0;
		compute_errors();
		for(t=1; t<T_MAX; t++) {
			DBGV(t);

			last_sender = gossip_choose_sender();
			node[last_sender].iteration();
			if(t % 100 == 0) {
				DBG("compute error");
				compute_errors();
				DBG("error : " << node[last_sender].cost);
			}
		}

		//////////////////////

		DBG("finished");
	}catch(const char* c) {DBG("ERROR : " << c );}
}
