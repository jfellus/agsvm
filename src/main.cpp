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

int T_MAX = get_config("T_MAX", 20000000*N);
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


string PREFIX = get_config_str("PREFIX", "");




//////////
// DATA //
//////////

Matrix X; // Sample
Matrix y; // Labels

Matrix X_test; // Test sample
Matrix y_test; // Test labels

float learningRate = 1;


int t = 0;
int nbgradients_evaluated = 0;

//////////////////////



#include <time.h>
#include <unistd.h>
#include <sys/time.h>
static struct timeval ts;
bool tic(int ms) {
	if(ms>0) {
	struct timeval tv;
	gettimeofday(&tv, 0);
	float dt = (tv.tv_sec-ts.tv_sec)*1000 + 0.001*(tv.tv_usec-ts.tv_usec);
	if(dt > ms) {
		ts = tv;
		return true;
	}
	return false;
	} else {
		return (nbgradients_evaluated)%(-ms*N)==0;
	}
}
int TICTIC = 1000;


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

	float cost;
	double weight;

	int iterations ;


	Node() {
		nbmessages_sent = 0;
		curbufsize = 0;
		mul = 1;
		curi = 0; shuffled_indices = 0;nbNeighbors = 0;neighbors = 0;isample = 0;
		this->id = __node_last_id++; iterations = 3; b=cost=n=0;
		weight = 0;
	}

	void init(Matrix& X, Matrix& y, int first, int n) {
		this->y.create_ref(&y[first], 1, n);
		this->X.create_ref(&X[first*X.width],X.width,n);

		w.create(D, 1); w.clear();
		w_avg.create(D, 1); w.clear();
		b = 0;

		this->n = n;
		if(SHUFFLE_DATASET) {
			shuffled_indices = new int[n];
			shuffle_dataset();
		}
	}

	~Node() {}

	void shuffle_dataset() {
		for(int i=0; i<n; i++) { shuffled_indices[i] = i;}
		for (int i = n - 1; i>=1; i--) {
			int j = rand()%i;
			int tmp = shuffled_indices[j];
			shuffled_indices[j]  = shuffled_indices[i];
			shuffled_indices[i]  = tmp;
		}
		isample = 0;
	}

	int nbmessages_sent;
	void iteration() {
		if(N==1) {
			optimize();
		} else {
			if(NB_MESSAGES<1) {
				if(nbmessages_sent>=1/NB_MESSAGES) {
					last_receiver = gossip_choose_receiver();
					send(node[last_receiver]);
					nbmessages_sent=0;
				} else {
					optimize_gossip();
					nbgradients_evaluated++;
					nbmessages_sent++;
				}
			} else {
				if(nbmessages_sent>=NB_MESSAGES) {
					optimize_gossip();
					nbgradients_evaluated++;
					nbmessages_sent = 0;
				} else {
					last_receiver = gossip_choose_receiver();
					if(send(node[last_receiver])) nbmessages_sent++;
				}
			}
		}
	}


	inline float eval_classifier(float* w, float* x) {
		return vector_ps_float(w,x,D);
	}

	// Single hinge loss = max(0,1-yf(x))
	inline float hinge_loss(float* x, float y, float* w) {
		float r = 1-y*eval_classifier(w,x);
		return MAX(r,0);
	}

	// Overall Hinge loss = sum_{x,y} max(0,1-y<w,x>)
	inline float hinge_loss(float* w) {
		float E = 0;
		for(int i=0; i< ::n; i++) E += hinge_loss(::X.get_row(i), ::y[i], w);
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
		float* sample = X.get_row(i);

		// Gradient for l2-regularized hinge loss
//		if(TRICK_LAMBDA_MUL) trick_lambda_mul(sample, i);
//		else {
			if( hinge_loss(sample, y[i],w) > 0) {
				for(int d=0; d<D; d++) w[d] -= learningRate * (LAMBDA*w[d] - y[i]*sample[d]);
			} else {
				for(int d=0; d<D; d++) w[d] -= learningRate * LAMBDA * w[d];
			}
//		}
	}

	float mul;
	void trick_lambda_mul(float* sample, int i) {
		if( y[i] * vector_ps_float(w, sample, D) * mul < 1) {
			for(int d=0; d<D; d++) w[d] = w[d]*mul - learningRate * (LAMBDA*w[d] - y[i]*sample[d]);
			mul = 1;
		} else {
			mul *= (1 - learningRate * LAMBDA);
		}
	}


	void pegasos() {
		SGD(1/(LAMBDA*iterations)); // Except learning rate, its a classical SGD (projection ?)
	}

	int* shuffled_indices;
	int isample;
	int draw_sample() {
		if(SHUFFLE_DATASET) {
			int i=shuffled_indices[isample++];
			if(isample>=n) shuffle_dataset();
			return i;
		}
		return rand()%n;
	}


	Matrix gradientsMemory;
	Matrix averagedGradient;
	void SAG(float learningRate) {
		if(!gradientsMemory) {
			//gradientsMemory.create(n, 1); gradientsMemory.clear();
			gradientsMemory.create(D, n); gradientsMemory.clear();
			averagedGradient.create(D, 1); averagedGradient.clear();
		}

		int i = draw_sample();
		float* sample = X.get_row(i);

		// Update averaged gradient
//		for(int d=0; d<D; d++) averagedGradient[d] -= gradientsMemory[i] * y[i] * sample[d];
//		gradientsMemory[i] = (y[i]*vector_ps_float(w, sample, D)) < 1 ? -1 : 0;
//		for(int d=0; d<D; d++) averagedGradient[d] += gradientsMemory[i] * y[i] * sample[d];

		// Update averaged gradient
		vector_sub_float(averagedGradient,gradientsMemory.get_row(i), D);
		if( hinge_loss(sample, y[i],w) > 0) {
			for(int d=0; d<D; d++) gradientsMemory[i*D + d] = (LAMBDA*w[d] - y[i]*sample[d]);
		} else {
			for(int d=0; d<D; d++) gradientsMemory[i*D + d] = LAMBDA * w[d];
		}
		vector_add_float(averagedGradient,gradientsMemory.get_row(i), D);
		if(curbufsize < n) curbufsize++;

		// Learn
	//	for(int d=0; d<D; d++) w[d] *= (1 - learningRate * LAMBDA);
		for(int d=0; d<D; d++) w[d] -= learningRate / curbufsize * averagedGradient[d];
	}

	void SAG_exact_regul(float learningRate) {
		if(!gradientsMemory) {
			//gradientsMemory.create(n, 1); gradientsMemory.clear();
			gradientsMemory.create(D, n); gradientsMemory.clear();
			averagedGradient.create(D, 1); averagedGradient.clear();
		}

		int i = draw_sample();
		float* sample = X.get_row(i);

		// Update averaged gradient
//		for(int d=0; d<D; d++) averagedGradient[d] -= gradientsMemory[i] * y[i] * sample[d];
//		gradientsMemory[i] = (y[i]*vector_ps_float(w, sample, D)) < 1 ? -1 : 0;
//		for(int d=0; d<D; d++) averagedGradient[d] += gradientsMemory[i] * y[i] * sample[d];

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
		for(int d=0; d<D; d++) w[d] = (1-learningRate*LAMBDA)*w[d] - learningRate/curbufsize*averagedGradient[d];
	}


	void SAG_gossip_exact_regul(float learningRate) {
		if(!gradientsMemory) {
			//gradientsMemory.create(n, 1); gradientsMemory.clear();
			gradientsMemory.create(D, n); gradientsMemory.clear();
			weight = 0;
		}

		if(weight>0.00001) for(int d=0; d<D; d++) w[d] = (1-learningRate*LAMBDA)*w[d] - learningRate/weight*averagedGradient[d];

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

		if(curbufsize < n) { curbufsize++; weight++; }
	}

	void SAG_gossip(float learningRate) {
		SAG_gossip_exact_regul(learningRate);
	}

	int curi;
	int curbufsize;
	void STAG(float learningRate) {
		if(!gradientsMemory) {
			curbufsize=0;
			curi = 0;
			gradientsMemory.create(D, STAG_BUFFER_SIZE); gradientsMemory.clear();
			averagedGradient.create(D, 1); averagedGradient.clear();
		}

		int i = draw_sample();
		float* sample = X.get_row(i);

		int lasti = (curi+STAG_BUFFER_SIZE+1)%STAG_BUFFER_SIZE;

		// Update averaged gradient
		for(int d=0; d<D; d++) averagedGradient[d] -= gradientsMemory[lasti*D+d];
		if( hinge_loss(sample, y[i],w) > 0) {
			for(int d=0; d<D; d++) gradientsMemory[curi*D + d] = (LAMBDA*w[d] - y[i]*sample[d]);
		} else {
			for(int d=0; d<D; d++) gradientsMemory[curi*D + d] = LAMBDA * w[d];
		}
		for(int d=0; d<D; d++) averagedGradient[d] += gradientsMemory[curi*D+d];

		curi = (curi+1)%STAG_BUFFER_SIZE;
		if(curbufsize < STAG_BUFFER_SIZE) curbufsize++;

		// Learn
		for(int d=0; d<D; d++) w[d] -= averagedGradient[d] * learningRate / curbufsize;
	}

	Matrix indicesMemory;
	void STAG_exact_regul(float learningRate) {
		if(!indicesMemory) {
			curbufsize=0;
			curi = 0;
			indicesMemory.create(1, STAG_BUFFER_SIZE); indicesMemory.clear();
			for(uint i = 0; i<STAG_BUFFER_SIZE; i++) indicesMemory[i] = -1;
			averagedGradient.create(D, 1); averagedGradient.clear();
		}

		int i = draw_sample();
		float* sample = X.get_row(i);

		int lasti = (curi+STAG_BUFFER_SIZE+1)%STAG_BUFFER_SIZE;

		// Update averaged gradient
		if(indicesMemory[lasti]>=0) for(int d=0; d<D; d++) averagedGradient[d] += y[indicesMemory[lasti]]*X.get_row(indicesMemory[lasti])[d];
		if( hinge_loss(sample, y[i],w) > 0) {
			indicesMemory[curi] = i;
			for(int d=0; d<D; d++) averagedGradient[d] -= y[i]*X.get_row(i)[d];
		} else {
			indicesMemory[curi] = -1;
		}

		curi = (curi+1)%STAG_BUFFER_SIZE;
		if(curbufsize < STAG_BUFFER_SIZE) curbufsize++;

		// Learn
		for(int d=0; d<D; d++) w[d] -= learningRate*(LAMBDA*w[d] + averagedGradient[d]/curbufsize);
	}


	void STAG_gossip(float learningRate) {
		if(!gradientsMemory) {
			curbufsize=0;
			curi = 0;
			gradientsMemory.create(D, STAG_BUFFER_SIZE); gradientsMemory.clear();
		}

		// Learn
		if(weight>0.00001) for(int d=0; d<D; d++) w[d] = (1-learningRate*LAMBDA)*w[d] - learningRate/weight*averagedGradient[d];

		int i = draw_sample();
		float* sample = X.get_row(i);

		int lasti = (curi+STAG_BUFFER_SIZE+1)%STAG_BUFFER_SIZE;

		// Update averaged gradient
		for(int d=0; d<D; d++) averagedGradient[d] -= gradientsMemory[lasti*D+d];
		if( hinge_loss(sample, y[i],w) > 0) {
			for(int d=0; d<D; d++) gradientsMemory[curi*D + d] = - y[i]*sample[d];
		} else {
			for(int d=0; d<D; d++) gradientsMemory[curi*D + d] = 0;
		}
		for(int d=0; d<D; d++) averagedGradient[d] += gradientsMemory[curi*D+d];

		curi = (curi+1)%STAG_BUFFER_SIZE;
		if(curbufsize < STAG_BUFFER_SIZE) {
			curbufsize++;
			weight++;
		}
	}

	void DA(float learningRate) {
		if(!averagedGradient) {averagedGradient.create(D, 1); averagedGradient.clear();}

		int i = draw_sample();
		float* sample = X.get_row(i);

		if( hinge_loss(sample, y[i],w) > 0) {
			for(int d=0; d<D; d++) averagedGradient[d] += learningRate*(LAMBDA*w[d] - y[i]*sample[d]);
		} else {
			for(int d=0; d<D; d++) averagedGradient[d] += learningRate*LAMBDA * w[d];
		}

		//dump(averagedGradient, D);

		// Learn
		for(int d=0; d<D; d++) w[d] = learningRate * averagedGradient[d];
	}

	//////////
	// FLOW //
	//////////

	void optimize() {

		// Choose algorithm here
		if(ALGO=="STAG") {
			if(EXACT_REGUL) STAG_exact_regul(LEARNING_RATE);
			else STAG(LEARNING_RATE);
		}
		else if(ALGO=="SAG") {
			if(EXACT_REGUL) SAG_exact_regul(LEARNING_RATE);
			else SAG(LEARNING_RATE);
		}
		else if(ALGO=="DA") DA( LEARNING_RATE /(LAMBDA*iterations) );
		else if(ALGO=="pegasos") pegasos();
		else if(ALGO=="SGD") SGD(LEARNING_RATE);

		iterations++;

		// Compute average w (for Dual Averaging etc...)
		//for(int d = 0; d < D; d ++) w_avg[d] = w_avg[d] * (iterations-1.0)/iterations + 1.0/iterations * w[d];
	}

	void optimize_gossip() {
		if(ALGO=="STAG") {
			STAG_gossip(LEARNING_RATE);
		} else {
			SAG_gossip(LEARNING_RATE);
		}
	}

	////////////
	// GOSSIP //
	////////////

	void init_gossip() {
		averagedGradient.create(D, 1); averagedGradient.clear();
		weight = 0;
	}

	int send(Node& node) {
		//if(weight<10e-7) return 0;
		weight *= 0.5;
		for(int d=0; d<D; d++) averagedGradient[d] *= 0.5;
		node.receive(id);
		return 1;
	}

	void receive(int sender) {
		weight += node[sender].weight;
		for(int d=0; d<D; d++) averagedGradient[d] += node[sender].averagedGradient[d];
	}


	void compute_estimate() {
		this->cost = //ALGO=="DA" ? compute_cost_avg_w() :
				compute_cost();
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

string fE, fEstddev;
std::ofstream ffE, ffEstddev;

void dump_classifier() {
	static int i=0;
	shell(TOSTRING("mkdir -p data/w/" << ALGO << "/N" << N << "/l" << LEARNING_RATE << "/"));
	node[0].w.write(TOSTRING("data/w/" << ALGO << "/N" << N << "/l" << LEARNING_RATE << "/" << fmt("%05d",i++) << ".fvec").c_str());
}

void compute_errors() {

	if(N==1) {
		for(int i=0; i<N; i++) node[i].compute_estimate();
		fappend(fE, fmt("%u %f\n", t, node[0].cost));
		setenv("GSVM_E_", node[0].cost);
	} else {
		double avgcost = 0;
		double cost2 = 0;
		int N = ::N > 1 ? 1 : ::N;
		for(int i=0; i<N; i++) node[i].compute_estimate();
		for(int i=0; i<N; i++) {
			avgcost += node[i].cost;
			cost2 += node[i].cost * node[i].cost;
		}

		avgcost /= N;
		cost2 /= N;

		ffE << ((float)nbgradients_evaluated/::N) << " " << avgcost << "\n";
	//	ffEstddev << ((float)nbgradients_evaluated/::N) << " " << sqrt(cost2 - avgcost*avgcost) << "\n";
		ffE.flush();
	//	ffEstddev.flush();
	//	setenv("GSVM_E_", avgcost);
	}
	dump_classifier();
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

	//	system("rm -rf data/*");
	shell("rm -rf plots/*");
	if(ALGO=="STAG" || ALGO=="STAGR") {
		if(N==1) {
		fE = newfilename(fmt("data/E_%s_N%u_%u_%f_%%u.txt", ALGO.c_str(),N, STAG_BUFFER_SIZE, LEARNING_RATE));
		fEstddev = newfilename(fmt("data/DEV_%s_%u_%f_%%u.txt", ALGO.c_str(), STAG_BUFFER_SIZE, LEARNING_RATE));
		} else {
			shell(fmt("mkdir -p data/N%u", N));
			shell(fmt("mkdir -p data/N%u/stddev", N));
			fE = newfilename(fmt("data/N%u/E_%s_N%u_%u_%f_%%u.txt", N, ALGO.c_str(), N, STAG_BUFFER_SIZE, LEARNING_RATE));
			fEstddev = newfilename(fmt("data/N%u/stddev/DEV_%s_N%u_%u_%f_%%u.txt", N, ALGO.c_str(), N, STAG_BUFFER_SIZE, LEARNING_RATE));
		}
	}
	else {
		fE = newfilename(fmt("data/E_%s_N%u_%f_M%f_%%u.txt", ALGO.c_str(), N, LEARNING_RATE, NB_MESSAGES));
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
//		if(CATEGORY != -1) {
//			for(int i=0; i<y.height; i++) y[i] = y[i]==CATEGORY ? 1 : -1;
//		}
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


int main(int argc, char **argv) {
	srand(clock());

	try {
	DBG_START("Init ");
	if(argc<=1) {
		init();
	} else {
//		DBG("with : " << argv[1]);
//		dataset = argv[1];
//		LEARNING_RATE = atof(argv[1]);
		PREFIX = argv[1];
		init();
		//chdir(dirname(argv[1]));
	}

	if(system("mkdir -p data")) {}

	DBG_END();

	////////////////////////

	if(N!=1) {for(int i=0; i<N; i++) node[i].init_gossip();}
//
//	setenv("GSVM_N_", N);
//	setenv("GSVM_STAG_BUFFER_SIZE_", STAG_BUFFER_SIZE);
//	setenv("GSVM_LEARNING_RATE_", LEARNING_RATE);


	t = 0;
	compute_errors();
	for(t=1; t/N<T_MAX; t++) {
		last_sender = gossip_choose_sender();
		node[last_sender].iteration();
		if(tic(TICTIC)) {
			FILE* f = fopen("TICTIC", "r"); if(f) {fscanf(f, "%u", &TICTIC); fclose(f);}
			DBG("t=" << ((float)nbgradients_evaluated/N));
			compute_errors();
		}
	}

	//////////////////////

	DBG("finished");
	}catch(const char* c) {DBG("ERROR : " << c );}
}
