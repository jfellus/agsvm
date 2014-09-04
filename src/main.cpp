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
#include <vector>
#include <string>
#include <libgen.h>

using namespace std;


template <class T> T get_config(const char* what, T default_val) {
	FILE* f = fopen("config.properties", "r");
	char line[512];
	char* c = 0;
	T v = default_val;
	while ( fgets (line , 512 , f) != NULL ) {
	      if(!strncmp(what, line, strlen(what))) {
	    	  c = line + strlen(what);
	    	  while(*c==' ' || *c=='=') c++;
	    	  v = (T) atof(c);
	      }
	}
	fclose(f);
	return v;
}


////////////
// PARAMS //
////////////

int LIMIT_NDATA = get_config("LIMIT_NDATA", -1);

bool LATE_PCA = get_config("LATE_PCA", 0);
bool USE_ENERGY = get_config("USE_ENERGY", 0);
float KEEP_ENERGY = get_config("KEEP_ENERGY", 0.99);
int q = get_config("NB_PROJS", 10);

int NB_NEIGHBORS = get_config("NB_NEIGHBORS", -1);

int N = get_config("NB_NODES", 100);
int T_MAX = N*40;
int n;
int D;


int NB_EDGES = 0;


//////////
// DATA //
//////////

Matrix X;

Matrix global_mu;
Matrix FULL_COVARIANCE;
Matrix FULL_COVARIANCE_th;

int t = 0;
int last_sender=-1;
int last_receiver=-1;

int last_nb_projs_sent;
int	last_msg_size;
long total_msg_size = 0;

//Matrix mat_plotall;

//////////////////////




#include "algebra.cpp"
int __node_last_id = 0;

//////////
// NODE //
//////////

class Node;
Node* node;

class Node {
public:
	int id;

	Matrix X;
	int n;
	Matrix mu;
	Matrix C;
	Matrix C_rec;
	Matrix U,L;

	float nb_projs_kept;

	float w;

	Node() {this->id = __node_last_id++;n=0;nb_projs_kept=0;w=0;}
	void init(Matrix& X, int first, int n) {
		this->n = n;
		//ni[i]=rand()%(n/N*2-2)+1;
		this->X.create_ref(X.get_row(first),D,n);
		this->mu.create(D,1);
		this->U.create(D,D);
		this->L.create(D,1);
		this->C_rec.create(D,D);
		this->C.create(D,D);
	}

	~Node() {}


	// 1) LOCAL PCA

	void local_pca() {
			X.mean_row(mu);
			C.covariance(X);
			C.eigendecompose(U,L);
			if(USE_ENERGY) nb_projs_kept = discard_after_energy(U,L,KEEP_ENERGY);
			else nb_projs_kept = discard_after(U,L,q);
//			C.reconstruct(U,L);
		//	X.save(fmt("data/x/x_%d.txt", id));
	}


	// 2) GOSSIP

	void init_gossip() {
		w = n; 		mu *= w;	C *= w;
	}

	int send(Node& node) {
		w /= 2; 	mu /= 2;	C /= 2;

		C.eigendecompose(U,L);
		if(USE_ENERGY) nb_projs_kept = discard_after_energy(U,L,KEEP_ENERGY);
		else nb_projs_kept = discard_after(U,L,q);
		C.reconstruct(U,L);

		node.receive(U,L,mu,w);
		return nb_projs_kept;
	}

	void receive(Matrix& U2, Matrix& L2, Matrix& mu2, float w2) {
		Matrix cr(D,D);
		cr.reconstruct(U2,L2);

		w += w2;	mu += mu2;		C += cr;
	}


	int send_perfect(Node& node) {
		w /= 2; 		C/=2;		 mu/=2;
		node.receive_perfect(C, mu, w);
		return C.width;
	}

	void receive_perfect(Matrix& C2, Matrix& mu2, float w2) {
		w += w2; 		mu += mu2;		C += C2;
	}

	void compute_estimate() {
		C_rec = C;
		C_rec /= w;

		Matrix mumut(D,D);
		mumut.AAt(mu);
		mumut /= w*w;
		C_rec -= mumut;

		if(LATE_PCA) {
			C_rec.eigendecompose(U,L);
			discard_after(U,L,q);
			C_rec.reconstruct(U,L);
		}
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

__multithread__(local_pca) (int i) {	node[i].local_pca(); }

void gossip() {
	last_sender = gossip_choose_sender();
	last_receiver = node[last_sender].gossip_choose_receiver();

	if(LATE_PCA) last_nb_projs_sent = node[last_sender].send_perfect(node[last_receiver]);
	else last_nb_projs_sent = node[last_sender].send(node[last_receiver]);

	last_msg_size = (last_nb_projs_sent+1)*(D+1);
	total_msg_size += last_msg_size;

	fappend("data/stats/total_msg_size.txt",fmt("%d\n",total_msg_size));
	fappend("data/stats/nb_projs_sent.txt",fmt("%d\n",last_nb_projs_sent));
	fappend("data/stats/total_msg_size_N.txt",fmt("%f %f\n",(float)t/N, (float)total_msg_size/N));
}


//////////
// DUMP //
//////////


void plot_Lmean() {
	Matrix Lmean(D,1);
	for(int i=0; i<N; i++) Lmean += node[i].L;
	Lmean.normalize_l1();
	plot(Lmean, fmt("plots/spectrum_%08d.jpg",t), 0,1);
}

/*void plot_mat_plotall() {
	for(int i=0; i<N; i++) {
		int cellx = i%(int)(sqrt(N)+0.5);
		int celly = i/(int)(sqrt(N)+0.5);
		mat_plotall.set_part(cellx*(D+1), celly*(D+1), node[i].C_rec);
		mat_plotall.set((cellx+1)*(D+1)-1, celly*(D+1), 0);
		mat_plotall.clear_row((celly+1)*(D+1)-1);
	}
	mat_plotall.set_part(N%(int)(sqrt(N)+0.5)*(D+1), N/(int)(sqrt(N)+0.5)*(D+1), FULL_COVARIANCE);

	mat_plotall.save(fmt("data/mat_plotall/%08d.txt",t));
	plot(mat_plotall, fmt("plots/mat_plotall/%08d.jpg",t), -100, 100);
}*/



void compute_errors() {
	if(last_sender!=-1) node[last_sender].compute_estimate();
	if(last_receiver!=-1) node[last_receiver].compute_estimate();

	// Relative Error to Consensus
	float FCnorm = FULL_COVARIANCE.n2();
	float REC = 0;
	for(int i=0; i<N; i++) REC += FULL_COVARIANCE.l2(node[i].C_rec);
	REC /= (N*FCnorm);

/*	// Pair-Wise Error
	float PWE = 0;
	for(int i=0; i<N; i++) for(int j=0; j<N; j++) PWE += node[i].C_rec.l2(node[j].C_rec);
	PWE /= N*N;
*/

	fappend("data/E.txt", fmt("%f\n", REC));
//	fappend("data/stats/PWE.txt", fmt("%f\n", PWE));

//	Matrix Lmean(D,1);
//	for(int i=0; i<N; i++) Lmean += node[i].L;
//	Lmean.normalize_l1();
//	Lmean.save(fmt("data/Lmean/%d.txt",t));
//	plot_Lmean();
//	plot_mat_plotall();
}

void compute_theoretic_covariance() {
	FULL_COVARIANCE_th.clear();
	Matrix mumut(D,D);
	Matrix mu(D,1);
	for(int i=0;i<N; i++) {
		FULL_COVARIANCE_th.sadd((float)node[i].n/n, node[i].C);
		mu.sadd((float)node[i].n/n, node[i].mu);
	}
	mumut.AAt(mu);
	FULL_COVARIANCE_th -= mumut;

	FULL_COVARIANCE_th.save("data/cov_th.txt");
	DBG("ERROR TH : " << FULL_COVARIANCE.l2(FULL_COVARIANCE_th));
}


Matrix Ufull;
Matrix Lfull;
__multithread__(dump_latepca_errors) (int nbprojs) {
	Matrix C(D,D);
	Matrix L(D,1); L = Lfull; for(int i=nbprojs; i<D; i++) L[i]=0;
	C.reconstruct(Ufull,L);
	float err = FULL_COVARIANCE.l2(C) / FULL_COVARIANCE.n2();
	CRITICAL_BEGIN();
		fappend("data/stats/latepca_errors.txt", fmt("%d %f\n", nbprojs, err));
	CRITICAL_END();
}

int main(int argc, char **argv) {
	try {
	DBG_START("Init ");
	if(argc<=1) {
		init("/home/jfellus/Documents/These/prog/bases/mnist/train.idx3-ubyte");
	} else {
		DBG("with : " << argv[1]);
		init(argv[1]);
		//chdir(dirname(argv[1]));
	}

//	if(system("mkdir -p plots/mat_plotall")) {}
//	if(system("mkdir -p data/mat_plotall")) {}
//	if(system("mkdir -p data/x")) {}
//	if(system("mkdir -p data/ul")) {}
//	if(system("mkdir -p data/cov")) {}
//	if(system("mkdir -p data/stats")) {}
//	if(system("mkdir -p data/Lmean")) {}

	DBG_END();

	DBG_START("Compute global covariance");
	// Centering X and compute full covariance
	X.mean_row(global_mu);
//	global_mu.save("data/gmu.txt");
//	X.save("data/x.txt");
	X -= global_mu;
//	X.save("data/xc.txt");

	FULL_COVARIANCE.covariance(X);
//	plot(FULL_COVARIANCE, "data/cov.jpg", -0.01, 0.01);
//	FULL_COVARIANCE.save("data/cov_full.txt");

	DBG_END();

	////////////////////////

	for(int i=0; i<N; i++) {DBG_PERCENT(i/N);node[i].local_pca();}
	compute_errors();

	compute_theoretic_covariance();

	for(int i=0; i<N; i++) node[i].init_gossip();

	for(t=0; t<T_MAX; t++) {
		DBGV(t);
		gossip();
		compute_errors();
	}

	//////////////////////


//	Ufull.create(D,D); Lfull.create(D,1);
//	FULL_COVARIANCE.eigendecompose(Ufull, Lfull);
//	dump_latepca_errors(D);

	DBG("finished");
	}catch(const char* c) {DBG("ERROR : " << c );}
}
