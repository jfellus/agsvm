/*
 * misc.cpp
 *
 *  Created on: 13 nov. 2013
 *      Author: jfellus
 */




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




void init(const char* datafile) {
	DBG("INIT");
	DBGV(NBTHREADS);
	system("rm -rf data/*");
	system("rm -rf plots/*");

	X.load(datafile);
	if(LIMIT_NDATA!=-1 && X.height > LIMIT_NDATA) X.height = LIMIT_NDATA;
	n = X.height;
	D = X.width;

	global_mu.create(D,1);

	FULL_COVARIANCE.create(D,D);
	FULL_COVARIANCE_th.create(D,D);


	node = new Node[N];
	create_network();

	int ndo = 0;
	for(int i=0; i<N-1; i++) {
		node[i].init(X, ndo, n/N);
		ndo += n/N;
	}
	node[N-1].init(X,ndo, n-ndo);

//	mat_plotall.create((D+1)*(int)(sqrt(N)+1),(D+1)*(int)(sqrt(N)+1));

	DBGV(LIMIT_NDATA);
	DBGV(N);
	DBGV(D);
	DBGV(n);
	if(USE_ENERGY) {DBGV(KEEP_ENERGY);}
	else {DBGV(q);}
}




void deinit() {
	delete[] node;
}
