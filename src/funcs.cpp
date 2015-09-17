/*
 * funcs.cpp
 *
 *  Created on: 16 sept. 2015
 *      Author: jfellus
 */



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
