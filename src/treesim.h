 
/* Standard neutral coalescent simulator */

/* code written by Jonathan Marchini (Department of Statistics, Oxford, 2003) */

using namespace std;

extern "C" { 
  void samp(int, int, int*);
}

int rmultinom1(vector<double> weights);

struct node {
	
	double branch_length;
	vector<double> mutations;
	vector<int> as;
	int upper;
	int lower;
	double recomb_rate;
	
	int status; /* (1 = coalesces, 2 = recombinant) */
	node *parent_left;
	node *parent_right;
	node *left;
	node *right;

};

class ctree
{
 public:
  int n, Sn;
  double theta, tmp1;

  node *MRCA;
  vector<node*> leaves;
  vector<node*> cl;
  vector<node*> node_list;

  vector<double> recomb_vec;
  double rate_recomb;
  double rate_coalesce;

  vector<int> tas;
  vector<double> weights;

  set<double> mutations;
  vector<vector<int> > sample;

  double tree_time;
  double total_time;
  double seg_sites;

  ctree();
  ctree(int n, double theta, int Sn1, double *recomb_vec1);
  ~ctree();
  void destroy_tree();
  void generate_topology();
  void coalesce();
  void coalesce(int allele1, int allele2, int n);
  void recombine();
  void put_mutations();
  void put_mutation(node* allele);
  void make_sample();
  void get_mutations(node* leaf, int lower, int upper, vector<int>* vec, vector<int>* ans);
  void gen_tree();
};


ctree::ctree()
{
  MRCA = NULL;
}


ctree::~ctree()
{
	destroy_tree();
}

void ctree::destroy_tree() {
	
	int i;
	for(i = 0; i < (int) node_list.size(); i++) 
		delete node_list[i];
	
	sample.clear();
	mutations.clear();
	leaves.clear();
	cl.clear();
	node_list.clear();
	recomb_vec.clear();
	weights.clear();	
}


ctree::ctree(int n1, double theta1, int Sn1, double *recomb_vec1) {

	int i;
	n = n1;
	theta = theta1;
	Sn = Sn1;
	MRCA = NULL;
	tree_time = 0.0;
	total_time = 0.0;
	
	for(i = 0; i < Sn; i++) {
		recomb_vec.push_back(recomb_vec1[i]);
	}

	rate_recomb = n * recomb_vec[Sn - 1] / 2;
	rate_coalesce = exp(lchoose((double) n, 2.0));
	
	tas.assign(Sn, n);
	
	for(i = 0; i < n; i++) {
		
		leaves.push_back(new node);
		(leaves[i])->right = NULL;
		(leaves[i])->left = NULL;
		(leaves[i])->parent_right = NULL;
		(leaves[i])->parent_left = NULL;
		(leaves[i])->branch_length = 0.0;
		(leaves[i])->lower = 1;
		(leaves[i])->upper = Sn;
		(leaves[i])->recomb_rate = recomb_vec[Sn - 1] / 2;
		(leaves[i])->as.assign(Sn , 1);
		
		cl.push_back(leaves[i]);
		weights.push_back((leaves[i])->recomb_rate);
		node_list.push_back(leaves[i]);		
	}


}

void ctree::generate_topology() {

	while (cl.size() > 1) {
	
		if(runif(0.0, 1.0) < rate_coalesce / (rate_coalesce + rate_recomb))
			coalesce();
		else
			recombine();
		
	}
	MRCA = cl[0];	
}
	
void ctree::coalesce() {
	
	int vec[2], tmp;
	samp(cl.size(), 2, vec);
	if(vec[1] < vec[0]) {
		tmp = vec[0];
		vec[0] = vec[1];
		vec[1] = tmp;
	}
	
	coalesce(vec[0] - 1, vec[1] - 1, cl.size());
}

void ctree::coalesce(int edge1, int edge2, int n) {

	float av = 1.0 / (float) (rate_coalesce + rate_recomb);
	double time = rexp(av);
	int flag = 1, i, j;
	
	//Rprintf("coalesce\n");

	for(i = 0; i <  (int) cl.size(); ++i) 
		(cl[i])->branch_length += time;
	
	tree_time += (n * time);
	total_time += time;
	
	(cl[edge1])->status = 1;
	(cl[edge2])->status = 1;
	
	(cl[edge1])->parent_left = new node; 
	(cl[edge2])->parent_left = (cl[edge1])->parent_left;
	
	(cl[edge1])->parent_right = NULL;
	(cl[edge2])->parent_right = NULL;
	
	(cl[edge1])->parent_left->left = cl[edge1];
	(cl[edge1])->parent_left->right = cl[edge2];
	(cl[edge1])->parent_left->branch_length = 0.0;
	(cl[edge1])->parent_left->as.resize(Sn);
	(cl[edge1])->parent_left->parent_left = NULL;
	(cl[edge1])->parent_left->parent_right = NULL;

	(cl[edge1])->parent_left->upper = 0;
	(cl[edge2])->parent_left->lower = 0;

	
	for(i = 0; i < Sn; i++) {
		
		if( ((cl[edge1])->as[i] == 1) || ((cl[edge2])->as[i] == 1) ) {
			(cl[edge1])->parent_left->as[i] = 1;  
			
			if(flag == 1) {
				(cl[edge1])->parent_left->lower = i + 1;
				flag = 0;
			}
			(cl[edge1])->parent_left->upper = i + 1;
		}
		else
			(cl[edge1])->parent_left->as[i] = 0;
		
		if( ((cl[edge1])->as[i] == 1) && ((cl[edge2])->as[i] == 1) ) 
			tas[i] -= 1;
	}
	cl.push_back((cl[edge1])->parent_left);
	node_list.push_back((cl[edge1])->parent_left);
	cl.erase(cl.begin() + edge2);
	cl.erase(cl.begin() + edge1);
	weights.resize(cl.size());
	
	for(j = 0; j <  (int) cl.size(); ++j) {
		flag = 1;
		
		(cl[j])->lower = 1;
		(cl[j])->upper = 1;
		
		for(i = 0; i < Sn; i++) {
			if(tas[i] == 1) 
				(cl[j])->as[i] = 0;  
			
			
			if((cl[j])->as[i] == 1) {
				if(flag == 1) {
					(cl[j])->lower = i + 1;
					flag = 0;  
				}
				(cl[j])->upper = i + 1;
			}
		}  
		(cl[j])->recomb_rate = (recomb_vec[(cl[j])->upper - 1] - recomb_vec[(cl[j])->lower - 1]) / 2;
		weights[j] = (cl[j])->recomb_rate;
		
	}
	rate_recomb = accumulate(weights.begin(), weights.end(), 0.0);
	rate_coalesce = exp(lchoose((double) cl.size(), 2.0));
}


void ctree::recombine() {
  
	int i, edge, split, flag1, flag2;
	double r, tot;

	//Rprintf("recombine\n");


	edge = rmultinom1(weights);
	
	r = runif(0.0, 1.0);
	tot = recomb_vec[(cl[edge - 1])->upper - 1] - recomb_vec[(cl[edge - 1])->lower - 1];
	       
	split = (cl[edge - 1])->lower;
	while(r > ((recomb_vec[split] - recomb_vec[(cl[edge - 1])->lower - 1]) / tot)) {
		split++;
	}
	float av = 1.0 / (float) (rate_coalesce + rate_recomb);
	double time = rexp(av);
	
	for(i = 0; i < (int) cl.size(); ++i) 
		(cl[i])->branch_length += time;
	
	tree_time += (cl.size() * time);
	total_time += time;

	(cl[edge - 1])->status = 2;

	(cl[edge - 1])->parent_left = new node; 
	(cl[edge - 1])->parent_right = new node; 
  
	(cl[edge - 1])->parent_right->left = cl[edge - 1];
	(cl[edge - 1])->parent_right->right = NULL;
	(cl[edge - 1])->parent_right->parent_right = NULL;
	(cl[edge - 1])->parent_right->parent_left = NULL;
	(cl[edge - 1])->parent_left->right = cl[edge - 1];
	(cl[edge - 1])->parent_left->left = NULL;
	(cl[edge - 1])->parent_left->parent_left = NULL;
	(cl[edge - 1])->parent_left->parent_right = NULL;
	(cl[edge - 1])->parent_left->branch_length = 0.0;
	(cl[edge - 1])->parent_right->branch_length = 0.0;
	(cl[edge - 1])->parent_left->as.resize(Sn);
	(cl[edge - 1])->parent_right->as.resize(Sn);
	
	(cl[edge - 1])->parent_right->upper = 0;
	(cl[edge - 1])->parent_right->lower = 0;
	(cl[edge - 1])->parent_left->upper = 0;
	(cl[edge - 1])->parent_left->lower = 0;

	flag1 = 1;
	flag2 = 1;
	for(i = 0; i < Sn; i++) {
		if( (i + 1) <= split ) {			
			if( (cl[edge - 1])->as[i] == 1 )
				(cl[edge - 1])->parent_left->as[i] = 1;
			else
				(cl[edge - 1])->parent_left->as[i] = 0;	
			if((cl[edge - 1])->parent_left->as[i] == 1) {

				if(flag1 == 1) {
					(cl[edge - 1])->parent_left->lower = i + 1;
					flag1 = 0;
				}
				(cl[edge - 1])->parent_left->upper = i + 1;
			}					
		}
		else {
			if( (cl[edge - 1])->as[i] == 1 )
				(cl[edge - 1])->parent_right->as[i] = 1;
			else
				(cl[edge - 1])->parent_right->as[i] = 0;
			if((cl[edge - 1])->parent_right->as[i] == 1) {

				if(flag2 == 1) {
					(cl[edge - 1])->parent_right->lower = i + 1;
					flag2 = 0;
				}
				(cl[edge - 1])->parent_right->upper = i + 1;
			}				
		}
	}


	(cl[edge - 1])->parent_left->recomb_rate = (recomb_vec[(cl[edge - 1])->parent_left->upper - 1] - recomb_vec[(cl[edge - 1])->parent_left->lower - 1]) / 2;	
	(cl[edge - 1])->parent_right->recomb_rate = (recomb_vec[(cl[edge - 1])->parent_right->upper - 1] - recomb_vec[(cl[edge - 1])->parent_right->lower - 1]) / 2;

	cl.push_back((cl[edge - 1])->parent_left);
	cl.push_back((cl[edge - 1])->parent_right);
	node_list.push_back((cl[edge - 1])->parent_left);
	node_list.push_back((cl[edge - 1])->parent_right);
	cl.erase(cl.begin() + edge - 1);
	  
	weights.resize(cl.size());

	rate_recomb = 0.0;
	for(i = 0; i < (int) cl.size(); ++i) {
	  weights[i] = (cl[i])->recomb_rate;
	  rate_recomb += (cl[i])->recomb_rate;
	}

	rate_coalesce = exp(lchoose((double) cl.size(), 2.0));
}



void ctree::put_mutations() {
	
	int i, j, num, interv;
	float mean = 0.0;
	double mut;	
	
	for(i = 0; i < (int) (node_list.size() - 1); i++) {
		mean = (theta / 2) * ((node_list[i])->branch_length);
		num = (int) rpois(mean);
		
		for(j = 0; j < num; j++) {
			mut = runif(0.0, 1.0);
			interv = (int) floor(Sn * mut);
			if((node_list[i])->as[interv] == 1) {
				((node_list[i])->mutations).push_back(mut);
				mutations.insert(mut);
			}
		}
				
	}

			
	seg_sites = (double) mutations.size();

}
		

void ctree::make_sample() {

  int i;	
  vector<int> vec((int) seg_sites, 0);
  vector<int> ans(Sn, 1);

  for(i = 0; i < (int) leaves.size(); ++i) {
	  get_mutations(leaves[i], 1, Sn, &vec, &ans);
	  sample.push_back(vec);
	  vec.clear();
	  vec.resize((int) seg_sites);
	  ans.clear();
	  ans.resize(Sn, 1);
  }
}
  


void ctree::get_mutations(node* leaf, int lower, int upper, vector<int>* vec, vector<int>* ans) {
	
	node *leaf1 = leaf;
	int j, pos, interv, sum;
	vector<int> ans_copy(Sn);
	
	sum = 0;
	for(j = 0; j < Sn; j++) {
		(*ans)[j] *= leaf1->as[j];
		ans_copy[j] = (*ans)[j];
		sum += (*ans)[j];
	}
	if(sum != 0) {
		
		for(j = 0; j < (int) (leaf1->mutations).size(); ++j) {
			interv = (int) floor(Sn * (leaf1->mutations)[j]);
			if((*ans)[interv] == 1) {
				pos = distance(mutations.begin(), mutations.find((leaf1->mutations)[j]));
				(*vec)[pos] = 1;
			}
		}
		if(leaf1->status == 1) {
			if(leaf1->parent_left != MRCA) {
				
				get_mutations(leaf1->parent_left, lower, upper, vec, ans);
			}
		}
		if(leaf1->status == 2) {	
			get_mutations(leaf1->parent_left, max(lower, leaf1->lower), min(upper, leaf1->upper), vec, ans);
			get_mutations(leaf1->parent_right, max(lower, leaf1->lower), min(upper, leaf1->upper), vec, &ans_copy);
		}
	}	
}

void ctree::gen_tree() {
	
  generate_topology();
  put_mutations();
  make_sample();

}
	

  
int rmultinom1(vector<double> weights) {
	
	int i, j, n = weights.size();
	float sum = 0.0, prob;
	
	
	for(i = 0; i < (n - 1); i++)
		sum += weights[i];
	
	for(i = 0; i < (n - 1); i++) {
		prob = (float) weights[i] / sum;
		j = (int) rbinom(1.0, prob);
		if(j == 1) return (i + 1);
		sum -= weights[i];
	}
	
	return n;
	
}


int rmultinom_unif(int n) {
	
	int i, j;
	double sum = 1.0, *probs, prob;
	
	probs = (double *) R_alloc(n, sizeof(double));
	
	for(i = 0; i < n; i++) 
		probs[i] = 1.0 / ((double) n);
	
	for(i = 0; i < (n - 1); i++) {
		prob = probs[i] / sum;
		j = (int) rbinom(1.0, prob);
		if(j == 1) return (i + 1);
		sum -= probs[i];
	}
	
	return n;
	
}

  
  
