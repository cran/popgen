
/* Standard neutral coalescent simulator */

/* code written by Jonathan Marchini (Department of Statistics, Oxford, 2003) */

#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector> 
#include <set> 
#include <algorithm>
#include <numeric>
#include <cmath> 
#include <math.h>
#include "treesim.h"


extern "C" {
  
  using namespace std;
  
  void samp(int N, int p, int *vec);
  
  
  void treesim(int *n, int *Sn, double *theta, double *rho_vec,
	       int *num_seg_sites, double *tree_time, double *total_time, 
	       int *mut_flag, double *mutations, int *samp_flag, int *sample){
    
    int i, j;
    
    /////////////
    // SET RNG //
    /////////////
    
    GetRNGstate();
    
    //===============//
    // Generate tree //
    //===============//
    
    ctree tree(*n, *theta, *Sn, rho_vec);
    tree.gen_tree();
    
    //===============//
    // Store output  //
    //===============//
    
    *num_seg_sites = (int) tree.seg_sites;
    *tree_time = tree.tree_time;
    *total_time = tree.total_time;
    
    if(*mut_flag == 1) {
      set<double>::iterator pos;
      j = 0;
      for(pos = tree.mutations.begin(); pos != tree.mutations.end(); ++pos) { 
	mutations[j] = *pos;
	j++;
      }
    }
    
    if(*samp_flag == 1) {
      for(i = 0; i < *n; i++) {
	for(j = 0; j < *num_seg_sites; j++) {
	  *(sample + i * *num_seg_sites + j) = tree.sample[i][j];
	}
      }
    }
    
    ////////////////
    // UN-SET RNG //
    ////////////////
    
    PutRNGstate();
    
  }
  
  
  void samp(int N, int p, int *vec) {
    
    /* random sampling without replacement from 1:N */
    
    int i, j, tmp, *pool;
    
    pool = (int *) R_alloc(N, sizeof(int));
    for(i = 0; i < N; i++) pool[i] = i + 1;
    
    for(i = 0; i < p; i++) {
      tmp = rmultinom_unif(N - i);
      vec[i] = pool[tmp - 1];
      for(j = (tmp - 1); j < (N - i - 1); j++) pool[j] = pool[j + 1];
    }
    
  }
}

