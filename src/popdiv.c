/* Implementation of Dirichlet-Multinomial model for assessing population  */
/* differentiation and isolation from unlinked genetic marker data.  */
/* For details see Nicholson et al. (2002), JRSS(B), 64, 1-21  */
/* and the discussion by J. L. Marchini and L. R. Cardon. */

/* code written by Jonathan Marchini (Department of Statistics, Oxford, 2002) */


#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>


double log_Lik_bottleneck_constant();
double log_Lik_bottleneck(double *p, double *c);
double log_Lik_bottleneck_p(double *p, double *c, int popn);
double log_Lik_bottleneck_l(double *p, double *c, int locus);
void update_p(double *p, double *c, double *new_p);
void jump_dirichlet(double *vec, int *n, double mult, double *log_ratio_dens);
void update_c(double *p, double *c, double *new_c);
void jump_normal(double *mu, double *sd, double *s, double *ratio_dens, double low_lim);
void rdirichlet(double *ans, double *probs, int *l);
double ddirichlet(double *x, double *pars, int *l);

void sequential_mean(double*, double*, double);
void sequential_variance(double*, double*, double*, double);
void sequential_sd(double*, double*, double*, double);
void sequential_mean_vec(double*, double*, double, int);
void sequential_variance_vec(double*, double*, double*, double, int);
void sequential_sd_vec(double*, double*, double*, double, int);


double *temp1, *temp2, *temp3, *d, *tmp_p, *tmp_c;
int C_ACCEPT_COUNT, P_ACCEPT_COUNT;
int *NUMA, L, P, iter, burnin, max_alleles, s = 1;
double *X, *N, p_mult, c_sd;



void popdiv(int *LL, int *PP, int *NUMAA, double *NN, double *XX, int *Burnin, int *Iter, int *Max_alleles,
	    double *pmult, double *csd, int *outc, double *c_out, 
	    double *mu_c, double *sd_c, double *mu_p, double *sd_p, 
	    double *p_acc_rate, double *c_acc_rate) {
	
	
	
	int i, j, l, it;
	double **p1, **c1, log_lik_const;
	
  
  /*=========*/
  /* set RNG */
  /*=========*/

  GetRNGstate();


  /*======================*/
  /* Initialize variables */
  /*======================*/

  L = *LL;
  P = *PP;
  burnin = *Burnin;
  iter = *Iter;
  max_alleles = *Max_alleles;
  p_mult = *pmult;
  c_sd = *csd;
  P_ACCEPT_COUNT = 0;
  C_ACCEPT_COUNT = 0;

  /*===================*/
  /* set storage space */
  /*===================*/

  NUMA = (int *) R_alloc(L, sizeof(int));
  N = (double *) R_alloc(L * P, sizeof(double));
  X = (double *) R_alloc(L * P * max_alleles, sizeof(double));

  for(l = 0; l < L; l++) NUMA[l] = NUMAA[l];
  for(l = 0; l < L * P; l++) N[l] = NN[l];
  for(l = 0; l < L * P * max_alleles; l++) X[l] = XX[l];
   
  d = (double *) R_alloc(P, sizeof(double));

  temp1 = (double *) R_alloc(max_alleles, sizeof(double));
  temp2 = (double *) R_alloc(max_alleles, sizeof(double));
  temp3 = (double *) R_alloc(max_alleles, sizeof(double));

  tmp_p = (double *) R_alloc(L * max_alleles, sizeof(double));
  tmp_c = (double *) R_alloc(P, sizeof(double));
  
  p1 = (double **) R_alloc(2, sizeof(double *));
  c1 = (double **) R_alloc(2, sizeof(double *));
  for (i = 0; i < 2; i++) {
	  p1[i] = (double *) R_alloc(L * max_alleles, sizeof(double));	
	  c1[i] = (double *) R_alloc(P, sizeof(double));
  }
  
  for(i = 0; i < L; i++) 
	  for(j = 0; j < NUMA[i]; j++) 
		  p1[0][i * max_alleles + j] = 1.0 / (double) NUMA[i];
  
  for(j = 0; j < P; j++) 
	  c1[0][j] = 0.01;
  
  /*=========*/
  /* Burn-in */
  /*=========*/
  
  Rprintf("Burnin...");
  for(it = 0; it < burnin / 2; it++) {
	  
	  /* updates */	
	  update_p(p1[0], c1[0], p1[1]);
	  update_c(p1[1], c1[0], c1[1]);
	  update_p(p1[1], c1[1], p1[0]);
	  update_c(p1[0], c1[1], c1[0]);		
	  
  }  
  Rprintf("finished\n");

  /*=================*/
  /* Main iterations */
  /*=================*/
  
  log_lik_const = log_Lik_bottleneck_constant();

  Rprintf("Sampling...");
  for(it = 0; it < iter / 2; it++) {
	  

	  /* updates */
	  update_p(p1[0], c1[0], p1[1]);
	  update_c(p1[1], c1[0], c1[1]);
	  update_p(p1[1], c1[1], p1[0]);
	  update_c(p1[0], c1[1], c1[0]);
	  
	  /* calculate mean and variance of c sequentially */	
	  
	  sequential_mean_vec(mu_p, p1[1], (double) (2 * it + 1), L * max_alleles);
	  sequential_sd_vec(sd_p, mu_p, p1[1], (double) (2 * it + 1), L * max_alleles);
	  
	  sequential_mean_vec(mu_p, p1[0], (double) (2 * it + 2), L * max_alleles);
	  sequential_sd_vec(sd_p, mu_p, p1[0], (double) (2 * it + 2), L * max_alleles);
	  
	  sequential_mean_vec(mu_c, c1[1], (double) (2 * it + 1), P);
	  sequential_sd_vec(sd_c, mu_c, c1[1], (double) (2 * it + 1), P);
	  
	  sequential_mean_vec(mu_c, c1[0], (double) (2 * it + 2), P);
	  sequential_sd_vec(sd_c, mu_c, c1[0], (double) (2 * it + 2), P);
	  
	  /* write samples to files */
	  

	  if(*outc == 1) {	  

		  for(j = 0; j < P; j++) {
			  *(c_out + 2 * it * P + j) = c1[1][j];
			  *(c_out + 2 * it * P + P + j) = c1[0][j];
		  }
	  }
		
	  
	  
  }
  Rprintf("finished\n");

  *p_acc_rate = (double) P_ACCEPT_COUNT;
  *c_acc_rate = (double) C_ACCEPT_COUNT;

  /*=========*/
  /* set RNG */
  /*=========*/

  PutRNGstate();

  
}

double log_Lik_bottleneck_constant() {
	
	int i, j, k;
	double tot;

	tot = 0.0;
	for(i = 0; i < L; i++) {
		for(j = 0; j < P; j++) {
		  tot += lgammafn(*(N + i * P + j) + 1.0);
		  for(k = 0; k < NUMA[i]; k++) 
		    tot -= lgammafn(*(X + i * P * max_alleles + j * NUMA[i] + k) + 1.0);
		}
	}
	
	return tot;
}
	

double log_Lik_bottleneck(double *p, double *c) {
	
	int i, j, k;
	double tot;

	for(j = 0; j < P; j++) d[j] = (1 - c[j]) / c[j];

	tot = 0.0;
	for(i = 0; i < L; i++) {
		for(j = 0; j < P; j++) {
		  tot += lgammafn(*(d + j));
		  for(k = 0; k < NUMA[i]; k++) 
		    tot += lgammafn(*(X + i * P * max_alleles + j * NUMA[i] + k) + *(p + i * max_alleles + k) * *(d + j));

		  tot -= lgammafn(*(N + i * P + j) + *(d + j));	

		  for(k = 0; k < NUMA[i]; k++) 
		    tot -= lgammafn(*(p + i * max_alleles + k) * *(d + j));	
		}
	}
	

	return tot;
}

double log_Lik_bottleneck_p(double *p, double *c, int popn) {
	
  int i, j, k;
  double tot;
  
  for(j = 0; j < P; j++) d[j] = (1 - c[j]) / c[j];

  tot = 0.0;
  for(i = 0; i < L; i++) { 
    tot += lgammafn(*(d + popn));
    for(k = 0; k < NUMA[i]; k++) 
      tot += lgammafn(*(X + i * P * max_alleles + popn * NUMA[i] + k) + *(p + i * max_alleles + k) * *(d + popn));
    
    tot -= lgammafn(*(N + i * P + popn) + *(d + popn));	
    
    for(k = 0; k < NUMA[i]; k++) 
      tot -= lgammafn(*(p + i * max_alleles + k) * *(d + popn));	
  }

  
  return tot;
}

double log_Lik_bottleneck_l(double *p, double *c, int locus) {
	
  int j, k;
  double tot;
  
  for(j = 0; j < P; j++) d[j] = (1 - c[j]) / c[j];
  
  tot = 0.0;	
  for(j = 0; j < P; j++) {
    tot += lgammafn(*(d + j));
    for(k = 0; k < NUMA[locus]; k++) 
      tot += lgammafn(*(X + locus * P * max_alleles + j * NUMA[locus] + k) + *(p + locus * max_alleles + k) * *(d + j));
    
    tot -= lgammafn(*(N + locus * P + j) + *(d + j));	
    
    for(k = 0; k < NUMA[locus]; k++) 
		    tot -= lgammafn(*(p + locus * max_alleles + k) * *(d + j));	
  }
  
  
  return tot;
}

void update_p(double *p, double *c, double *new_p) {
	
	int i, k;
	double MH_prob, log_trans_ratio;


	for(i = 0; i < L; i++) {

	  for(k = 0; k < NUMA[i]; k++) *(tmp_p + i * max_alleles + k) = *(p + i * max_alleles + k);

		jump_dirichlet((tmp_p + i * max_alleles), (NUMA + i), p_mult, &log_trans_ratio);

		MH_prob = log_Lik_bottleneck_l(tmp_p, c, i) - log_Lik_bottleneck_l(p, c, i);	 
 		MH_prob = exp(MH_prob + log_trans_ratio); 
		
		if(runif(0.0, 1.0) < MH_prob) {
		  for(k = 0; k < NUMA[i]; k++) *(new_p + i * max_alleles + k) = *(tmp_p + i * max_alleles + k);
		  P_ACCEPT_COUNT++;
		}
		else {
		  for(k = 0; k < NUMA[i]; k++) *(new_p + i * max_alleles + k) = *(p + i * max_alleles + k);
		}
	}
}


void jump_dirichlet(double *vec, int *n, double mult, double *log_ratio_dens) {
  
	int i;
	double p1, p2;
 	

	for(i = 0; i < *n; i++) { 
	  temp1[i] = vec[i]; 
	  temp2[i] = mult * vec[i];
	}

	rdirichlet(vec, temp2, n);

	p1 = ddirichlet(vec, temp2, n);

	for(i = 0; i < *n; i++) 
	  temp3[i] = mult * vec[i];
       
	p2 = ddirichlet(temp1, temp3, n);
	

	*log_ratio_dens = p2 - p1;
	
}

void update_c(double *p, double *c, double *new_c) {
	
	int j;
	double MH_prob, trans_ratio;


	for(j = 0; j < P; j++) {

		jump_normal((c + j), &c_sd, (tmp_c + j), &trans_ratio, 0.0);
		
		MH_prob = log_Lik_bottleneck_p(p, tmp_c, j) - log_Lik_bottleneck_p(p, c, j);
		MH_prob = exp(MH_prob) * trans_ratio;
		
		if(runif(0.0, 1.0) < MH_prob) {
		  new_c[j] = tmp_c[j];		  
		  C_ACCEPT_COUNT++;
		}
		else {
		  new_c[j] = c[j];		  
		}	  
	}
}


void jump_normal(double *mu, double *sd, double *s, double *ratio_dens, double low_lim) {

	double p1, p2;
 	

	do 
		*s = rnorm(*mu, *sd);
	while ((*s > 1.0) || (*s < low_lim));

	p1 = (pnorm(1.0, *mu, *sd, 1L, 0L) - pnorm(0.0, *mu, *sd, 1L, 0L));	
	p2 = (pnorm(1.0, *s, *sd, 1L, 0L) - pnorm(0.0, *s, *sd, 1L, 0L));
	p1 *= dnorm(*s, *mu, *sd, 0L);	
	p2 *= dnorm(*mu, *s, *sd, 0L);

	*ratio_dens = p2 / p1;
	
}


void rdirichlet(double *ans, double *probs, int *l){
	
	/* generate a single draw from a Dirichlet(probs) distribution */

	int i;
	double t = 0.0;
	for(i = 0; i < *l; i++) {
		*(ans + i) = rgamma(*(probs + i), 1.0);
		t += *(ans + i);
	}
	for(i = 0; i < *l; i++) *(ans + i) /= t;
}

double ddirichlet(double *x, double *pars, int *l) {

	/* returns the log-density of x from a Dirichlet distribution with parameter vector pars */
	int i;
	double t1 = 0.0, log_dens = 0.0;
		
	for(i = 0; i < *l; i++) {
		t1 += *(pars + i);
		log_dens += (*(pars + i) - 1) * log(*(x + i)) - lgammafn(*(pars + i));
	}
	
	log_dens += lgammafn(t1);
			
	return(log_dens);
}

void sequential_mean(double *mu, double *x, double m) {

	/* calculate a mean sequentially */

	*mu = (*x + (m - 1) * *mu) / m;

}

void sequential_variance(double *var, double *mu, double *x, double m) {

	/* calculate the variance sequentially */

	if(m == 1) *var = *x;
	if(m == 2) *var = ((*var - *x) / 2.0) * (*var - *x);
	if(m > 2) *var = (m - 2.0) * *var / (m - 1.0) + m * (*x - *mu) * (*x - *mu) / ((m - 1.0) * (m - 1.0)); 
}

void sequential_sd(double *sd, double *mu, double *x, double m) {

	/* calculate the variance sequentially */

	if(m == 1) *sd = *x;
	if(m == 2) *sd = sqrt(((*sd * *sd - *x) / 2.0) * (*sd * *sd - *x));
	if(m > 2) *sd = sqrt((m - 2.0) * *sd * *sd / (m - 1.0) + m * (*x - *mu) * (*x - *mu) / ((m - 1.0) * (m - 1.0))); 
}


void sequential_mean_vec(double *mu, double *x, double m, int vec_length) {

	int i;

	for(i = 0; i < vec_length; i++) sequential_mean(mu + i, x + i, m);
}

void sequential_variance_vec(double *var, double *mu, double *x, double m, int vec_length) {

	int i;

	for(i = 0; i < vec_length; i++) sequential_variance(var + i, mu + i, x + i, m); 
}

void sequential_sd_vec(double *sd, double *mu, double *x, double m, int vec_length) {

	int i;

	for(i = 0; i < vec_length; i++) sequential_sd(sd + i, mu + i, x + i, m); 
}
