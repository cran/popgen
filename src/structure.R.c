
#include <R.h>
#include <Rmath.h>
#define min(a,b) ((a) <= (b) ? (a) : (b))

void logLik_calc(double*, int*, double*, double*, int*, int*, int*, int*);
void update_P(double*, double*, int*, int*, double*, double*, int*, int*, int*, int*, int*, int*, double*, double*);
void update_Q(double*, double*, int*, int*, int*, int*, double*); 
void update_Alpha(double*, double*, double*, int*, int*, double*); 		
void update_Z(int*, int*, double*, double*, int*, int*, int*, int*, double*, int*);
void update_pi(double*, double *, double*, double*, double*, double*, double*, int*, int*, int*, int*, double*, double*, double*);
void update_cc(double*, double *, double*, double*, double*, double*, int*, int*, int*, int*, double*, double*, double*, double*, double*);
void jump_pi(double*, double*, int*, double*, int*, int*, int*, double*, double*);
void jump_c(double*, double*, int*, double*, int*, double*, double);
void rmultinom_int_JM(int*, double*);
void rdirichlet_JM(double*, double*, int*);
double ddirichlet_JM(double*, double*, int*);

void structure_JM(int *X, 
		  int *J, 
		  int *maxJ, 
		  int *n, 
		  int *L, 
		  int *K, 
		  int *burn_in, 
		  int *num_sample, 
		  int *thin, 
		  int *print_step, 
		  int *sample_type, 
		  double *alpha_start, 
		  double *alpha_sd, 
		  int *alpha_step, 
		  double *max_alpha, 
		  int *fix_alpha, 
		  double *beta, 
		  int *use_pi_and_c, 
		  int *z_init, 
		  int *z_init_mat, 
		  int *fix_z, 
		  double *pi_init_mat, 
		  int *fix_pi, 
		  double *c_init_mat, 
		  double *c_prior_mu, 
		  double *c_prior_sd, 
		  double *c_sd, 
		  double *alpha_sample,  
		  double *q_sample,  
		  double *p_sample,
		  int *z_sample, 
		  double *pi_sample, 
		  double *c_sample, 
		  double *logLik_sample) {
	
	int i, j, k, l, it, num, store0, *Z, samp;
	double *P, *N, *Q, *M, *logLik, Alpha;
	double *pi, *pi_temp, *c, *c_temp;
	double *store1, *store2, *store3, *store4, *store5, *store6;
	
	/* set random number generator */
	GetRNGstate(); 
	
	/* allocate memory for small workspace vectors */
	store1 = (double *) R_alloc(*maxJ, sizeof(double));
	store2 = (double *) R_alloc(*maxJ, sizeof(double));
	store3 = (double *) R_alloc(*maxJ, sizeof(double));
	store4 = (double *) R_alloc(*maxJ, sizeof(double));
	store5 = (double *) R_alloc(*maxJ, sizeof(double));
	store6 =  (double *) R_alloc(*K, sizeof( double));
	
	/* allocate memory for main workspace arrays */
	P =  (double *) R_alloc (((*maxJ) * (*K) * (*L)), sizeof(double));
	N =  (double *) R_alloc (((*maxJ) * (*K) * (*L)), sizeof( double));
	Q = (double *) R_alloc (((*n) * (*K)), sizeof(double));
	M = (double *) R_alloc (((*n) * (*K)), sizeof(double));
	Z = (int *) R_alloc ((2 * (*n) * (*L)), sizeof(int));
	pi = (double *) R_alloc (((*maxJ) * (*L)), sizeof(double));
	c = (double *) R_alloc ((*K), sizeof(double));
	
	
	/* allocate memory for storage arrays */
	pi_temp = (double *) R_alloc (((*maxJ) * (*L)), sizeof(double));
	c_temp = (double *) R_alloc ((*K), sizeof(double));
	logLik = (double *) R_alloc (*num_sample, sizeof(double));
	for(i = 0; i < *num_sample; i++) *(logLik + i) = 0.0;
	
	/* initialization of Z */
	if(z_init == 0) {
		for(k = 0; k < *K; k++) *(store6 + k) = 1.0 / ((double) (*K)); 
		for(i = 0; i < (2 * (*n) * (*L)); i++) rmultinom_int_JM((Z + i), store6);
	}
	else {
		for(i = 0; i < (2 * (*n) * (*L)); i++) *(Z + i) = *(z_init_mat + i);
	}	
	
	/* initialisaztion of pi */
	for(l = 0; l < *L; l++) {
		for(j = 0; j < *maxJ; j++) *(pi + l * *maxJ + j) = *(pi_init_mat + l * *maxJ + j);
	}
	
	/* initialisaztion of c */
	for(k = 0; k < *K; k++) *(c + k) = *(c_init_mat + k);
	
	/* initialisaztion of Alpha and num */
	Alpha = *alpha_start;
	num = 0;
	
	/* main loop */
	Rprintf("Burn-in ... ");
	Rprintf("[%07d]", 2);
	for(it = -(*burn_in); it <= (*num_sample * *thin) ;it++){
		
		if(it < 0) Rprintf("\b\b\b\b\b\b\b\b\b[%07d]", (it + *burn_in) + 1);
		if(it == 0) Rprintf(" finished\n");	
		if(it == 0) Rprintf("Iterations completed ... [%07d]", 0);
		if(it > 0) Rprintf("\b\b\b\b\b\b\b\b\b[%07d]", it);
		
		/*----------*/
		/* update P */
		/*----------*/
		update_P(P, N, X, Z, pi, c, n, K, L, use_pi_and_c, J, maxJ, store1, store2);
		
		/*----------*/
		/* update Q */
		/*----------*/
		update_Q(Q, M, Z, n, K, L, &Alpha);
		
		/*-------------*/
		/* update Alpha*/
		/*-------------*/
		if((*fix_alpha == 0) && (*alpha_step * (it / *alpha_step) == it)) 
			update_Alpha(&Alpha, alpha_sd, max_alpha, n, K, Q);
		
		/*----------*/
		/* update Z */
		/*----------*/
		if(*fix_z == 0)
			update_Z(Z, X, P, Q, n, L, K, maxJ, store6, &store0);
		
		/*--------------------------*/
		/* calculate Log Likelihood */
		/*--------------------------*/
		if(it == ((num + 1) * (*(thin)))) 
			logLik_calc((logLik + num), X, P, Q, n, L, K, maxJ);
		
		/*-----------*/
		/* update pi */
		/*-----------*/
		if(*use_pi_and_c == 1) {
			if(*fix_pi == 0) {
				update_pi(P, pi, pi_temp, c, c_temp, c_sd, beta, K, L, J, maxJ, store3, store4, store5);
			}
		}
		/*----------*/
		/* update c */
		/*----------*/
		if(*use_pi_and_c == 1) 
			update_cc(P, pi, pi_temp, c, c_temp, c_sd, K, L, J, maxJ, store3, store4, store5, c_prior_mu, c_prior_sd);
		
		/*-------------------*/
		/* sample parameters */
		/*-------------------*/
		samp =  (num + 1) * *thin;
		if(samp == it) {
			*(alpha_sample + num) = Alpha;	
			if(*sample_type == 1) 
				for(i = 0; i < ((*maxJ) * (*K) * (*L)); i++) *(p_sample + i) += *(P + i);			
			if(*(sample_type + 1) == 1) 
				for(i = 0; i < ((*n) * (*K)); i++) *(q_sample + i) += *(Q + i);
			if(*(sample_type + 2) == 1)
				for(i = 0; i < (2 * (*n) * (*L)); i++) *(z_sample + i) += *(Z + i);
			if(*use_pi_and_c == 1) {
				if(*(sample_type + 3) == 1)  
					for(i = 0; i < (*L) * (*maxJ); i++) *(pi_sample + num * (*L) * (*maxJ) + i) = *(pi + i);
				if(*(sample_type + 4) == 1)
					for(i = 0; i < (*K); i++) *(c_sample + num * (*K) + i) = *(c + i);
			}
			num += 1;
		}
	}
	
	Rprintf("\n");
	for(i = 0; i < *num_sample; i++) *(logLik_sample + i) = *(logLik + i);
	
	/* close down random number generator */
	PutRNGstate();
		
}

void update_P(double *P, double *N, int *X, int *Z, double *pi, double *c, int *n, int *K, int *L, int *type, int *J, int *maxJ, double *store1, double *store2) {
	
	int i, j, k, l, a;
	double t;
	
	for(k = 0; k < *K; k++) { 
		for(l = 0; l < *L; l++) {
					
			/* Calculate N */
			for(j = 0; j < *(J + l); j++) {
				
				t = 0.0;
				for(i = 0; i < *n; i++) { 	
					for(a = 0; a < 2; a++) {
						if( ( *(Z + i*2*(*L) + a*(*L) + l) == (k + 1)) && 
						    ( *(X + i*2*(*L) + a*(*L) + l) == (j + 1))) t += 1.0;
					}
				}
				if(*type == 1) {
					*(N + k*(*maxJ)*(*L) + j*(*L) + l) = t + *(pi + l * *maxJ + j) * (1.0 / *(c + k) - 1.0);
				}
				else {
					*(N + k*(*maxJ)*(*L) + j*(*L) + l) = t + 1.0;
				} 
			}
			
			for(j = 0; j < *(J + l); j++) *(store1 + j) = *(N + k*(*maxJ)*(*L) + j*(*L) + l);
			
			rdirichlet_JM(store2, store1, (J + l));
			
			for(j = 0; j < *(J + l); j++) *(P + k*(*maxJ)*(*L) + j*(*L) + l) = *(store2 + j); 
		}	
	}
}      

void update_Q(double *Q, double *M, int *Z, int *n, int *K, int *L, double *Alpha) {
	
	int i, k, l, a;
	double t;
	
	for(i = 0; i < *n; i++) {
		for(k = 0; k < *K; k++) { 
			
			/* calculate M */
			t = 0.0;	
			for(a = 0; a < 2; a++) { 	
				for(l = 0; l < *L; l++) {				
					if((*(Z + i*2*(*L) + a*(*L) + l)) == (k + 1)) t += 1.0;
				}
			}
			*(M + i*(*K) + k) = t + *Alpha;	
		}
		
		rdirichlet_JM((Q + i*(*K)), (M + i*(*K)), K);
	}
}

void update_Alpha(double *Alpha, double *alpha_sd, double *max_alpha, int *n, int *K, double *Q) {
	
	int i, k;
	double Alpha_new, q, tot1, p1, p2, tmp1, MH_accept_prob;
	
	Alpha_new = rnorm(*Alpha, *alpha_sd);
	if((Alpha_new > 0.01) && (Alpha_new < *max_alpha)) {
		
		q = 0.0;
		tot1 = 1.0;
		for(i = 0; i < *n; i++) {
			for(k = 0; k < *K; k++) {
				
				tot1 *= *(Q + i*(*K) + k);
				
				if (tot1 < 1e-200) {
					q += log(tot1);  
					tot1 = 1.0;
				}
			}
		}
		q += log(tot1);
		
		p1 = (Alpha_new - 1.0) * q + (lgammafn(*K * Alpha_new) - *K * lgammafn(Alpha_new)) * *n; 
		p2 = (*Alpha - 1.0) * q + (lgammafn(*K * *Alpha) - *K * lgammafn(*Alpha)) * *n; 
		
		tmp1 = p1 - p2;	
		
		MH_accept_prob = min(1.0, exp(tmp1));
		if(runif(0.0, 1.0) < MH_accept_prob) *Alpha = Alpha_new;
	}
}


void update_Z(int *Z, int *X, double *P, double *Q, int *n, int *L, int *K, int *maxJ, double *store6, int *store0) {
	
	int i, k, l, a, temp;
	double t;

	for(i = 0; i < *n; i++) {
		for(l = 0; l < *L; l++) {
			for(a = 0; a < 2; a++) {
				
				temp = *(X + i*2*(*L) + a*(*L) + l);
				t = 0.0;	
				if(temp != 999) {
					for(k = 0; k < *K; k++) {
						*(store6 + k) = *(P + k*(*maxJ)*(*L) + (temp - 1)*(*L) + l) * *(Q + i*(*K) + k); 
						t += *(store6 + k);					     	      
					}
				}
				if(temp == 999) {
					for(k = 0; k < *K; k++) {
						*(store6 + k) = *(Q + i*(*K) + k); 
						t += *(store6 + k);					     	      
					}
				}
				for(k = 0; k < *K; k++) *(store6 + k) /= t;
				rmultinom_int_JM(store0, store6);
				*(Z + i*2*(*L) + a*(*L) + l) = *store0;
			}
		}	
	}
}

void logLik_calc(double *logLik_ans, int *X, double *P, double *Q, int *n, int *L, int *K, int *maxJ) {
	
	int i, k, l, a, temp;
	double t, tot1 = 1, tot2 = 0;
	
	for(i = 0; i < *n; i++) {
		for(a = 0; a < 2; a++) {
			for(l = 0; l < *L; l++) {
				temp = *(X + i*2*(*L) + a*(*L) + l);
				if(temp != 999) {
					t = 0.0;
					for(k = 0; k < *K; k++)
						t += (*(P + k*(*maxJ)*(*L) + (temp - 1)*(*L) + l)) * (*(Q + i*(*K) + k));
					tot1 *= t;
					
					if(tot1 < 1e-200) {
						tot2 += log(tot1);
						tot1 = 1;
					}
				}
			}
		}
	}
	tot2 += log(tot1);
	*logLik_ans = tot2;
}



void update_pi(double *P, double *pi, double *pi_temp, double *c, double *c_temp, double *c_sd, double *beta, int *K, int *L, int *J, int *maxJ,
	       double *store1, double *store2, double *store3) {
	
	int j, k, l;
	double t = 0.0, MH_p;
	
	for(l = 0; l < *L; l++) {
		
		t = 0.0;
		
		/* sample a new pi using draws from Dirichlet distributions */
		jump_pi(pi, pi_temp, &l, &t, L, J, maxJ, store1, beta); /* t = log( p(pi | pi_temp) / p(pi_temp | pi) )  */
		
		/* calculate log( p(P | pi_temp) / p(P | pi) ) */
		for(k = 0; k < *K; k++) {
			for(j =0; j < *(J + l); j++) {
				*(store1 + j) = *(pi + l * *maxJ + j) * ((1.0 / *(c + k)) - 1.0);
				*(store2 + j) = *(pi_temp + l * *maxJ + j) * ((1.0 / *(c + k)) - 1.0);
				*(store3 + j) = *(P + k*(*maxJ)*(*L) + j*(*L) + l); /* need to do this because P isn't in the right order */
			} 
			t -= ddirichlet_JM(store3, store1, (J + l));
			t += ddirichlet_JM(store3, store2, (J + l));
		}
		
		/* calculate Metropolis-Hastings ratio */
		MH_p = min(1.0, exp(t));
		
		/* decide whether to move or not */
		if(runif(0.0, 1.0) < MH_p) {
			for(j = 0; j < *maxJ; j++) {
				*(pi + l * *maxJ + j) = *(pi_temp + l * *maxJ + j);
			}
		}
	}
	
}

void jump_pi(double *pi, double *pi_temp, int *locus, double *log_ratio_dens, int *L, int *J, int *maxJ, double *store_temp, double *beta) {
	
	int j;
	
	for(j = 0; j < *(J + *locus); j++) *(store_temp + j) = *beta * *(pi + *locus * *maxJ + j);

	rdirichlet_JM((pi_temp + *locus * *maxJ), store_temp, (J + *locus));

	*log_ratio_dens -= ddirichlet_JM((pi_temp + *locus * *maxJ), store_temp, (J + *locus));

	for(j = 0; j < *(J + *locus); j++) *(store_temp + j) = *beta * *(pi_temp + *locus * *maxJ + j);

	*log_ratio_dens += ddirichlet_JM((pi + *locus * *maxJ), store_temp, (J + *locus));	
}

void update_cc(double *P, double *pi, double *pi_temp, double *c, double *c_temp, double *c_sd, int *K, int *L, int *J, int *maxJ,
	      double *store1, double *store2, double *store3, double *c_prior_mu, double *c_prior_sd) {
	
	int j, k, l;
	double t = 0.0, MH_p, alpha, beta;
	
	for(k = 0; k < *K; k++) {
		
		t = 0.0;
		
		/* sample a new c using draws from truncated Normal distributions */
		jump_c(c, c_temp, &k, c_sd, K, &t, 0.001);  /* t = log( p(c | c_temp) / p(c_temp | c) )  */
		
		/* calculate log( p(P | c_temp) / p(P | c) ) */
		
		for(l = 0; l < *L; l++) {
			for(j =0; j < *(J + l); j++) {
				*(store1 + j) = *(pi + l * *maxJ + j) * ((1.0 / *(c + k)) - 1.0);
				*(store2 + j) = *(pi + l * *maxJ + j) * ((1.0 / *(c_temp + k)) - 1.0);
				*(store3 + j) = *(P + k*(*maxJ)*(*L) + j*(*L) + l); /* need to do this because P isn't in the right order */
			}
			t -= ddirichlet_JM(store3, store1, (J + l));
			t += ddirichlet_JM(store3, store2, (J + l));
		}
		
		/* prior on c */
		
		alpha = (*c_prior_mu * *c_prior_mu) / (*c_prior_sd * *c_prior_sd);
		beta = (*c_prior_sd * *c_prior_sd) / *c_prior_mu;
		t += dgamma(c_temp[k], alpha, beta, 1) - dgamma(c[k], alpha, beta, 1);
		
		/* calculate Metropolis-Hastings ratio */
		MH_p = min(1.0, exp(t));
		
		/* decide whether to move or not */
		if(runif(0.0, 1.0) < MH_p) *(c + k) = *(c_temp + k);		
	}
}

void jump_c(double *c, double *c_temp, int *k, double *c_sd, int *K, double *log_ratio_dens, double low_lim) {
	
	double t = 1.0;
	
	do 
		*(c_temp + *k) = rnorm(*(c + *k), *c_sd);
	while ((*(c_temp + *k) > 1.0) || (*(c_temp + *k) < low_lim));
	
	t *= dnorm(*(c_temp + *k), *(c + *k), *c_sd, 0);
	t *= (pnorm(1, *(c + *k), *c_sd, 1, 0) - pnorm(0, *(c + *k), *c_sd, 1, 0));
	
	t /= dnorm(*(c + *k), *(c_temp + *k), *c_sd, 0);
	t /= (pnorm(1, *(c_temp + *k), *c_sd, 1, 0) - pnorm(0, *(c_temp + *k), *c_sd, 1, 0));
	
	*log_ratio_dens = log(t);
}

void rmultinom_int_JM(int *ans, double *probs){

	/* generate a single draw from a Multinomial(probs) distribution */
	int i;
	double unif;
	unif = runif(0.0, 1.0);
	i = 0;
	while(unif > 0) {
		*ans = (i + 1);
		unif -= *(probs + i);
		i += 1;
	}
}

void rdirichlet_JM(double *ans, double *probs, int *l){
	
	/* generate a single draw from a Dirichlet(probs) distribution */
	int i;
	double t = 0.0;
	for(i = 0; i < *l; i++) {
		*(ans + i) = rgamma(*(probs + i), 1.0);
		t += *(ans + i);
	}
	for(i = 0; i < *l; i++) *(ans + i) /= t;
}

double ddirichlet_JM(double *x, double *pars, int *l) {

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

