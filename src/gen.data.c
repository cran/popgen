#include <R.h>
#include <Rmath.h>

void rmultinom_int_JM1(int*, double*);

void rmultinom_int_JM1(int *ans, double *probs){

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


void gen_data(double *freq, int *n, int *L, int *new_data) {

  int i, j, tmp;
  
  /* set random number generator */
  GetRNGstate(); 
  
  for(i = 0; i < *n; i++) {
    for(j = 0; j < *L; j++) {
      rmultinom_int_JM1(&tmp, freq + j * 3);
      *(new_data + i * *L + j) = tmp - 1;
    }
  }

  /* close down random number generator */
  PutRNGstate();
}
  
