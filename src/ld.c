#include <R.h>
#include <Rmath.h>
#include <math.h>

void ld_JM(int *n, int *L, int *flag, int *data, double *mat) {

  int i, j, l, l1, l2;
  double f00, f01, f10, f11, fA1, fA0, fB1, fB0, DAB, nn;
  double tab[3][3], p, q, r, s, V, xx, yy, ww, zz;
  double delta, x, y, z, w, dprime = 0.0, r2;

  if(*flag == 0) { /*  haplotype data */

    for(i = 0; i < (*L - 1); i++) {
      for(j = i + 1; j < *L; j++) {

	f00 = 0.0; 	f01 = 0.0; 	f10 = 0.0; 	f11 = 0.0;
	fA1 = 0.0;	fA0 = 0.0; 	fB1 = 0.0;	fB0 = 0.0; 
	DAB = 0.0; 	nn = 0.0;       dprime = 0.0;   r2 = 0.0;
	
	for(l = 0; l < *n; l++) {
	  
	  if(*(data + l * *L + i) == 0) {
	    
	    if(*(data + l * *L + j) == 0) {
	      fA0 += 1.0;
	      fB0 += 1.0;
	      f00 += 1.0;
	      nn += 1.0;
	    } else {
	      fA0 += 1.0;
	      fB1 += 1.0;
	      f01 += 1.0;
	      nn += 1.0;
	    }
	  } else {
	    if(*(data + l * *L + j) == 0) {
	      fA1 += 1.0;
	      fB0 += 1.0;
	      f10 += 1.0;
	      nn += 1.0;
	    } else {
	      fA1 += 1.0;
	      fB1 += 1.0;
	      f11 += 1.0;
	      nn += 1.0;
	    }
	  }
	}

	f00 /= nn; f01 /= nn; f10 /= nn; f11 /= nn;
	fA0 /= nn; fA1 /= nn; fB0 /= nn; fB1 /= nn;
	DAB = f11 - fA1 * fB1;
	
	if (DAB < 0) {
	  dprime = -1.0 * DAB / fmin2(fA1 * fB1, fA0 * fB0);
	}
	if (DAB > 0) {
	  dprime = DAB / fmin2(fA0 * fB1, fA1 * fB0);
	}
	if (DAB == 0) {
	  dprime = 0;
	  if (fmin2(fA1 * fB1, fA0 * fB0) == 0.0)
	    dprime = 1.0;
	} 
	
	r2 = (DAB * DAB) / (fA1 * fB1 * fA0 * fB0);

	*(mat + i * *L + j) = dprime;
	*(mat + j * *L + i) = r2;
      }
    }

  } else { /* genotype data */

    for(i = 0; i < (*L - 1); i++) {
      for(j = i + 1; j < *L; j++) {
	
	for(l1 = 0; l1 < 3; l1++) 
	  for(l2 = 0; l2 < 3; l2++) 
	    tab[l1][l2] = 0.0;
	
	for(l = 0; l < *n; l++) {

	  if(*(data + l * *L + i) == 0) {
	    if(*(data + l * *L + j) == 0) tab[0][0] += 1.0;
	    if(*(data + l * *L + j) == 1) tab[1][0] += 1.0;
	    if(*(data + l * *L + j) == 2) tab[2][0] += 1.0;
	  }
	  if(*(data + l * *L + i) == 1) {
	    if(*(data + l * *L + j) == 0) tab[0][1] += 1.0;
	    if(*(data + l * *L + j) == 1) tab[1][1] += 1.0;
	    if(*(data + l * *L + j) == 2) tab[2][1] += 1.0;
	  }
	  if(*(data + l * *L + i) == 2) {
	    if(*(data + l * *L + j) == 0) tab[0][2] += 1.0;
	    if(*(data + l * *L + j) == 1) tab[1][2] += 1.0;
	    if(*(data + l * *L + j) == 2) tab[2][2] += 1.0;
	  }
	}

	nn = tab[0][0] + tab[0][1] + tab[0][2] + tab[1][0] + tab[1][1] + tab[1][2] + tab[2][0] + tab[2][1] + tab[2][2];

	V = tab[1][1];
	
	xx = 2.0 * tab[0][0] + tab[0][1] + tab[1][0];
	yy = tab[0][1] + 2.0 * tab[0][2] + tab[1][2];
	ww = tab[1][0] + 2.0 * tab[2][0] + tab[2][1];
	zz = tab[1][2] + tab[2][1] + 2.0 * tab[2][2];
	
	p = (tab[1][0] + tab[1][1] + tab[1][2] + 2.0 * (tab[2][0] + tab[2][1] + tab[2][2])) / (2.0 * nn);
	q = (tab[0][1] + tab[1][1] + tab[2][1] + 2.0 * (tab[1][2] + tab[0][2]  + tab[2][2])) / (2.0 * nn);
	
	delta = 0.0;

	for(l1 = 0; l1 < 100; l1++) {
	  
	  r = ((1.0 - p) * (1.0 - q) + delta) * (p * q + delta);
	  r = r / ( ((1.0 - p) * (1.0 - q) + delta) * (p * q + delta) + ((1.0 - p) * q - delta) * (p * (1.0 - q) - delta) );
	  s = 1.0 - r;
	  
	  x = xx + r * V;
	  y = yy + s * V;
	  w = ww + s * V;
	  z = zz + r * V;
	  p = (y + w) / (2.0 * nn);
	  q = (z + w) / (2.0 * nn);
	  delta = w / (2.0 * nn) - p * q;
	}
	

	r2 = (delta * delta) / (q * p * (1.0 - q) * (1.0 - p));
	
			
	if(delta > 0) dprime = fabs(delta) / fmin2(q * (1.0 - p), (1.0 - q) * p);
	if(delta < 0) dprime = fabs(delta) / fmin2((1.0 - q) * (1.0 - p), q * p);
	if(delta == 0.0) dprime = 0.0;

	*(mat + i * *L + j) = dprime;
	*(mat + j * *L + i) = r2;
      }
    }
  }
      
}













