/*
  Eigenvalues of a square matrix using LAPACK DGEEV
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ii 127773
#define jj 2147483647
#define i1 16807
#define i2 2836
#define p 4.656612875e-10

/* A simple random number generator */

double ran(int seed) {
  static int ix;
  int  k1;
  if (seed!=-1 && seed!=0) ix=seed;
  k1 = ix/ii;
  ix = i1*(ix-k1*ii)-k1*i2;
  if (ix<0) ix = ix+jj;
  return ix*p;
}


void fill_matrix(double *A, int n, int seed) {
  int i;
  double x;
  x=ran(seed); /* This is only for setting the seed */
  for (i=0;i<n*n;i++) A[i]=ran(0);
  return;
}


int main(int argc, char *argv[])
{
  int n,i,j,*pivot,info,lwork;
  double *A,*wr,*wi,*vl,*vr,*work;
  char jobvl,jobvr;
  int seed;
  FILE *fp;
  
  /* Comamnd line arguments: matrix size, RNG seed */

  if (argc!=3) {
    fprintf(stderr,"usage: %s n seed\n",argv[0]);
    return -1;
  }
  n=atoi(argv[1]);
  seed=atoi(argv[2]);

  /* Allocate space. Note that matrix A is essentially a 1D array */

  A=(double *)malloc((size_t)n*n*sizeof(double));
  pivot=(int *)malloc((size_t)n*sizeof(int));
  wr=(double *)malloc((size_t)n*sizeof(double));
  wi=(double *)malloc((size_t)n*sizeof(double));
  vl=(double *)malloc((size_t)n*n*sizeof(double));
  vr=(double *)malloc((size_t)n*n*sizeof(double));
  lwork=10*n;
  work=(double *)malloc((size_t)lwork*sizeof(double));
  jobvl='N';
  jobvr='N';

  /* Fill the matrix with random numbers */

  fill_matrix(A,n,seed);
  
  /* ---- The eigenvalue calculation proper ---- */

  dgeev_(&jobvl,&jobvr, &n, A, &n,  wr, wi, vl, &n,   vr, &n,   work, &lwork, &info);

  /* Print eigenvalues to file "evalues.datc" */

  fp=fopen("evalues.datc","w");
  for (i=0;i<n;i++) fprintf(fp,"%12.8g %12.8g\n", wr[i],wi[i]);	
  fclose(fp);

  return 0;
}  

