/*
   Compute eigenvalues of a real matrix by calling LAPACK from C.
   Matrix is read in from stdin.
   The 1st input line should give the problem size.
   For example a 3x3 problem:

3
0.60379247919382  0.27218792496996  0.19881426776106
0.74678567656443  0.44509643228795  0.93181457846166
0.41864946772751  0.84622141782432  0.52515249630517

   A. Kuronen, 2007
   Modified from: http://www.nacse.org/demos/lapack/codes/eigen-c.html 
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main()
{
  int n,i,j,*pivot,info,lwork;
  double *A,*wr,*wi,*vl,*vr,*work;
  char jobvl,jobvr;
  
  scanf("%d",&n);
  printf("\nn %d\n",n);
  A=(double *)malloc((size_t)n*n*sizeof(double));
  pivot=(int *)malloc((size_t)n*sizeof(int));
  wr=(double *)malloc((size_t)n*sizeof(double));
  wi=(double *)malloc((size_t)n*sizeof(double));
  vl=(double *)malloc((size_t)n*n*sizeof(double));
  vr=(double *)malloc((size_t)n*n*sizeof(double));
  lwork=10*n;
  work=(double *)malloc((size_t)lwork*sizeof(double));
  
  for (i=0;i<n;i++) for (j=0;j<n;j++) scanf("%lg",&A[j*n+i]);

  printf("A\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) printf("%12.8g ",A[j*n+i]);
    printf("\n");
  }
  printf("\n");

  jobvl='N';
  jobvr='N';
  
  dgeev_(&jobvl,&jobvr, &n, A, &n,  wr, wi, vl, &n,   vr, &n,   work, &lwork, &info);
  /*     JOBVL, JOBVR,  N,  A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO ) */

  
  printf("INFO: %d\n",info);
  printf("Eigenvalues\n");
  for (i=0;i<n;i++) printf("%12.8g %12.8g\n", wr[i],wi[i]);	
  printf("\n");

  return 0;
}  

