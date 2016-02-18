/*
   Solving the matrix equation A*x=b by calling LAPACK from C.
   Matrix A and the RHS vector b are read in from stdin.
   The 1st input line should give the problem size.
   For example a 3x3 problem:

3
0.60379247919382  0.27218792496996  0.19881426776106
0.74678567656443  0.44509643228795  0.93181457846166
0.41864946772751  0.84622141782432  0.52515249630517
0.42889236534100  0.30461736686939  0.18965374754717

   A. Kuronen, 2007
   Modified from: http://www.nacse.org/demos/lapack/codes/linear-c.html 
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
  int n,i,j,c1,c2,*pivot,ok;
  double *A,*b;
  
  scanf("%d",&n);
  printf("\nn %d\n",n);
  A=(double *)malloc((size_t)n*n*sizeof(double));
  b=(double *)malloc((size_t)n*sizeof(double));
  pivot=(int *)malloc((size_t)n*sizeof(int));

  for (i=0;i<n;i++) for (j=0;j<n;j++) scanf("%lg",&A[j*n+i]);
  for (i=0;i<n;i++) scanf("%lg",&b[+i]);

  printf("A\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) printf("%12.8g ",A[j*n+i]);
    printf("\n");
    }
  printf("\nb\n");
  for (i=0;i<n;i++) printf("%12.8g ",b[i]);
  printf("\n");

  c1=n;
  c2=1;
  
  dgesv_(&c1, &c2, A, &c1, pivot, b, &c1, &ok);      
  
  printf("x\n");
  for (i=0;i<n;i++) printf("%12.8g ", b[i]);	
  printf("\n");

  return 0;
}  

