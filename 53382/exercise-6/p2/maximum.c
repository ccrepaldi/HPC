#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NMAX 20000000
#define M 20
#define K 2001

double  a[NMAX],x;

int main() {

  int i;
  clock_t t1,t2;

  for (i=0;i<NMAX;i++) 
    a[i]=exp((double)((i-NMAX)%M/(10.0*M)))+sin((double)(i%1000*M)/K);

  x=-1e20;
  t1=clock();
  for (i=0;i<NMAX;i++)
    if (a[i]>=x) x=a[i];
  t2=clock();

  printf("%20.10g %20.10g\n",((double)(t2-t1))/CLOCKS_PER_SEC,x);

  return 0;

}
