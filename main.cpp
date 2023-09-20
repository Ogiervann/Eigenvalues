#include <iostream>
#include <cstdio>
#include <ctime>
#include "eigenvalues.h"


int main(int argc, char* argv[]){
  int n, m, k, res, its = 0;
  char* filename = 0;
  double eps, res1 = -1, res2 = -1, t1 = 0, t2 = 0;
  unsigned int start, end;

  double *a;
  double *x;
  double *y;
  double *z;


  if(!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n) == 1 && sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%lf", &eps) == 1 && sscanf(argv[4], "%d", &k) == 1) || (k == 0 && argc == 5)){
    printf("Usage %s n m eps k (file)\n", argv[0]);
    return 1;
  }

  if(k == 0){
    filename = argv[5];
  }


  a = new double[n*n];
  if(a == 0){
    printf("Not enough memory\n");
    return 1;
  }
  res = initMatrix(a, n, k, filename);

  if(res == 1){
    delete[] a;
    printf("File not found.\n");
    return 1;
  }
  if(res == 2){
    delete[] a;
    printf("Error in file.\n");
    return 2;
  }

  printMatrix(a, n, n, m);
  printf("\n");


  x = new double[n];
  y = new double[n];
  z = new double[n];

  if(x == 0 || y == 0 || z == 0){
    printf("Not enough memory\n");
    delete[] a;
    if(x != 0)delete[] x;
    if(y != 0)delete[] y;
    if(z != 0)delete[] z;
    return 1;
  }

  start = clock();
  res = TridiagonalAlgo(n, a, x, y, z);
  end = clock();
  t1 = (double)(end - start)/CLOCKS_PER_SEC;

  //res = ??; обработка ошибок


  start = clock();
  its = EigenvaluesAlgo(n, a, x, eps);
  end = clock();
  t2 = (double)(end - start)/CLOCKS_PER_SEC;
  //res = ??; обработка ошибок
  printMatrix(x, 1, n, m);

  initMatrix(a, n, k, filename);
  res1 = Residual1(n, a, x);
  res2 = Residual2(n, a, x);

  printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d \
  Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
  argv[0], res1, res2, its, its / n, t1, t2);
  delete[] a;
  delete[] x;
  delete[] y;
  delete[] z;
  return 0;
}
