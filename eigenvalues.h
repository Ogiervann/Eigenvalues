#include <cmath>
#include "matrinit.h"

#define EPS 1e-16

int TridiagonalAlgo(int n, double* a, double* x, double*  y, double* z);
int EigenvaluesAlgo(int n, double* a, double* x,double eps);

int signChNum(double m, double* a, int n, double eps);

int kEigenvalue(int k, int n, double* A, double *x, double eps, double norm, int* its);


int tdStep(int k, int n, double* a, double eps, double* x, double* y, double* z);


double dotProduct(double* a, double* b, int n);
double normMatr(double* a, int n);

double Residual1(int n, double* a, double* x);
double Residual2(int n, double* a, double* x);
