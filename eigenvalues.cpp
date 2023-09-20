#include "eigenvalues.h"

int TridiagonalAlgo(int n, double* a, double* x, double*  y, double* z){
  int k;

  double norm = normMatr(a, n);

  for(k = 1; k < n-1; k++){
    tdStep(k, n, a, EPS * norm, x, y, z);
  }
  //printMatrix(a, n, n, n);
  //printf("\n");

  return 0;


}

int tdStep(int k, int n, double* a, double eps, double* x, double* y, double* z){
  double sk = 0;
  //double tr1 = 0, tr2 = 0;
  int i, j;
/*
  for(i = 0; i < n; i++){
    tr1 += a[i*n+i];
  }
*/
  for(i = k+1; i < n; i++){
    sk += a[i*n+k-1] * a[i*n+k-1];
  }
double norm_a1;
  if ((sk + a[k*n+k-1] * a[k*n+k-1])>0)
  {
    norm_a1 = sqrt(sk + a[k*n+k-1] * a[k*n+k-1]);
  }
  else
  {
    return 0;
  }

  if(norm_a1 < eps){
    return 0;
  }
  x[0] = a[(k-1)*n + k] - norm_a1;
  for(i = 1; i < n-k; i++){
    x[i] = a[(k-1)*n+k+i];
  }
  double norm_x;
  if ((sk+x[0]*x[0])>0)
  {
    norm_x = sqrt(sk + x[0] * x[0]);
  }
  else
  {
    return 0;
  }


  if(norm_x < eps){
    return 0;
  }

  for(i = 0; i < n-k; i++){
    x[i] /= norm_x;
  }



  for(i = 0; i < n-k; i++){
    y[i] = 0;
    for(j = 0; j < n-k; j++){
      y[i] += a[(i+k)*n + j+k] * x[j];
    }
  }

  double dot = dotProduct(x, y, n-k);

  for(i = 0; i < n-k; i++){
    z[i] = 2*(y[i] - dot * x[i]);
  }

  for(i = 0; i < n-k; i++){
    for(j = 0; j < n-k; j++){
      a[(i+k)*n+j+k] -= x[i]*z[j] + x[j]*z[i];
    }
  }


  a[k*n+(k-1)] = norm_a1;
  a[(k-1)*n+k] = norm_a1;

  for(i = k+1; i < n; i++){
    a[(k-1)*n+i] = 0;
    a[i*n+(k-1)] = 0;
  }
/*
  for(i = 0; i < n; i++){
    tr2 += a[i*n+i];
  }
*/
  //printf("%lf %lf\n", tr1, tr2);

  return 0;

}


int EigenvaluesAlgo(int n, double* A, double* x,double eps){
  int k=1, res;
  int its = 0;
  double norm = normMatr(A, n);
  while(k <= n){
    res=kEigenvalue(k, n, A, x, eps, norm, &its);
    if(res <= 0){
      k += n;
    }
    k+=res;
  }

  return its;
}


int kEigenvalue(int k, int n, double* A, double *x, double eps, double norm, int* its){
  double b = 2*normMatr(A, n), a = -b, c = (b+a)/2;
  int na = signChNum(a, A, n, norm*EPS), nb = signChNum(b, A, n, norm*EPS), nc = signChNum(c, A, n, norm*EPS);
  int i, it = 0;
  //printf("%lf %lf %lf\n", a, c, b);
  //printf("%d %d %d\n", na, nc, nb);

  while(b - a > eps*norm){
    it++;
    if(b-c < eps*norm || c-a < eps*norm){
    //  printf("Here\n");
      a = b;
      continue;
    }
    //printf("%10.3e %10.3e %10.3e\n", b-a, b-c, c-a);
    //printf("%lf %lf %lf\n", a, c, b);
    //printf("%d %d %d\n", na, nc, nb);

    if(nc < k){
      a = c;
      na = nc;
    }
    else{
      b = c;
      nb = nc;
    }
    c = a/2 + b/2;
    nc = signChNum(c, A, n, norm * EPS);
  }
  int res = nb-na;
  if(res <= 0){
    printf("Dosadno %lf %lf %d %d\n", b, a, nb, na);
    return res;
  }
  for(i = k-1; i < k+res-1; i++){
    x[i] = c;
  }
  *its += it;
  return res;
}


int signChNum(double m, double* a, int n, double eps){
  //p.23, p.96

  double l = a[0]-m;
  int res = (l < 0 ? 1 : 0), i;
  //printf("%lf:\n0 %lf\n",m, l);
  for(i = 1; i < n; i++){
    if(fabs(l) < eps) l = eps;
    l = a[i*n+i]-m -a[i*n+i-1]*a[(i-1)*n+i]/l;
    //printf("%d %lf\n", i, l);
    if(l < 0) res++;
  }

  return res;
}

double dotProduct(double* a, double* b, int n){
  double res = 0;
  for(int i = 0; i < n; i++){
    res += a[i] * b[i];
  }
  return res;
}

double normMatr(double* a, int n){
  double norm = 0, tmp = 0;
  int i, j;
  for(i = 0; i < n; i++){
    tmp = 0;
    for( j = 0; j < n; j++){
      tmp += fabs(a[i*n+j]);
    }
    norm = (tmp > norm? tmp: norm);
  }
  return norm;
}

double Residual1(int n, double* a, double* x){
  double sum = 0, norm = normMatr(a, n);
  int i;

  for(i = 0; i < n; i++){
    sum += a[i*n+i] - x[i];
  }

  if(fabs(norm) > 0) return fabs(sum)/fabs(norm);
  return -1;
}

double Residual2(int n, double* a, double* x){
  double norm = normMatr(a, n);

  double normA = 0, normx =0;
  int i, j;
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      normA += a[i*n+j]*a[j*n+i];
    }
    normx += x[i] * x[i];
  }
  if ((normA>0)&&(normx>0))
  {
    normA = sqrt(normA);
    normx = sqrt(normx);
  }


  if(fabs(norm) > 0)return fabs(normA-normx)/norm;
  return -1;
}
