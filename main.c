#include "stdio.h"
#include "f.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#define EPSILON 1e-10
#define AUREA_EPSILON 1e-3
#define AUREA_RO 1e-2
#define ITER_MAX 10000


#ifdef DEBUG
#define debug(a) printf a
#else
#define debug(a) (0)
#endif

inline void printm(double * mat, int n, int m) {
  #ifdef DEBUG
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      printf("%g ", mat[i*m+j]);
    }
    printf("\n");
  }
  printf("\n\n");
  #endif
}

// Gradientes
void gf(double * x, double * g){
  g[0] = dfx1(x[0], x[1]);
  g[1] = dfx2(x[0], x[1]);
}

void ggf(double * x, double * g){
  g[0] = ddfx1x1(x[0], x[1]);
  g[1] = ddfx1x2(x[0], x[1]);
  g[2] = ddfx2x1(x[0], x[1]);
  g[3] = ddfx2x2(x[0], x[1]);
}

// Helpers
void scalar_mul(double * a, double c, int n, double * o){
  for(int i = 0; i < n; i++){
    o[i] *= c;
  }
}
void add(double * a, double * b, int n, double * o){
  for(int i = 0; i < n; i++){
    o[i] = a[i] + b[i];
  }
}
void id(double * o, int n, int m){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      if(i==j){
        o[i*m + j] = 1;
      } else {
        o[i*m + j] = 0;
      }
    }
  }
}
void clone(double * a, int n, double * o){
  for(int i = 0; i < n; i++){
    o[i] = a[i];
  }
}
void mul(double * a, double * b, int n, int m, int p, double * o){
  // a: n x m
  // aik = a[i*m + k]
  // b: m x p
  // bkj = a[k*p + j]
  // o: n x p
  // oij = [i*p + j]
  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      o[i*p + j] = 0;
      for(int k = 0; k < m; k++){
        o[i*p + j] += a[m*i + k] * b[p*k + j];
      }
    }
  }
}

void invert(double * im, double * om){
	// Define the dimension n of the matrix
	// and the signum s (for LU decomposition)
	int n = 2;
	int i,j,s;

	// Define all the used matrices
	gsl_matrix * m = gsl_matrix_alloc (n, n);
	gsl_matrix * inverse = gsl_matrix_alloc (n, n);
	gsl_permutation * perm = gsl_permutation_alloc (n);

  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      gsl_matrix_set(m, i,j, im[i*2+j]);
    }
  }

	// Make LU decomposition of matrix m
	gsl_linalg_LU_decomp (m, perm, &s);

	// Invert the matrix m
	gsl_linalg_LU_invert (m, perm, inverse);

  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      om[i*2+j] = gsl_matrix_get(inverse, i,j);
    }
  }
}

inline double lanterna_verde(double * xv, double * d, double t) {
  return f(xv[0] + t*d[0], xv[1] + t*d[1]);
}

double aurea(double * x, double * d){
  double r = sqrt( pow(d[0],2) + pow(d[1],2) );
  double a = 0, s = AUREA_RO, b = 2*s;
  double t1 = (sqrt(5) - 1)/2;
  double t2 = 1 - t1;

  while ( lanterna_verde(x,d,a) < lanterna_verde(x,d,b) ){
    // We are on the crescent side already
    if(b > 100)
      return 0;
    a = s;
    s = b;
    b = 2*b;
  }

  double u = a + t1*(b-a);
  double v = a + t2*(b-a);

  while (b-a > AUREA_EPSILON) {
    if(lanterna_verde(x,d,u) < lanterna_verde(x,d,v)) {
      b = v;
      v = u;
      u = a + t1*(b-a);
    } else {
      a = u;
      u = v;
      v = a + t2*(b-a);
    }
  }
  // Implementar busca da seção áurea
  return (u+v)/2;
}

// Métodos
long grad_method(double * x){

  long iter = 0;
  // g = Gradiente
  double g[2] = {0,0};
  double t = 0;
  while(1){
    iter++;
    gf(x,g);
    debug(("g:\n"));
    printm(g,1,2);
    debug(("X:\n"));
    printm(x,1,2);

    if(fabs(g[0]) < EPSILON && fabs(g[1]) < EPSILON || iter > ITER_MAX)
      return iter;

    // get t
    g[0] = - g[0];
    g[1] = - g[1];
    t = aurea(x, g);

    // x' = x + d*t = x - g*t
    x[0] = x[0] + g[0]*t;
    x[1] = x[0] + g[1]*t;
  }
}

long newt_method(double * x){
  long iter = 0;
  double g[2] = {0,0};
  double gg[4] = { 0, 0,
                   0, 0 };

  double d[2] = {0,0};
  double t;
  while(1){
    iter++;
    gf(x,g);

    if(fabs(g[0]) < EPSILON && fabs(g[1]) < EPSILON || iter > ITER_MAX)
      return iter;

    ggf(x,gg);
    invert(gg,gg);
    mul(gg,g,2,2,1,d);

    d[0] = -d[0];
    d[1] = -d[1];
    t = aurea(x, d);

    x[0] = x[0] + t*d[0];
    x[1] = x[1] + t*d[1];
  }
}

long bfgs(double * x){
  long iter = 0;
  double g[2];
  // H0 = I
  double h[4];

  double q[2];
  double p[2];
  double w1[4];
  double w2[4];
  double r;

  double d[2] = {0,0};
  gf(x,g);
  ggf(x,h);
  invert(h,h);
  double t;
  while(1){
    iter++;
    if(fabs(g[0]) < EPSILON && fabs(g[1]) < EPSILON || iter > ITER_MAX)
      return iter;

    mul(h, g, 2, 2, 1, d);
    scalar_mul(d, -1, 2, d);
    t = aurea(x, d);

    // d = t*d
    scalar_mul(d, t, 2, d);

    // p = -xk
    scalar_mul(x,-1, 2, p);
    // q = -gk
    scalar_mul(g,-1, 2, q);

    // xk = xk + tk*dk
    add(x, d, 2, x);
    // gk = \/f
    gf(x,g);

    // pk+1 = pk+1 - pk
    add(x,p,2,p);
    // qk+1 = gk+1 - gk
    add(g,q,2,q);
    // r = 1/qTp
    mul(q,p,1,2,1,w1);
    r = w1[0];
    // w1 = I - r*p*qT
    mul(p,q,2,1,2,w1);
    scalar_mul(w1,-1 * r,4,w1);
    id(w2,2,2);
    add(w2,w1,4,w1);

    // h = (I - r*p*qT)*H
    mul(w1,h,2,2,2,w2);
    clone(w2,4,h);

    // w1 = I - r*q*pT
    mul(q,p,2,1,2,w1);
    scalar_mul(w1,-1 * r,4,w1);
    id(w2,2,2);
    add(w2,w1,4,w1);

    // h = (I - r*p*qT)*H*(I - r*q*pT)
    mul(h,w1,2,2,2,w2);
    clone(w2,4,h);

    // w1 = r*p*pT
    mul(p,p,2,1,2,w1);
    scalar_mul(w1,r,4,w1);

    // h = (I - r*p*qT)*H*(I - r*q*pT) + r*p*pT
    add(h,w1,4,h);
  }
}

void run(long (*runner)(double * x)){
  double xs[10][2] = {
    {-0.6,0.7},
    {1,0},
    {0.001,0.001},
    {0.5, 0.5},
    {0.4, 0.6},
    {0.6, 0.4},
    {0.2, -0.2},
    {0.1, 0.1},
    {-0.1, -0.5},
    {0.3, 0.01}
  };

  for(int i = 0; i < 10; i ++){
    printf("x0 = (%.10e,%.10e)\n", xs[i][0], xs[i][1]);
    printf("f(x0) = (%.10e)\n", f(xs[i][0], xs[i][1]));
    long iter = runner(xs[i]);
    printf("iter = %d\n", iter);
    printf("x* = (%.10e,%.10e)\n", xs[i][0], xs[i][1]);
    printf("f(x*) = (%.10e)\n", f(xs[i][0], xs[i][1]));
    printf("---------\n");
  }
}

int main(){
  printf("Gradiente\n");
  run(grad_method);
  printf("\n\n\nNewton\n");
  run(newt_method);
  printf("\n\n\nBFGS\n");
  run(bfgs);

  printf("\n=============END=============\n\n");

  return 0;
}
