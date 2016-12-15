#include "stdio.h"
#include "f.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#define EPSILON 1e-100
#define AUREA_EPSILON 1e-5
#define AUREA_RO 1e-2

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
void mul(double * m, double * v, double * o){
  // Only 2x2 * 2x1, don't bother implementing algorithm for this hahaha
  o[0] = m[0]*v[0] + m[1] * v[1];
  o[1] = m[2]*v[0] + m[3] * v[1];
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

inline double lanterna_verde(double * xv, double * d, double x) {
  return f(xv[0] + x*d[0], xv[1] + x*d[1]);
}

double aurea(double * x, double * d){
  double r = sqrt( pow(d[0],2) + pow(d[1],2) );
  double ud[2] = {
    d[0]/r,
    d[1]/r
  };
  double a = 0, s = AUREA_RO, b = 2*s;
  double t1 = (sqrt(5) - 1)/2;
  double t2 = 1 - t1;

  while ( lanterna_verde(x,d,a) < lanterna_verde(x,d,b) ){
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
void grad_method(double * x){
  // g = Gradiente
  double g[2] = {0,0};
  double t = 0;
  while(1){
    gf(x,g);

    if(fabs(g[0]) < EPSILON && fabs(g[1]) < EPSILON)
      return;

    // get t
    g[0] = - g[0];
    g[1] = - g[1];
    t = aurea(x, g);

    // x' = x + d*t = x - g*t
    x[0] = x[0] + g[0]*t;
    x[1] = x[0] + g[1]*t;
  }
}

void printm(double * mat, int m, int n) {
  for (int i = 0; i < m; i++){
      for (int j = 0; j < n; j++){
      printf("%g ", mat[i*m+j]);
    }
    printf("\n");
  }
  printf("\n\n");
}

void newt_method(double * x){
  double g[2] = {0,0};
  double gg[4] = { 0, 0,
                   0, 0 };

  double d[2] = {0,0};
  double t;
  while(1){
    gf(x,g);

    if(fabs(g[0]) < EPSILON && fabs(g[1]) < EPSILON)
      return;

    ggf(x,gg);
    invert(gg,gg);
    mul(gg,g,d);

    d[0] = -d[0];
    d[1] = -d[1];
    t = aurea(x, d);

    x[0] = x[0] + t*d[0];
    x[1] = x[1] + t*d[1];
  }
}

void run(void (*runner)(double * x)){
  double x[2] = {0.4,0.4};
  runner(x);
  printf("x* = (%.10e,%.10e)\n", x[0], x[1]);
  printf("f(x*) = (%.10e)\n", f(x[0], x[1]));
}

int main(){

  printf("Gradiente\n");
  run(grad_method);
  printf("Newton\n");
  run(newt_method);

  return 0;
}
