#include "math.h"
#include "f.h"

inline double fexp(double x1, double x2){
  return exp(-1*pow(x1,2) - pow(pow(x1,2)-x2, 2));
}

double f(double x1, double x2) {
  return -1*fexp(x1,x2);
}

// Primeiras derivadas
double dfx1(double x1, double x2) {
  return
    2 *
    (2*(pow(x1,2) - x2)*x1 + x1) *
    fexp(x1,x2)
    ;
}

double dfx2(double x1, double x2){
  return -2*(pow(x1,2) - x2)*fexp(x1,x2);
}

// Segundas derivadas
double ddfx1x1(double x1, double x2){
  return
    fexp(x1,x2)*
    (
     (-4)*
     pow(2*(pow(x1,2) - x2)*x1 + x1,2)+
     2*
     (6*pow(x1,2) - 2*x2 + 1)
   );
}
double ddfx1x2(double x1, double x2){
  return
    fexp(x1,x2) *
    (
     4 * (2*(pow(x1,2) - x2)*x1 + x1)* (pow(x1,2) - x2)
     -
     4*x1
   );
}
double ddfx2x1(double x1, double x2){
  return ddfx1x2(x1,x2);
}
double ddfx2x2(double x1, double x2){
  return fexp(x1,x2)*(-4*pow(pow(x1,2) - x2,2) + 2);
}
