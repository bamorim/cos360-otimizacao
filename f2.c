#include "math.h"
#include "f.h"
#define DELTA 3e-15

// Primeiras derivadas
// Derivadas não são continuas, usar diferenças continuas
double _dfx1(double x1, double x2){
  return (f(x1 + DELTA, x2) - f(x1 - DELTA, x2))/(DELTA*2);
}

double _dfx2(double x1, double x2){
  return (f(x1, x2 + DELTA) - f(x1, x2 - DELTA))/(DELTA*2);
}

// Segundas derivadas
double _ddfx1x1(double x1, double x2){
  return (dfx1(x1 + DELTA, x2) - dfx1(x1 - DELTA, x2))/(DELTA*2);
}
double _ddfx1x2(double x1, double x2){
  return (dfx1(x1, x2 + DELTA) - dfx1(x1, x2 - DELTA))/(DELTA*2);
}
double _ddfx2x1(double x1, double x2){
  return (dfx2(x1 + DELTA, x2) - dfx2(x1 - DELTA, x2))/(DELTA*2);
}
double _ddfx2x2(double x1, double x2){
  return (dfx2(x1, x2 + DELTA) - dfx2(x1, x2 - DELTA))/(DELTA*2);
}

double f(double x1, double x2) {
  return sqrt(pow(pow(x1,2) - x2,2) + pow(x1,2));
}

// Primeiras derivadas
double dfx1(double x1, double x2) {
  double real = (2*(pow(x1,2) - x2)*x1 + x1)/f(x1,x2);
  if(real != real) {
    return _dfx1(x1,x2);
  } else {
    return real;
  }
}

double dfx2(double x1, double x2){
  double real = 0-(pow(x1,2) - x2)/f(x1,x2);
  if(real != real){
    return _dfx2(x1, x2);
  } else {
    return real;
  }
}

// Segundas derivadas
double ddfx1x1(double x1, double x2){
  double real = -pow(2*(pow(x1,2) - x2)*x1 + x1,2)/pow(pow(pow(x1,2) - x2,2) + pow(x1,2), 3/2) + (6*pow(x1,2) - 2*x2 + 1)/f(x1,x2);
  if(real != real){
    return _ddfx1x1(x1, x2);
  } else {
    return real;
  }
}
double ddfx1x2(double x1, double x2){
  double real = (2*(pow(x1,2) - x2)*x1 + x1)*(pow(x1,2) - x2)/pow( pow(pow(x1,2) - x2,2) + pow(x1,2),3/2) - 2*x1/f(x1,x2);

  if(real != real){
    return _ddfx1x2(x1, x2);
  } else {
    return real;
  }
}
double ddfx2x1(double x1, double x2){
  return ddfx1x2(x1,x2);
}
double ddfx2x2(double x1, double x2){
  double real = -pow(pow(x1,2) - x2,2)/pow(pow(pow(x1,2) - x2,2) + pow(x1,2),3/2) + 1/f(x1,x2);
  if(real != real){
    return _ddfx2x2(x1,x2);
  } else {
    return real;
  }
}
