#include <math.h>
#include "taylor5a.h"

extern double jet_(double *x, double *yn, double *h);

double taylor5a(double x0, double y0, double xn, int N){
  double h = (xn - x0)/N;
  double yn = y0;
  int n;
  for (n =0; n<N; n++) {
    xn = x0 + n*h;
    yn = jet_(&xn, &yn, &h);

  }
  return yn;
}
