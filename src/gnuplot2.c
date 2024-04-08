#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define DIV 200
#define SEC 1000000
#define WAIT 10 * SEC

int main(int argc, char **argv) {
  int i;
  double x, y;
  FILE *GPLT;

  GPLT = popen("gnuplot", "w");
  fprintf(GPLT, "plot ’-’\n");
  for (i = 0; i < DIV; i++) {
    x = i / (double)DIV;
    y = cos(4 * M_PI * x);
    fprintf(GPLT, "%f %f\n", x, y);
  }
  fprintf(GPLT, "e\n");
  fflush(GPLT);
  usleep(WAIT);
  pclose(GPLT);

  return 0;
}
