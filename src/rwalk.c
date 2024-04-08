// gcc -o rwalk rwalk.c -lm
// ./rwalk [maxiter] [interval]
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#define L 512

int main(int argc, char **argv)
{
  int i, k, interval = 10000, maxiter = 100;
  double x, y, S=sqrt(1.0/L);
  FILE *GP;

  if (argc > 1)  maxiter = atoi(argv[1]);
  if (argc > 2)  interval = atoi(argv[2]);
  srand48(time(NULL));

  GP = popen("gnuplot ","w");
  fprintf(GP,
    "set term qt\n"
    "set xrange[-1:1]\n"
    "set yrange[-1:1]\n"
    "set nokey\n"
    "set title '2D Random Walk' font 'Helvetica,24'\n");
  fflush(GP);
  usleep(100000);

  for (k = 0; k < maxiter; k++){
    fprintf(GP, "plot '-' with lines lt 1\n");
    fflush(GP);
    usleep(interval);
    x=0; y=0;
    for (i = 0; i < L; i++){
      fprintf( GP, "%f %f\n",
              x += S*(drand48()-0.5),  y += S*(drand48()-0.5) );
    }
    fprintf(GP,"e\n");
    fflush(GP);
    usleep(interval);
  }
  fclose(GP);
  return 0;
}
