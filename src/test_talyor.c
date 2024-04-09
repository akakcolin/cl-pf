#include <stdio.h>
#include <stdlib.h>
#include "taylor5a.h"

int main(){
  int N;
  printf("%5s %16s\n", "N", "y");
  for (N=8; N<256; N*=2){
    double y = taylor5a(0,0,10, N);
    printf("%5d %16.11f\n", N, y);
  }
  exit(0);
}
