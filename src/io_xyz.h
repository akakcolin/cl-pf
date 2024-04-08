#ifndef IO_XYZ_H
#define IO_XYZ_H

#include "defines.h"

#define LINE_SIZE 300

void xyz_write(FILE* outfile, Crystal *c);
void xyz_read(FILE* infile, Crystal *c);

FILE *open_output_xyzfile(int my_rank);

#endif
