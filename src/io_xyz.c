#include "io_xyz.h"

FILE *open_output_xyzfile(int my_rank) {
  FILE *out_file = NULL;

  if (my_rank == 0) {
    out_file = fopen("gen.xyz", "w");
    if (!out_file) // check permissions
    {
      printf("***ERROR: cannot create geometry.out \n");
      exit(EXIT_FAILURE);
    }
  }

  return out_file;
}

void xyz_read(FILE* infile, Crystal *c){
  int i,coda;
  double *dptr;
  char buffer[LINE_SIZE+1],*ptr,*ptr2;
  int num_atoms;
  double latt[9];
  printf("read crystal'Z will be 1 \n");
  if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
  if (sscanf(buffer,"%d", &num_atoms)!=1)
    error_exit("Unable to read number of atoms");
  if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");

  /* We insist on a lattice being present, either in a Lattice string
   * on this line, or at the end of the file with "%PBC" here
   */
  coda=0;
  if (!strncmp(buffer,"%PBC",4))
    coda=1;
  else{
    ptr=strstr(buffer,"Lattice=");
    if (!ptr)
      error_exit("XYZ file does not include lattice specification");
    ptr+=strlen("Lattice=")+1;
    dptr=(double*)latt;
    if (sscanf(ptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",dptr,dptr+1,
           dptr+2,dptr+3,dptr+4,dptr+5,dptr+6,dptr+7,dptr+8)!=9)
      error_exit("Failed to parse lattice");
  }
  c->natoms = num_atoms;

  memcpy(&c->lattice_vectors[0][0], latt, 9*sizeof(double));
  c->coords=malloc(num_atoms*3*sizeof(double));
  c->atoms = malloc(num_atoms*sizeof(int));

  for(i=0;i<num_atoms;i++){
    if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
    ptr=buffer;
    while(isspace(*ptr)) ptr++;
    if (isdigit(*ptr)){
      if (sscanf(ptr,"%d %lf %lf %lf",c->atoms_num[i], c->coords[i*3], c->coords[i*3+1], c->coords[i*3+2])!=4)
          error_exit("Error parsing atom line");
      char *symbol = get_pte_label(c->atoms_num[i]);
      strcpy(&c->atoms[i*2], symbol);
    }
    else{
      ptr2=ptr;
      while((*ptr2)&&(!isspace(*ptr2))) ptr2++;
      *ptr2=0;
      c->atoms_num[i]=atsym2no(ptr);
      if (sscanf(ptr2+1,"%lf %lf %lf", c->coords[i*3],
         c->coords[i*3+1], c->coords[i*3+2])!=3)
    error_exit("Error parsing atom line");
    }
  }

  if (coda){
    i=0;
    while(1){
      if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
      ptr=buffer;
      while(isspace(*ptr)) ptr++;

      if (!(strncasecmp(ptr,"vec",3))){
        ptr2=ptr;
        while((*ptr2)&&(!isspace(*ptr2))) ptr2++;
        if (sscanf(ptr2+1,"%lf %lf %lf",latt[i*3],latt[i*3+1], latt[i*3+2])!=3)
          error_exit("Error parsing lattice vector");
        i++;
        if (i==3)
          {
            memcpy(&c->lattice_vectors[0][0], latt, 9*sizeof(double));
            break;
          }
      }
    }
  }
}


void xyz_write(FILE* outfile, Crystal *c){
  int i;
  char *fmt;
  int HIPREC=1;
  int num_atoms = c->natoms;


  if (HIPREC)
    fmt="%3s % 19.14f % 19.14f % 19.14f\n";
  else
    fmt="%3s % 14.8f % 14.8f % 14.8f\n";

  fprintf(outfile,"%d\n",num_atoms);

  fprintf(outfile,"Lattice=\"");
  for (i=0;i<3;i++)
    fprintf(outfile,"%.8f %.8f %.8f%s",c->lattice_vectors[i][0],c->lattice_vectors[i][1],
              c->lattice_vectors[i][2],(i==2)?"":" ");

  fprintf(outfile,"\" Properties=species:S:1:pos:R:3\n");

  for(i=0;i<num_atoms;i++)
  {

    fprintf(outfile,fmt, atno2sym(c->atoms_num[i]), c->coords[i*3], c->coords[i*3+1], c->coords[i*3+2]);
  }
//            atno2sym(m->atoms[i].atno),m->atoms[i].abs[0],
//                     m->atoms[i].abs[1],m->atoms[i].abs[2]);
  fflush(outfile);
}
