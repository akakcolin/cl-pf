// https://github.com/StructuralBioinformaticsLab/libmol/blob/6d1c3f158a7dccadf6ee60098f1f1dbaf423be22/mol.0.0.6/sdf.c
//


#include "define.h"
#include <unistd.h>

void generate_random_rotation_matrix(double rot[3][3], uint64_t *rng) {
  double phi = 2 * PI * rng_uniform(rng);
  // theta, angle with x-axis. this should acos(u), -1<u<1
  double u = 2 * rng_uniform(rng) - 1;
  double theta = acos(u);
  // compute random axis with theta and phi
  double axis[] = {cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
  double psi = 2 * PI * rng_uniform(rng);

  double cosp = cos(psi);
  double vercosp = 1 - cosp; // Thanks S.L. Loney
  double sinp = sin(psi);

  // using Rodrigues rotation formula
  rot[0][0] = cosp + axis[0] * axis[0] * vercosp;
  rot[0][1] = axis[0] * axis[1] * vercosp - axis[2] * sinp;
  rot[0][2] = axis[0] * axis[2] * vercosp + axis[1] * sinp;

  rot[1][0] = axis[0] * axis[1] * vercosp + axis[2] * sinp;
  rot[1][1] = cosp + axis[1] * axis[1] * vercosp;
  rot[1][2] = axis[1] * axis[2] * vercosp - axis[0] * sinp;

  rot[2][0] = axis[0] * axis[2] * vercosp - axis[1] * sinp;
  rot[2][1] = axis[1] * axis[2] * vercosp + axis[0] * sinp;
  rot[2][2] = cosp + axis[2] * axis[2] * vercosp;
}



int **malloc_Mat2d_i(int n1, int n2){
    int i;
    int **array =(int **)malloc(n1*sizeof(int*));
    for(i=0;i<n1;i++){
        array[i]=(int *)malloc(n2*sizeof(int));
    }
    return array;
}

void free_Mat2d_i(int **array, int n1){
  int i;
  for( i=0; i< n1; i++){
    free(array[i]);
  }
  free(array);
}


void free_rgroup(RgroupRecord *rg){
  //  free(rg->name);
  free(rg->num_rgp);
  free(rg->rgp);
  free_Mat2d_i(rg->atomid, rg->num_mol);
  free_Mat2d_i(rg->groupid, rg->num_mol);
  free(rg);
}

void free_meta(Sdf_MetaData *meta){
  for(int i=0;i<meta->num_data; i++){
    Sdf_Ctab *ctab= meta->ctabs + i;
    free(ctab->atoms);
    free(ctab->bonds);
  }
  free(meta->ctabs);
  free(meta->counts);
  free(meta);
}

void copy_ctab_atom(Sdf_Ctab_Atom *a, Sdf_Ctab_Atom *b){
  a->pos[0] = b->pos[0];
  a->pos[1] = b->pos[1];
  a->pos[2] = b->pos[2];
  a->mass_diff = b->mass_diff;
  strcpy(a->symbol, b->symbol);
  a->valence = b->valence;
  a->charge = b->charge;
  a->global_id = b->global_id;
  a->atom_num = b->atom_num;
  a->mol_id = b->mol_id;
}

void copy_ctab_atoms(Sdf_Ctab_Atom *a, Sdf_Ctab_Atom *b, int num_atoms){
  for(int i=0; i<num_atoms; i++){
    if(b[i].global_id > -1){
    copy_ctab_atom(&a[i], &b[i]);
   }
  }
}

int check_cycle_pair(RgroupRecord *rg, int *cycle, int num_cycle){
  int i, j, k;
  int fatom, satom;
  int fpair, spair;
  int fi, si;
  int num_check=0;
  for(i=0; i<num_cycle; i++){
    fi = cycle[2*i]-1;
    si = cycle[2*i+1]-1;
    for(j=0; j<rg->num_rgp[fi]; j++){
      fpair = rg->groupid[fi][j];
      fatom = rg->atomid[fi][j];
      for(k=0; k<rg->num_rgp[si]; k++){
        spair = rg->groupid[si][k];
        satom = rg->atomid[si][k];
        if(fpair == spair){
          num_check++;
        }
      }
    }
  }
  if(num_check == num_cycle) return 1;
  else
    return 0;
};

int gen_reorder_pair(int **pair_record,
                 int num_remove,
                 int *cycle,
                 int num_cycle,
                 RgroupRecord *rg,
                 int **reorder_pair){
  int i, j,k;
  int fpair, spair;
  int fi, si;

  int num_check=0;
  int find_same=0;
  for(j=0; j< num_remove; j++){

    for(i=0; i< num_cycle; i++){
      fi = cycle[2*i]-1;
      si = cycle[2*i+1]-1;
      if(pair_record[j][0] == fi && pair_record[j][2] == si){
        reorder_pair[num_check][0] = pair_record[j][0];
        reorder_pair[num_check][1] = pair_record[j][1];
        reorder_pair[num_check][2] = pair_record[j][2];
        reorder_pair[num_check][3] = pair_record[j][3];
        num_check++;
      }
    }
  }
  for(j=0; j< num_remove; j++){
    find_same=0;
    for(i=0; i< num_check; i++){ // carefule num_check
      fi = reorder_pair[i][0];
      si = reorder_pair[i][2];
      if(pair_record[j][0] == fi && pair_record[j][2] == si){
        find_same=1;
        break;
      }
    }
    if(find_same==1) continue;
    else{ // find new pair
      reorder_pair[num_check][0] = pair_record[j][0];
      reorder_pair[num_check][1] = pair_record[j][1];
      reorder_pair[num_check][2] = pair_record[j][2];
      reorder_pair[num_check][3] = pair_record[j][3];
      num_check++;
    }
  }
  if(num_check == num_remove) return 1;
  else return 0;
}

int get_pair( int **pair_record, int *max_total_atoms,
              int *max_total_bonds,
              RgroupRecord *rg,
              Sdf_Ctab_Counter *counters){
  int i,j,k,l;
  int fatom, satom;

  for(i=0; i<rg->num_mol; i++){
    max_total_atoms[0] += counters[i].num_atoms;
    max_total_bonds[0] += counters[i].num_bonds;
  }

  int num_remove=0;

  int fpair, spair;
  for(i=0; i<rg->num_mol; i++){
    for(j=0; j<rg->num_rgp[i]; j++){
      fpair = rg->groupid[i][j];
      fatom = rg->atomid[i][j];
      for(k=0; k<rg->num_mol, k!=i; k++){
        //if(k==i) continue;
        for(l=0; l<rg->num_rgp[k]; l++){
          spair = rg->groupid[k][l];
          satom = rg->atomid[k][l];
          if(fpair == spair){
            if (i > k) {
              pair_record[num_remove][0] = k;
              pair_record[num_remove][1] = satom;
              pair_record[num_remove][2] = i;
              pair_record[num_remove][3] = fatom;
              printf("find connect pair index %d:%d-%d:%d %d\n", k, l, i, j,
                     spair);
            } else {
              pair_record[num_remove][0] = i;
              pair_record[num_remove][1] = fatom;
              pair_record[num_remove][2] = k;
              pair_record[num_remove][3] = satom;
              printf("find connect pair index %d:%d-%d:%d %d\n", i, j, k, l,
                     spair);
            }
            num_remove++;
          }
        }
      }
    }
  }
  printf("find pair is %d\n", num_remove);
  return num_remove;
}


// normalize any number of degree to the range (-180, 180)
double normalize_angle(double theta){
  int f = (180-theta)/360;
  if(f<0){
    f = -1*f;
  }
  return (theta + f*360.0);
}

void vec4f_x_mat4(const double *a, const double *b, const int n, double *c)
{
  int t = 0;
  for (int i = 0; i < n; ++i)
  {
    t = i * 4;
    c[t] = a[t] * b[0] + a[t + 1] * b[4] + a[t + 2] * b[8] + a[t + 3] * b[12];
    c[t + 1] =
        a[t] * b[1] + a[t + 1] * b[5] + a[t + 2] * b[9] + a[t + 3] * b[13];
    c[t + 2] =
        a[t] * b[2] + a[t + 1] * b[6] + a[t + 2] * b[10] + a[t + 3] * b[14];
    c[t + 3] =
        a[t] * b[3] + a[t + 1] * b[7] + a[t + 2] * b[11] + a[t + 3] * b[15];
  }
}

void C1(const double *a, const int n, double *c)
{
  int t = 0;
  int m = 0;
  if (n == 3)
  {
    c[0] = a[0];
    c[1] = a[1];
    c[2] = a[2];
    c[3] = 1.0f;
    c[4] = a[3];
    c[5] = a[4];
    c[6] = a[5];
    c[7] = 1.0f;
    c[8] = a[6];
    c[9] = a[7];
    c[10] = a[8];
    c[11] = 1.0f;
    c[12] = 1.0f;
    c[13] = 1.0f;
    c[14] = 1.0f;
    c[15] = 1.0f;
  }
  else
  {
    for (int i = 0; i < n; ++i)
    {
      t = i * 3;
      c[t + m] = a[t];
      c[t + m + 1] = a[t + 1];
      c[t + m + 2] = a[t + 2];
      c[t + m + 3] = 1.0f;
      m += 1;
    }
  }
}

//TrimAcos: trim the acos argument to be within (-1, 1).
double TrimAcos(double acosarg)
{
  if (acosarg <= -1.0 + 1.0e-12 || acosarg >= 1.0 - 1.0e-12)
  {
    if (acosarg < -1.00001 || acosarg > 1.00001)
    {
      printf("TRIM_ACOS >> Error: ACOS argument is %12.4lf\nTRIM_ACOS >> "
             "This is erroneously large!\n",
             acosarg);
      exit(1);
    }
    if (acosarg >= 1.0 - 1.0e-12)
    {
      return 1.0 - 1.0e-12;
    }
    else
    {
      return -1.0 + 1.0e-12;
    }
  }
  return acosarg;
}


double Dihedral(const double *alocs, const double *blocs, const double *clocs, const double *dlocs)
{
  int i;
  double theta, acosarg;
  double ab[3], bc[3], cd[3], scr[3], crabbc[3], crbccd[3];

  for (i = 0; i < 3; ++i)
  {
    ab[i] = blocs[i] - alocs[i];
    bc[i] = clocs[i] - blocs[i];
    cd[i] = dlocs[i] - clocs[i];
  }
  crabbc[0] = ab[1] * bc[2] - ab[2] * bc[1];
  crabbc[1] = ab[2] * bc[0] - ab[0] * bc[2];
  crabbc[2] = ab[0] * bc[1] - ab[1] * bc[0];

  crbccd[0] = bc[1] * cd[2] - bc[2] * cd[1];
  crbccd[1] = bc[2] * cd[0] - bc[0] * cd[2];
  crbccd[2] = bc[0] * cd[1] - bc[1] * cd[0];

  double dmag_crabbc = sqrt(crabbc[0] * crabbc[0] + crabbc[1] * crabbc[1] + crabbc[2] * crabbc[2]);
  double dmag_crbccd = sqrt(crbccd[0] * crbccd[0] + crbccd[1] * crbccd[1] + crbccd[2] * crbccd[2]);
  double dotp = (crabbc[0] * crbccd[0] + crabbc[1] * crbccd[1] + crabbc[2] * crbccd[2]);

  acosarg = dotp / dmag_crabbc / dmag_crbccd;

  scr[0] = crabbc[1] * crbccd[2] - crabbc[2] * crbccd[1];
  scr[1] = crabbc[2] * crbccd[0] - crabbc[0] * crbccd[2];
  scr[2] = crabbc[0] * crbccd[1] - crabbc[1] * crbccd[0];

  double dotp_scr = (scr[0] * bc[0] + scr[1] * bc[1] + scr[2] * bc[2]);
  if (dotp_scr > 0.0)
  {
    theta = acos(TrimAcos(acosarg));
  }
  else
  {
    theta = -acos(TrimAcos(acosarg));
  }
  return theta;
}


double dot(double a[3], double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


double vector3_norm(double a[3]) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

void cross_vector3_vector3(double cross[3], double a[3], double b[3]) {
  cross[0] = a[1] * b[2] - a[2] * b[1];
  cross[1] = a[2] * b[0] - a[0] * b[2];
  cross[2] = a[0] * b[1] - a[1] * b[0];

  return;
}

void normalise_vector3(double a[3]) {
  double norm;
  norm = vector3_norm(a);
  a[0] /= norm;
  a[1] /= norm;
  a[2] /= norm;
}


void get_rotate_matrix_from_axis_angle(double rotmat[16], double axis[3], double angle){
  normalise_vector3(axis);
  double n1, n2, n3;
  double n12, n22, n32;
  double cosq, sinq;
  n1 = axis[0];
  n2 = axis[1];
  n3 = axis[2];
  n12 = n1 * n1;
  n22 = n2 * n2;
  n32 = n3 * n3;
  cosq = cos(angle);
  sinq = sin(angle);

  rotmat[0] = n12 + (1-n12)*cosq;
  rotmat[1] = n1 * n2 * (1-cosq) + n3 * sinq;
  rotmat[2] = n1 * n3 * (1-cosq) - n2 *sinq;
  rotmat[3] = 0.0;

  rotmat[4]= n1 * n2 * (1-cosq) - n3 * sinq;
  rotmat[5]= n22 + (1-n22)*cosq;
  rotmat[6]= n2 * n3 * (1-cosq) + n1 * sinq;
  rotmat[7]= 0.0;

  rotmat[8] = n1 * n3 * (1-cosq) + n2 * sinq;
  rotmat[9] = n2 * n3 * (1-cosq) - n1 * sinq;
  rotmat[10] = n32 + (1-n32)*cosq;
  rotmat[11] = 0.0;

  rotmat[12] = 0.0;
  rotmat[13] = 0.0;
  rotmat[14] = 0.0;
  rotmat[15] = 1.0;

}

void np_eye4(double *a)
{
  a[0] = 1.0;
  a[1] = 0.0;
  a[2] = 0.0;
  a[3] = 0.0;
  a[4] = 0.0;
  a[5] = 1.0;
  a[6] = 0.0;
  a[7] = 0.0;
  a[8] = 0.0;
  a[9] = 0.0;
  a[10] = 1.0;
  a[11] = 0.0;
  a[12] = 1.0;
  a[13] = 0.0;
  a[14] = 0.0;
  a[15] = 1.0;
}



double cart_dist(double p1[3], double p2[3]) {
  double c[3];
  vector3_subtract(p1, p2, c);
  return vector3_norm(c);
}


double getAngle(double a1[3], double a2[3], double a3[3]){
    double d1_2 = cart_dist(a1, a2); // distance atom1 to atom2
    double d2_3 = cart_dist(a2, a3); // distance atom2 to atom3
    double d3_1 = cart_dist(a3, a1); // distance atom3 to atom1

    double numerator = pow(d3_1, 2) - pow(d1_2, 2) - pow(d2_3, 2);
    double denominator = -2 * d1_2 * d2_3;

    double radians = acos(numerator / denominator);

    double degrees = radiansToDegrees(radians);

    if(degrees > 180){
        degrees -= 180;
    }

    return degrees;
}



// rotate a1 about the vector defined by a2 and a3
// axis_a  start point of axis
// axis_b  end point of axis
// theta  degree
// rotate points around axis
void rotate_along_axis(double *coords, int num_atoms,
                       double axis_a[3], double axis_b[3], double theta){
  int i;
  double axis[3];
  double t0[16];
  double rotmat[16];
  double *coord4s=NULL;
  double *t4=NULL;

  coord4s=(double *)malloc(sizeof(double)*4*num_atoms);
  t4=(double *)malloc(sizeof(double)*4*num_atoms);

  C1(coords, num_atoms, coord4s);

  // axis = axis_b - axis_a
  vector3_subtract(axis_b, axis_a, axis);
  np_eye4(t0);
  t0[12] = -axis_a[0];
  t0[13] = -axis_a[1];
  t0[14] = -axis_a[2];

  // 
  vec4f_x_mat4(coord4s, t0, num_atoms, t4);
  get_rotate_matrix_from_axis_angle(rotmat, axis, theta);
  vec4f_x_mat4(t4, rotmat, num_atoms, coord4s);
  t0[12] = axis_a[0];
  t0[13] = axis_a[1];
  t0[14] = axis_a[2];

  vec4f_x_mat4(coord4s, t0, num_atoms, t4);
  for(i=0;i<num_atoms;i++){
    coords[i*3] = t4[i*4];
    coords[i*3+1]= t4[i*4+1];
    coords[i*3+2] = t4[i*4+2];
  }

  free(coord4s);
  free(t4);
}


void sdf_gen_mol(Sdf_Ctab_Atom *atoms, int num_atoms, uint64_t *rng){
  int i;
  Sdf_Ctab_Atom *atom=NULL;
  double *coords=NULL;
  double rotmat[3][3];
  generate_random_rotation_matrix(rotmat, rng);

  coords=(double *)malloc(sizeof(double)*3*num_atoms);

  for(i=0;i<num_atoms; i++){
    atom = atoms+i;
    memcpy(&coords[i*3], atom->pos, 3*sizeof(double));
    vector3_mat3b3_multiply(rotmat, &coords[i*3], &coords[i*3]);
  }

  for(i=0;i<num_atoms;i++){
    atom = atoms+i;
    memcpy(atom->pos, &coords[i*3], 3*sizeof(double));
  }

  free(coords);

}
void sdf_rotate_along_axis(Sdf_Ctab_Atom *atoms, int num_atoms,
                       double axis_a[3], double axis_b[3], double theta){
  int i;
  double axis[3];
  double t0[16];
  double rotmat[16];
  double *coord4s=NULL;
  double *t4=NULL;
  double *coords=NULL;
  Sdf_Ctab_Atom *atom=NULL;

  coord4s=(double *)malloc(sizeof(double)*4*num_atoms);
  t4=(double *)malloc(sizeof(double)*4*num_atoms);
  coords=(double *)malloc(sizeof(double)*3*num_atoms);
  for(i=0;i<num_atoms; i++){
    atom = atoms+i;
    memcpy(&coords[i*3], atom->pos, 3*sizeof(double));
  }

  C1(coords, num_atoms, coord4s);

  // axis = axis_b - axis_a
  vector3_subtract(axis_b, axis_a, axis);
  np_eye4(t0);
  t0[12] = -axis_a[0];
  t0[13] = -axis_a[1];
  t0[14] = -axis_a[2];

  
  vec4f_x_mat4(coord4s, t0, num_atoms, t4);
  get_rotate_matrix_from_axis_angle(rotmat, axis, theta);
  vec4f_x_mat4(t4, rotmat, num_atoms, coord4s);
  t0[12] = axis_a[0];
  t0[13] = axis_a[1];
  t0[14] = axis_a[2];

  vec4f_x_mat4(coord4s, t0, num_atoms, t4);
  for(i=0;i<num_atoms;i++){
    atom = atoms+i;
    coords[i*3] = t4[i*4];
    coords[i*3+1]= t4[i*4+1];
    coords[i*3+2] = t4[i*4+2];
    memcpy(atom->pos, &coords[i*3], 3*sizeof(double));
  }
  free(coord4s);
  free(t4);
  free(coords);
}


void getNormalPlane(double normp[3], const Plane p){
      //normal vector = (point2 - point1) CROSS (point3 - point1)
    double ax = p.a2[0] - p.a1[0];
    double ay = p.a2[1] - p.a1[1];
    double az = p.a2[2] - p.a1[2];

    double bx = p.a3[0] - p.a1[0];
    double by = p.a3[1] - p.a1[1];
    double bz = p.a3[2] - p.a1[2];

    //calculate the cross product of atomA and atomB
    normp[0] = ay * bz - az * by;
    normp[1] = az * bx - ax * bz;
    normp[2] = ax * by - ay * bx;
}


double getDistance_bypoints(double a0[3], double a1[3], double a2[3], double a3[3]){
  double angle=getAngle(a1, a2, a3);
  double b1[3], b2[3], b0[3];
  double dotbb, denom;
  double t;
  double d2;
  if(fabs(angle)<0.0001){
    // determine if line pass through origin
    vector3_subtract(a1, a0, b1);
    vector3_subtract(a2, a1, b2);
    // dot product
    dotbb = dot(b1, b2);
    // find the square of the a1 and a2
    ARR3MUL(b2,b2);
    denom = b2[0] + b2[1] + b2[2];
    t = -1 * dotbb/denom;
    // find the distance^2 from the origin th the line generated by the atoms
    b0[0] = (a1[0] - a0[0]) + t*(a2[0] - a1[0]);
    b0[1] = (a1[1] - a0[1]) + t*(a2[1] - a1[1]);
    b0[2] = (a1[2] - a0[2]) + t*(a2[2] - a1[2]);
    d2 = dot(b0, b0);
    return sqrt(d2);
  }else
    {
      // a1 a2 a3 do not define a line
      return -1;
    }
}


Plane createPlane(double a1[3], double a2[3], double a3[3]){
    Plane p;
    ARR3SET(p.a1, a1);
    ARR3SET(p.a2, a2);
    ARR3SET(p.a3, a3);
    return p;
}


void getNormalPlane_bypoints(double a0[3], double a1[3], double a2[3], double a3[3]){
    //atoms are co-linear
    if(fabs(getAngle(a1, a2, a3)) < .00001){
      //will try to use 0,0,0 then 1,0,0 then 0,1,0 as normals
      //determine if line passes through origin
      a0[0]=0;
      a0[1]=0;
      a0[2]=0;
      double d2_origin = getDistance_bypoints(a0, a1, a2, a3);
      a0[0]=1;
      a0[1]=0;
      a0[2]=0;
      double d2_x1 = getDistance_bypoints(a0, a1, a2, a3);
      a0[0]=0;
      a0[1]=1;
      a0[2]=0;

      if (d2_origin > .00001) {  // use origin
        a0[0] = 0;
        a0[1] = 0;
        a0[2] = 0;
      } else if (d2_x1 > .00001) {  // use 1,0,0
        a0[0] = 1;
        a0[1] = 0;
        a0[2] = 0;
      } else {  // use 0,1,0
        a0[0] = 0;
        a0[1] = 1;
        a0[2] = 0;
      }
      double a[3];
      Plane p = createPlane(a0, a2, a1);
      getNormalPlane(a0, p);
    }
}


// rotate a1 about a2 in the plane defined by a1 a2 a3
void rotate_inplane(double *coords, int num_atoms,
                    double a2[3], double a3[3],
                       double theta){
  double a1[3];
  ARR3SET(a1, coords); // first atom in coords always is a1 !!
  Plane p1=createPlane(a1, a2, a3);
  double n[3];
  double angle=getAngle(a1, a2, a3);
  double vend[3];
  if(fabs(angle) < 0.0001){
    // find a line that is perpendicular to the line of atoms as the normal
    getNormalPlane_bypoints(n, a1, a2, a3);
  }else{
    getNormalPlane(n, p1);
  }
  vector3_add(a2, n, vend);
  // rotate about the normal vector
  rotate_along_axis(coords, num_atoms,
                    a2, vend, theta);
}


// rotate a1 about a2 in the plane defined by a1 a2 a3
void sdf_rotate_inplane(Sdf_Ctab_Atom *atoms,
                        int num_atoms,
                        double a1[3],
                        double a2[3],
                        double a3[3],
                        double theta){
  double *coords=NULL;
  coords=malloc(sizeof(double)*num_atoms*3);
  Sdf_Ctab_Atom *atom=NULL;
  for(int i=0; i< num_atoms; i++){
    atom = atoms+i;
    memcpy(&coords[i*3], atom->pos, 3*sizeof(double));
  }
  // first atom in coords always is a1 !!
  Plane p1=createPlane(a1, a2, a3);
  double n[3];
  double angle=getAngle(a1, a2, a3);
  double vend[3];
  if(fabs(angle) < 0.0001){
    // find a line that is perpendicular to the line of atoms as the normal
    getNormalPlane_bypoints(n, a1, a2, a3);
  }else{
    getNormalPlane(n, p1);
  }
  vector3_add(a2, n, vend);
  // rotate about the normal vector
  rotate_along_axis(coords, num_atoms,
                    a2, vend, theta);

  for(int i=0; i< num_atoms; i++){
    atom = atoms+i;
    memcpy(atom->pos, &coords[i*3], 3*sizeof(double));
  }
  free(coords);
}

double atom_dist2(Sdf_Ctab_Atom a, Sdf_Ctab_Atom b){
  double v[3];
  v[0] = a.pos[0] - b.pos[0];
  v[1] = a.pos[1] - b.pos[1];
  v[2] = a.pos[2] - b.pos[2];
  return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

double get_dist2(double a[3], double b[3]){
  double v[3];
  v[0] = a[0] - b[0];
  v[1] = a[1] - b[1];
  v[2] = a[2] - b[2];
  return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

int check_vec3_isNull(double a[3], double tol) {
  if (fabs(a[0]) < tol && fabs(a[1]) < tol && fabs(a[2]) < tol)
    return 1;
  else
    return 0;
}

void rotation_mat_around_axis(double rot[3][3], double axis[3], double psi) {
  double cosp = cos(psi);
  double vercosp = 1 - cosp; // Thanks S.L. Loney
  double sinp = sin(psi);

  // using Rodrigues rotation formula
  rot[0][0] = cosp + axis[0] * axis[0] * vercosp;
  rot[0][1] = axis[0] * axis[1] * vercosp - axis[2] * sinp;
  rot[0][2] = axis[0] * axis[2] * vercosp + axis[1] * sinp;

  rot[1][0] = axis[0] * axis[1] * vercosp + axis[2] * sinp;
  rot[1][1] = cosp + axis[1] * axis[1] * vercosp;
  rot[1][2] = axis[1] * axis[2] * vercosp - axis[0] * sinp;

  rot[2][0] = axis[0] * axis[2] * vercosp - axis[1] * sinp;
  rot[2][1] = axis[1] * axis[2] * vercosp + axis[0] * sinp;
  rot[2][2] = cosp + axis[2] * axis[2] * vercosp;
}

void rotation_matrix_from_vectors(double rotmat[3][3], double a[3],
                                  double b[3]) {
  // find axis and angle
  normalise_vector3(a);
  normalise_vector3(b);
  double axis[3];
  cross_vector3_vector3(axis, a, b);
  double angle;
  angle = acos(dot(a, b));

  // if and b are the same, return identity
  if (check_vec3_isNull(axis, TOL)) {
    rotmat[0][0] = 1;
    rotmat[0][1] = 0;
    rotmat[0][2] = 0;

    rotmat[1][0] = 0;
    rotmat[1][1] = 1;
    rotmat[1][2] = 0;

    rotmat[2][0] = 0;
    rotmat[2][1] = 0;
    rotmat[2][2] = 1;

    return;
  }
  // find rotation matrix from axis and angle
  normalise_vector3(axis);
  rotation_mat_around_axis(rotmat, axis, angle);
}

/* rodrigues_rotation()
 *  Applies rotation of theta about vector *v to vector *a (in situ).
 *  Vector *v should be from offset (eg if defined as bond, you need to
 *  subtract that offset from vector *a prior to applying rotation).
 * arg *a: Vector to be rotated.
 * arg *v: Vector to rotate about.
 * arg theta: Angle (in radians)
 * returns: 1 in all cases
 */
int rodrigues_rotation(double a[3], double v[3], const double theta) {
  double offset[3];
  ARR3SET(offset, a);

  double R[3][3];
  double ct = cos(theta);
  double st = sin(theta);
  R[0][0] = cos(theta) + powf(v[0], 2.) * (1 - ct);
  R[0][1] = v[0] * v[1] * (1 - ct) - v[2] * st;
  R[0][2] = v[1] * st + v[0] * v[2] * (1 - ct);

  R[1][0] = v[2] * st + v[0] * v[1] * (1 - ct);
  R[1][1] = ct + powf(v[1], 2.) * (1 - ct);
  R[1][2] = -v[0] * st + v[1] * v[2] * (1 - ct);

  R[2][0] = -v[1] * st + v[0] * v[2] * (1 - ct);
  R[2][1] = v[0] * st + v[1] * v[2] * (1 - ct);
  R[2][2] = ct + powf(v[2], 2.) * (1 - ct);

  a[0] = R[0][0] * offset[0] + R[0][1] * offset[1] + R[0][2] * offset[2];
  a[1] = R[1][0] * offset[0] + R[1][1] * offset[1] + R[1][2] * offset[2];
  a[2] = R[2][0] * offset[0] + R[2][1] * offset[1] + R[2][2] * offset[2];

  return 1;
}


Molecule *allocate_molecule(int num_atoms){
  Molecule *mol = malloc(sizeof(Molecule));
  mol->num_atoms = num_atoms;
  mol->names = (int *)malloc(sizeof(int) * num_atoms);
  mol->atom_nums =(int *)malloc(sizeof(int) * num_atoms);
  mol->coords = (double *)malloc(sizeof(double)* num_atoms * 3);
  return mol;
}

Molecule* to_mol(Sdf_Ctab_Atom *atoms, int num_atoms){
  int i;
  Molecule *a=NULL;

  int nl=0;
  Sdf_Ctab_Atom *atom;
  for(i=0; i<num_atoms; i++){
    atom = atoms+i;
    if(atom->symbol[0]== 'R' && atom->symbol[1] == '#') continue;
    nl++;
  }
  a = allocate_molecule(nl);
  a->num_atoms = nl;

  nl = 0;
  for(i=0; i<num_atoms; i++){
    atom = atoms+i;
    a->names[nl] = atom->global_id;
    if(atom->symbol[0]== 'R' && atom->symbol[1] == '#') continue;
    a->names[nl] = atom->global_id ;
    a->atom_nums[nl] =  atom->atom_num;
    memcpy(&a->coords[nl * 3], atom->pos, 3 * sizeof(double));
    //        printf("%d %d %f %f %f\n", a->names[nl], a->atom_nums[nl], atom->pos[0], atom->pos[1], atom->pos[2]);
    nl++;
  }

  return a;
}


void free_molecule(Molecule *mol){
  free(mol->coords);
  free(mol->names);
  free(mol->atom_nums);
  free(mol);
}


int check_site_overlap(Sdf_Ctab_Atom *atoms, int num_atoms, double pos[3], double cutoff_value){
  double r1[3];
  double dist2;
  double cutoff2= cutoff_value*cutoff_value;
  int i;
  Sdf_Ctab_Atom *atom;

   for (int i=0; i< num_atoms; i++){
    atom = atoms+i;
    if(atom->symbol[0]== 'R' && atom->symbol[1] == '#') continue;
    r1[0] = atom->pos[0] - pos[0];
    r1[1] = atom->pos[1] - pos[1];
    r1[2] = atom->pos[2] - pos[2];
    dist2= r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2];
    if(dist2 < cutoff2){   // too close
      //      printf(" %f %f \n", dist2, cutoff2);
      return 0;
    }
  }
  return 1;
}


int check_mol_overlap(Molecule *mol, Molecule *tmp_mol, double cutoff_scalar) {
  // number of atoms in a molecule
  int mn = mol->num_atoms;
  int cn = tmp_mol->num_atoms;
  int nn = mn + cn;
  double temp;
  int mol_idz = 0;
  double cutoff;
  double rr1[3], rr2[3];
  double dist;
  // Compute inv lattice
  double *radius = malloc(sizeof(double) * nn);

  int i=0;
  int j=0;
  for (i = 0; i < cn; i++) {
    radius[i] = get_pte_vdw(tmp_mol->atom_nums[i]);
  }
  for(i=0; i< mn; i++){
    radius[i+cn] = get_pte_vdw(mol->atom_nums[i]);
  }

   for (i = 0; i < mn; i++) {
    mol_idz = (i + cn);
    rr1[0] = mol->coords[i*3];
    rr1[1] = mol->coords[i*3+1];
    rr1[2] = mol->coords[i*3+2];
    for (j = 0; j < cn; j++) {
      rr2[0] = tmp_mol->coords[j * 3];
      rr2[1] = tmp_mol->coords[j * 3 + 1];
      rr2[2] = tmp_mol->coords[j * 3 + 2];
      cutoff = (radius[j]+radius[mol_idz])*cutoff_scalar;
      cutoff = cutoff * cutoff;
      dist=get_dist2(rr1, rr2);
      if(dist < cutoff){
        //                        printf("dist2 is %f < %f\n", dist, cutoff);
        free(radius);
        return 0;
      }
    }
   }
   free(radius);
   return 1;
}


void next_pos(double *prev, double *pos, double radius, uint64_t *rng){
  // Generate a random point on the sphere
  double dr, u, v, th, phi;
  double dv[3];
  dr = 0.8 + 0.4*rng_uniform(rng);
  dr = dr*radius;
  u = rng_uniform(rng);
  v = rng_uniform(rng);
  th = acos(2.0*v-1.0);
  phi = TWOPI*u;
  dv[0] = dr*sin(th)*cos(phi);
  dv[1] = dr*sin(th)*sin(phi);
  dv[2] = dr*cos(th);
  pos[0] = prev[0] + dv[0];
  pos[1] = prev[1] + dv[1];
  pos[2] = prev[2] + dv[2];
}



void print_mol2xyzfile(Molecule *mol, FILE *outfile) {
  static int counter = 1;

  char *fmt;
  int N = mol->num_atoms;
  fprintf(outfile, " %d\n", N);
  //fprintf(out_file,
  //        "Properties=species:S:1:pos:R:3 pbc=\"F F F\" mol_id=%d "
  //        "struct_number=%d \n",
  //        mol->mol_id, counter);
  fprintf(outfile," Lattice=\"15 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0\" Properties=species:S:1:pos:R:3 pbc=\"T T T\" \n");

  fmt = "%3s % 14.8f % 14.8f % 14.8f\n";
  for (int i = 0; i < N; i++) {
    fprintf(outfile, fmt, atno2sym(mol->atom_nums[i]), mol->coords[i * 3],
            mol->coords[i * 3 + 1], mol->coords[i * 3 + 2]);
  }

  fflush(outfile);
}


void read_cyc_data(Settings *set){
  FILE *fileptr;
  int i = 0;
  char buffer[LINE_SIZE + 1], *ptr, *ptr2;
  int num_pair;

  // find_number_of atoms
  fileptr = fopen(set->cyc_file, "r");
  // check if file exits
  if (!fileptr) {
    printf("***ERROR: no %s file \n", set->cyc_file);
    exit(EXIT_FAILURE);
  }

  if (!fgets(buffer, LINE_SIZE, fileptr)) error_exit("Unexpected EOF");
  if (sscanf(buffer, "%d", &num_pair) != 1)
    error_exit("Unable to read number of atoms and bonds");

  //  if (!fgets(buffer, LINE_SIZE, fileptr)) error_exit("Unexpected EOF");

  //int *bondinfo=NULL;
  //bondinfo =(int *) malloc(2*num_bonds*sizeof(int));
  set->cyc=(int *)malloc(2*num_pair*sizeof(int));
  for (i=0; i< num_pair; i++){
    if (!fgets(buffer, LINE_SIZE, fileptr)) error_exit("Unexpected EOF");
    ptr = buffer;
    while (isspace(*ptr)) ptr++;
    if (isdigit(*ptr)) {
      if (sscanf(ptr, "%d %d", &set->cyc[i*2], &set->cyc[i*2+1]) != 2){
        error_exit("Error parsing bond info line");
      }
      printf("%d %d \n", set->cyc[i*2], set->cyc[i*2+1]);
    }
  }
}


void read_sdf_v2000_better(RgroupRecord *rg, Sdf_MetaData *meta, Settings *set) {
  int i;
  FILE *fp =NULL;
  char buffer[LINE_SIZE + 1];
  char *ptr;
  char end_str[4];
  int tmp;
  int records_read=0;

  Sdf_Ctab *cc=NULL;
  Sdf_Ctab_Counter *counter=NULL;


  fp = fopen(set->sdf_file, "r");
  if(!fp){
    printf("***ERROR: error in %s file\n", set->sdf_file);
    exit(EXIT_FAILURE);
  }

  int nmodels = set->num_mols;

  rg->num_mol = nmodels; // number molecule
  rg->groupid = malloc_Mat2d_i(nmodels, MAX_CONNECTIONS);
  rg->total_num =0;
  rg->num_rgp = (int *)malloc(nmodels * sizeof(int));
  rg->rgp = malloc(nmodels * MAX_CONNECTIONS*sizeof(int));
  rg->atomid = malloc_Mat2d_i(nmodels, MAX_CONNECTIONS);

  meta->num_data = nmodels;

  meta->ctabs = malloc(nmodels*sizeof(Sdf_Ctab));
  meta->counts = malloc(nmodels*sizeof(Sdf_Ctab_Counter));

  // use to save finale without R group
  Sdf_MetaData *final_meta=NULL;
  final_meta = malloc(sizeof(Sdf_MetaData));
  final_meta->num_data = 0; // important
  final_meta->ctabs = malloc(nmodels*sizeof(Sdf_Ctab));
  final_meta->counts = malloc(nmodels*sizeof(Sdf_Ctab_Counter));


  for (int n = 0; n < nmodels; n++) {
    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");
    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");
    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");

    if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexped");

    cc = meta->ctabs+n;
    counter=meta->counts+n;

    records_read =sscanf(buffer, "%d %d %d %d %d %d %d %d %d %d %d %s\n",
                         &(counter->num_atoms), &(counter->num_bonds),
               &(counter->num_atom_lists), &tmp, /* for obsolete record */
               &(counter->flag_chiral), &(counter->num_stext_entries),
               &(counter->num_react_components), &(counter->num_reactants),
               &(counter->num_products), &(counter->num_intermediates),
               &(counter->num_additional_lines), counter->version);

    if (records_read != 10)
      error_exit("Error parsing sdf headline\n");

    counter->last_globalid = -1; // init
    cc->counter.num_atoms = counter->num_atoms;
    cc->counter.num_bonds = counter->num_bonds;

    int num_atoms = counter->num_atoms;
    int num_bonds = counter->num_bonds;
    printf(" index molec %d num_atoms is %d num_bonds is %d\n", n, num_atoms, num_bonds);

    cc->atoms = (Sdf_Ctab_Atom *)malloc(num_atoms*sizeof(Sdf_Ctab_Atom));
    cc->bonds = (Sdf_Ctab_Bond *)malloc(num_bonds*sizeof(Sdf_Ctab_Bond));

    int index_rg=0;
    Sdf_Ctab_Atom *atom=NULL;

    for (i = 0; i < num_atoms; i++) {
      if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF");
      atom = cc->atoms + i;
      ptr = buffer;
      while (isspace(*ptr)) ptr++;
      if (isdigit(*ptr)) {
        records_read =
            sscanf(ptr, "%lf %lf %lf %s %d %d %d %d %d %d %d %d %d %d %d %d",
                   &(atom->pos[0]), &(atom->pos[1]), &(atom->pos[2]), atom->symbol,
                   &(atom->mass_diff), &(atom->charge), &(atom->stereo_parity),
                   &(atom->hydrogen_count), &(atom->flag_stereo_care),
                   &(atom->valence), &(atom->h0_designator),
                   &(atom->type_react_component), &(atom->num_react_component),
                   &(atom->mapping_number), &(atom->flag_invertion),
                   &(atom->flag_exact_change));
        atom->global_id = -1;
        atom->atom_num = atsym2no(atom->symbol); //
        atom->mol_id = n; // biaoji

        if (records_read != 16) error_exit("Error parsing sdf ctab atom data");
        if (atom->symbol[0] == 'R' && atom->symbol[1] == '#') {
          rg->total_num++;
          index_rg++;
          rg->num_rgp[n] = index_rg;
        }
      }
    }
    printf("rg-num_rgp[%d] is %d\n", n, rg->num_rgp[n]);
    Sdf_Ctab_Bond *bond=NULL;

    for(i=0; i< num_bonds; i++){
      bond = cc->bonds+i;
      if (!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF");
      ptr = buffer;

      while(isspace(*ptr)) ptr++;
      if(isdigit(*ptr)){
        records_read = sscanf(ptr, "%d %d %d %d %d %d",
                              &(bond->fatom),
                              &(bond->satom),
                              &(bond->type),
                              &(bond->stereo),
                              &(bond->topology),
                              &(bond->reacting_center_status)
                              );
        if(records_read != 6) error_exit("Error parsing bond data");
      }
    }

    if(!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF");

    int num_rgp = rg->num_rgp[n]; // num rgroup in molecule
    int tmp_rgp[9] = {0,0,0,0,0,0,0,0,0};
    if(num_rgp == 1){
      records_read = sscanf(buffer, "M  RGP %d %d %d",
                            &(tmp_rgp[0]),
                            &(tmp_rgp[1]),
                            &(tmp_rgp[2]));
      if(records_read != 3) error_exit("Error parsing RGP number 1");

      if(tmp_rgp[0] != num_rgp) error_exit("Error RGroup number in num_rgp=1 \n");
      rg->atomid[n][0] = tmp_rgp[1];  // R's local index
      rg->groupid[n][0] = tmp_rgp[2]; // R's names
      printf("RGP is %d %d %d\n", num_rgp,rg->atomid[n][0], rg->groupid[n][0]);
    }else if(num_rgp == 2){
      records_read = sscanf(buffer, "M  RGP %d %d %d %d %d",
                            &(tmp_rgp[0]),
                            &(tmp_rgp[1]),
                            &(tmp_rgp[2]),
                            &(tmp_rgp[3]),
                            &(tmp_rgp[4]));
      if(records_read != 5) error_exit("Error parsing RGP number 2");

      if(tmp_rgp[0] != num_rgp) error_exit("Error RGroup number in num_rgp=2 \n");
      rg->atomid[n][0] = tmp_rgp[1];  // R's local index
      rg->groupid[n][0] = tmp_rgp[2]; // R's names
      rg->atomid[n][1] = tmp_rgp[3];  // R's local index
      rg->groupid[n][1] = tmp_rgp[4]; // R's names
      printf("RGP is %d %d %d %d %d \n", num_rgp, rg->atomid[n][0], rg->groupid[n][0], rg->atomid[n][1], rg->groupid[n][1]);
    }
    else if(num_rgp == 3){
       records_read = sscanf(buffer, "M  RGP %d %d %d %d %d %d %d",
                            &(tmp_rgp[0]),
                            &(tmp_rgp[1]),
                            &(tmp_rgp[2]),
                            &(tmp_rgp[3]),
                             &(tmp_rgp[4]),
                             &(tmp_rgp[5]),
                             &(tmp_rgp[6]));
      if(records_read != 7) error_exit("Error parsing RGP number 3");
      if(tmp_rgp[0] != num_rgp) error_exit("Error RGroup number in num_rgp=3 \n");
      rg->groupid[n][0] = tmp_rgp[2]; // R's names
      rg->groupid[n][1] = tmp_rgp[4];
      rg->groupid[n][2] = tmp_rgp[6];

      rg->atomid[n][0] = tmp_rgp[1]; // R's local index
      rg->atomid[n][1] = tmp_rgp[3];
      rg->atomid[n][2] = tmp_rgp[5];
      printf("RGP is %d %d %d %d %d %d %d \n", num_rgp,rg->atomid[n][0], rg->groupid[n][0],
             rg->atomid[n][1], rg->groupid[n][1], rg->atomid[n][2], rg->groupid[n][2]);
    }

    else if(num_rgp == 4){
       records_read = sscanf(buffer, "M  RGP %d %d %d %d %d %d %d %d %d",
                            &(tmp_rgp[0]),
                            &(tmp_rgp[1]),
                            &(tmp_rgp[2]),
                            &(tmp_rgp[3]),
                             &(tmp_rgp[4]),
                             &(tmp_rgp[5]),
                             &(tmp_rgp[6]),
                             &(tmp_rgp[7]),
                             &(tmp_rgp[8]));
      if(records_read != 9) error_exit("Error parsing RGP number 4");
      if(tmp_rgp[0] != num_rgp) error_exit("Error RGroup number in num_rgp=4 \n");
      rg->groupid[n][0] = tmp_rgp[2];  // Rgroup name is int format
      rg->groupid[n][1] = tmp_rgp[4];
      rg->groupid[n][2] = tmp_rgp[6];
      rg->groupid[n][3] = tmp_rgp[8];

      rg->atomid[n][0] = tmp_rgp[1];   //local atom index in mol
      rg->atomid[n][1] = tmp_rgp[3];
      rg->atomid[n][2] = tmp_rgp[5];
      rg->atomid[n][3] = tmp_rgp[7];

      printf("RGP is %d %d %d %d %d %d %d %d %d \n", num_rgp, rg->atomid[n][0], rg->groupid[n][0],
             rg->atomid[n][1], rg->groupid[n][1],
             rg->atomid[n][2], rg->groupid[n][2],
             rg->atomid[n][3], rg->groupid[n][3]);
    }
    else{
      printf("number RGP is %d  allowed number is %d not supported\n", num_rgp, MAX_CONNECTIONS);
      exit(0);
    }

    if(!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF"); // M  END
    if(!fgets(buffer, LINE_SIZE, fp)) error_exit("Unexpected EOF"); // $$$$

    records_read = sscanf(buffer, "%s", end_str);
    if(end_str[0] == '$'){
      printf("parser sdf finished\n");
    }
  }

  printf("finished reading sdf\n");

  int max_pair = (rg->num_mol-1)*2;
  int max_total_atoms=0;
  int max_total_bonds=0;
  int **pair_record=malloc_Mat2d_i(max_pair, 4);
  int num_remove_atoms = get_pair(pair_record,
                                  &max_total_atoms,
                                  &max_total_bonds,
                                  rg, meta->counts);

  int **reorder_pair=malloc_Mat2d_i(num_remove_atoms, 4);

  int is_right=gen_reorder_pair(pair_record, num_remove_atoms, set->cyc, set->num_cycle, rg, reorder_pair);

  if(is_right ==1) {
    printf("gen reorder_pair correctly\n");}

  printf("total_atoms is %d total_bond is %d will remove %d\n",
         max_total_atoms, max_total_bonds, num_remove_atoms);

  int **sindex=malloc_Mat2d_i(nmodels, max_total_atoms);

  Sdf_Ctab *tmp_ctab=malloc(sizeof(Sdf_Ctab));
  //Sdf_Ctab *final_ctab = final_meta->ctabs; // just a pointer

  for(int n=0; n<nmodels; n++){
    counter = meta->counts+n;
    int num_atoms = counter->num_atoms;
    for(i=0; i<num_atoms; i++){
      sindex[n][i] = -1;
    }
  }

  tmp_ctab->atoms = (Sdf_Ctab_Atom *)malloc(max_total_atoms*sizeof(Sdf_Ctab_Atom));
  tmp_ctab->bonds = (Sdf_Ctab_Bond *)malloc(max_total_bonds*sizeof(Sdf_Ctab_Bond));

  for(i=0; i<max_total_atoms; i++){
    tmp_ctab->atoms[i].global_id = -1;
  }

  int index=0;

  for(int n=0; n<nmodels; n++){
    counter = meta->counts+n;
    cc = meta->ctabs+n;
    int num_atoms = counter->num_atoms;
    Sdf_Ctab_Atom *atom =NULL;

    int nl =0;
    int start_index = index;
    for(i=0; i<num_atoms; i++){
      atom = cc->atoms + i;
      if(atom->symbol[0] == 'R' && atom->symbol[1]=='#') continue;
      //copy_ctab_atom(&(tmp_ctab->atoms[index]), atom);
      //tmp_ctab->atoms[index].global_id = index;
      atom->global_id = index; // notice this
      sindex[n][i] = index;
      index++;
      nl++;
    }
  } // finish atom list

  tmp_ctab->counter.num_atoms = index;

  //final_ctab->atoms = (Sdf_Ctab_Atom *)malloc(index*sizeof(Sdf_Ctab_Atom));
  //for(i=0; i<index; i++){
  //  copy_ctab_atom(&final_ctab->atoms[i], &tmp_ctab->atoms[i]);
  // }


  index =0;
  int find_rgroup_bond = 0;
  int ratom_id=0;
  for(int n=0; n<nmodels; n++){
    counter = meta->counts+n;
    cc = meta->ctabs+n;
    int num_bonds = counter->num_bonds;
    Sdf_Ctab_Bond *bond = NULL;
    int num_rgp = rg->num_rgp[n]; // num rgroup in molecule

    for(i=0; i<num_bonds; i++){
      bond = cc->bonds + i;
      find_rgroup_bond = 0;
      for(int nr=0; nr<num_rgp; nr++){
        ratom_id = rg->atomid[n][nr];
        if(ratom_id == bond->fatom || ratom_id == bond->satom){
          find_rgroup_bond=1;
          break; // break from nr loog
        }
      }

      if(find_rgroup_bond == 1) continue;
      else{
        // copy unchanged bond topology to tmp_ctab
        copy_ctab_atom(&(tmp_ctab->bonds[index].sca), &(cc->atoms[bond->fatom]));
        copy_ctab_atom(&(tmp_ctab->bonds[index].scb), &(cc->atoms[bond->satom]));
        tmp_ctab->bonds[index].fatom = cc->atoms[bond->fatom].global_id;
        tmp_ctab->bonds[index].satom = cc->atoms[bond->satom].global_id;
        tmp_ctab->bonds[index].type = bond->type;
        index++;
      }
    }
  }
  // build linker pair

  int pair[4];
  int atom1i, atom2i;
  Sdf_Ctab_Atom atom1, atom2;

  Sdf_Ctab *cc2=NULL;
  Sdf_Ctab_Counter *counter2=NULL;
  int find_rgroup_bond2=0;
  find_rgroup_bond = 0;

  int current_num_mol=0;
  int current_num_atoms=0;

  for(i=0; i<num_remove_atoms; i++){
    find_rgroup_bond =0;
    find_rgroup_bond2=0;

    pair[0] = pair_record[i][0]; // num_mol id
    pair[1] = pair_record[i][1] -1; // local atom id
    pair[2] = pair_record[i][2]; // mol2 id
    pair[3] = pair_record[i][3] - 1; // local atom id

    printf("bond connect index %d [%d:%d-%d:%d]\n", i, pair[0]+1, pair[1]+1, pair[2]+1, pair[3]+1);

    // mol1
    cc = meta->ctabs + pair[0];
    counter = meta->counts + pair[0];

    //mol2
    cc2 = meta->ctabs + pair[2];
    counter2 = meta->counts + pair[2];

    // now adjust cc2's position

    int num_bonds_1 = counter->num_bonds;
    int num_bonds_2 = counter2->num_bonds;

    Sdf_Ctab_Bond *bond=NULL;
    double trans[3];

    find_rgroup_bond = 0;

    for(int ii=0; ii<num_bonds_1; ii++){
      bond = cc->bonds + ii;
      if(bond->fatom == pair[1]){
        atom1i=bond->satom;
        find_rgroup_bond=1;
        break;
      }
      if(bond->satom == pair[1]){
        atom1i=bond->fatom;
        find_rgroup_bond=1;
        break;
      }
    }

    find_rgroup_bond2=0;
    for(int ii=0; ii<num_bonds_2; ii++){
      bond = cc2->bonds + ii;
      if(bond->fatom == pair[3]){
        atom2i=bond->satom;
        find_rgroup_bond2=1;
        break;
      }
      if(bond->satom == pair[3]){
        atom2i=bond->fatom;
        find_rgroup_bond2=1;
        break;
      }
    }

    if (find_rgroup_bond == 1 && find_rgroup_bond2 == 1) {
      copy_ctab_atom(&(tmp_ctab->bonds[index].sca), &(cc->atoms[atom1i]));
      copy_ctab_atom(&(tmp_ctab->bonds[index].scb), &(cc2->atoms[atom2i]));
      tmp_ctab->bonds[index].fatom = cc->atoms[atom1i].global_id;
      tmp_ctab->bonds[index].satom = cc2->atoms[atom2i].global_id;
      tmp_ctab->bonds[index].type = 1;
      // shift cc2 mol to close cc ?
      index++;
      // shift pair[2] molecule
    }
    // add atom_coords to
  }

  Molecule *tmp_mol = allocate_molecule(tmp_ctab->counter.num_atoms);

  FILE *out_file = fopen("tmp_mol.xyz", "w");
  FILE *out_file2 = fopen("mol2.extxyz", "w");

  int count_inner = 0;
  int verdict = 0;
  int new_group = 0;
  //  int fpair, fatom, spair, satom;

  int fatom, satom, l, j, k;
  int ratom1, ratom2;
  int fpair, spair;
  int seed;
  uint64_t rng[2];
  double ranpos[3];
  int num_bonds1, num_bonds2;

  int max_try = set->max_attempts;
  double cutoff_scalar = set->cutoff_scalar;
  double pos_cutoff = set->pos_cutoff;
  double site_cutoff = set->site_cutoff;

  srand((unsigned int)time(NULL));
  seed = rand();
  rng_seed(rng, seed);
  rng_next(rng);

  Sdf_Ctab_Bond *bond = NULL;

  find_rgroup_bond = 0;

  Sdf_Ctab_Atom *atom = NULL;

  tmp_mol->num_atoms = 0;
  current_num_atoms = 0;
  int num_atoms1 = 0;

  int num_atoms2 = 0;
  int find_pair = 0;
  int find_pos = 0;

  int *occupaied=malloc(sizeof(int)*num_remove_atoms);
  int *unoccupaied = malloc(sizeof(int)*num_remove_atoms);
  for(i=0; i<num_remove_atoms; i++) {
    occupaied[i]=-1;
    unoccupaied[i]=-1;
  }


  // int fpair, spair;
  for (i = 0; i < rg->num_mol; i++) {
    for (j = 0; j < rg->num_rgp[i]; j++) {
      fpair = rg->groupid[i][j];
      fatom = rg->atomid[i][j];

      printf("fatom is  %d rg->atomid[%d][%d]\n", fatom, i, j);
      //   find_pair = 0;

      // in mol1
      for (k = 0; k < rg->num_mol; k++) {
        if (k == i) break;
        //if (find_pair == 1) break;
        //  if(occupaied[find_pos] == 1) break;
        for (l = 0; l < rg->num_rgp[k]; l++) {
          spair = rg->groupid[k][l];
          satom = rg->atomid[k][l];
          printf("satom is  %d rg->atomid[%d][%d]\n", satom, k, l);
          int try_number = 0;

          if (fpair == spair) {
            //find_pair++;
            Sdf_Ctab *cc = meta->ctabs + i;
            Sdf_Ctab_Counter *counter = meta->counts + i;
            Sdf_Ctab *cc2 = meta->ctabs + k;
            Sdf_Ctab_Counter *counter2 = meta->counts + k;
            num_bonds1 = counter->num_bonds;
            num_bonds2 = counter2->num_bonds;
            num_atoms1 = counter->num_atoms;
            num_atoms2 = counter2->num_atoms;

            printf(
                " pair %d cc is %d cc2 is %d atoms1 %d atoms2 %d bonds1 %d bonds2 %d fatom %d "
                "satom %d\n",
                find_pair, i, k, num_atoms1,  num_atoms2, num_bonds1, num_bonds2, fatom, satom);

            for (int ii = 0; ii < num_bonds1; ii++) {
              bond = cc->bonds + ii;
              if (bond->fatom == fatom) {
                atom1i = bond->satom;
                ratom1 = bond->fatom;
                find_rgroup_bond = 1;
                break;
              }

              if (bond->satom == fatom) {
                atom1i = bond->fatom;
                ratom1 = bond->satom;
                find_rgroup_bond = 1;
                break;
              }
            }

            find_rgroup_bond2 = 0;

            for (int ii = 0; ii < num_bonds2; ii++) {
              bond = cc2->bonds + ii;
              if (bond->fatom == satom) {
                atom2i = bond->satom;
                ratom2 = bond->fatom;
                find_rgroup_bond2 = 1;
                break;
              }
              if (bond->satom == satom) {
                atom2i = bond->fatom;
                ratom2 = bond->satom;
                find_rgroup_bond2 = 1;
                break;
              }
            }

            printf("atom1i %d:%d-%d  atom2i %d:%d-%d \n", i, atom1i, ratom1, k,
                   atom2i, ratom2);
            //} else {

            if (tmp_mol->num_atoms == 0) {
              // first add atoms to tmp_mol
              for (int ii = 0; ii < num_atoms2; ii++) {
                atom = cc2->atoms + ii;
                if (atom->symbol[0] == 'R' && atom->symbol[1] == '#') continue;
                tmp_mol->names[current_num_atoms] = atom->global_id;
                tmp_mol->atom_nums[current_num_atoms] = atom->atom_num;
                memcpy(&tmp_mol->coords[current_num_atoms * 3], atom->pos,
                       3 * sizeof(double));
                current_num_atoms++;
              }
              tmp_mol->num_atoms = current_num_atoms;
              //              print_mol2xyzfile(tmp_mol, out_file);

            }
            // rotate loop

            printf("mol1 is %d:%d mol2 is %d:%d:%d \n", i, atom1i, k, atom2i, ratom2);

            double pos[3];
            pos[0] = cc2->atoms[atom2i-1].pos[0];
            pos[1] = cc2->atoms[atom2i-1].pos[1];
            pos[2] = cc2->atoms[atom2i-1].pos[2];

            double over1[3];
            while (try_number < max_try) {
              next_pos(pos, over1, pos_cutoff, rng);
              //printf(" cc2->atoms[%d] is %f %f %f over1 is %f %f %f\n", atom2i, pos[0], pos[1], pos[2], over1[0], over1[1], over1[2]);
              int isoverlap = check_site_overlap(
                  cc2->atoms, cc2->counter.num_atoms, over1, site_cutoff);
               if (isoverlap == 0) {
                 //printf("isoverlap\n");
                 continue;
               }

              sdf_gen_mol(cc->atoms, cc->counter.num_atoms, rng);

              over1[0] = over1[0] - cc->atoms[atom1i - 1].pos[0];
              over1[1] = over1[1] - cc->atoms[atom1i - 1].pos[1];
              over1[2] = over1[2] - cc->atoms[atom1i - 1].pos[2];

              for (int ii = 0; ii < cc->counter.num_atoms; ii++) {
                cc->atoms[ii].pos[0] += over1[0];
                cc->atoms[ii].pos[1] += over1[1];
                cc->atoms[ii].pos[2] += over1[2];
              }

              //printf(" cc->atoms will shift %f %f %f \n",over1[0], over1[1], over1[2]);

              Molecule *mol2 = NULL;
              mol2 = to_mol(cc->atoms, cc->counter.num_atoms);
              verdict = check_mol_overlap(mol2, tmp_mol, cutoff_scalar);
              if (verdict == 1) {
                printf("%d cc %d:%d will add to tmp_mol\n",find_pos, i, atom1i);
                occupaied[find_pos]=1;
                unoccupaied[find_pair]=-1;
                find_pair++;
                find_pos++;

                for (int ii = 0; ii < mol2->num_atoms; ii++) {
                  tmp_mol->names[current_num_atoms] = mol2->names[ii];
                  tmp_mol->atom_nums[current_num_atoms] = mol2->atom_nums[ii];
                  memcpy(&tmp_mol->coords[current_num_atoms * 3],
                         &mol2->coords[ii * 3], 3 * sizeof(double));
                  current_num_atoms++;
                }
                free_molecule(mol2);
                tmp_mol->num_atoms = current_num_atoms;
                //print_mol2xyzfile(tmp_mol, out_file);
                break;
              }
              try_number++;
              free_molecule(mol2);
            }
            if (try_number >= max_try) {
              printf("failed to generate s solid positions, l is %d\n", l);
              unoccupaied[find_pair]=1;
              l--;

            } else {
              printf("generate s solide position\n");
              break;
            }
          }
        }
      }
    }
  }
  printf("findpos is %d\n", find_pair);
  // check cycle?

  print_mol2xyzfile(tmp_mol, out_file);
  free(occupaied);
  free(unoccupaied);
  free_Mat2d_i(reorder_pair, num_remove_atoms);
}

int main(int argc, const char **argv){
  int n;

  FILE *sdf_file;
  Settings *set;
  set = (Settings *)malloc(sizeof(Settings));
  init_settings(set);
  read_settings(set);
  read_cyc_data(set);

  RgroupRecord *rg = (RgroupRecord *)malloc(sizeof(RgroupRecord));

  Sdf_MetaData *meta =(Sdf_MetaData *)malloc(sizeof(Sdf_MetaData));
  read_sdf_v2000_better(rg, meta, set);

  free_rgroup(rg);
  free_meta(meta);
  free_settings(set);
}
