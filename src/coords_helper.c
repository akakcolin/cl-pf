
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

inline double dist(const double a[3]){ return a[0]*a[0]+a[1]*a[1] + a[2]*a[2]; }

void matxmat33(const double A[3][3], const double B[3][3], double dot[3][3]) {
  dot[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
  dot[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
  dot[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];

  dot[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
  dot[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
  dot[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];

  dot[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
  dot[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
  dot[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];
}

void transposeMatrix(const double source_matrix[3][3], double matrix[3][3]) {
  int i, j;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      matrix[i][j] = source_matrix[j][i];
    }
  }
}

void get_cart_im(const double *fc, const double *lat, double *cart_im) {
  int i, j, k;
  for (j = 0; j < 3; ++j) {
    for (i = 0; i < 27; ++i) {
      cart_im[i * 3 + j] = 0.0;
    }
    for (k = 0; k < 3; ++k) {
      for (i = 0; i < 27; ++i) {
        cart_im[i * 3 + j] += fc[i * 3 + k] * lat[k * 3 + j];
      }
    }
  }
}

void dot_2d(const double *fc, const double *lat, const int fc_shape0,
            const int lat_shape1, const int fc_shape1, double *cart_fc) {
  int i, j, k;
  for (j = 0; j < lat_shape1; ++j) {
    for (i = 0; i < fc_shape0; ++i) {
      cart_fc[i * lat_shape1 + j] = 0.0;
    }
    for (k = 0; k < fc_shape1; ++k) {
      for (i = 0; i < fc_shape0; ++i) {
        cart_fc[i * lat_shape1 + j] +=
            fc[i * fc_shape1 + k] * lat[k * lat_shape1 + j];
      }
    }
  }
}

void dot_2d_mod(const double *fc, const double *lat, const int fc_shape0,
                const int lat_shape1, const int fc_shape1, double *cart_fc) {
  int i, j, k;
  double dmod;
  double integer;
  for (j = 0; j < lat_shape1; ++j) {
    for (i = 0; i < fc_shape0; ++i) {
      cart_fc[i * lat_shape1 + j] = 0.0;
    }
    for (k = 0; k < fc_shape1; ++k) {
      for (i = 0; i < fc_shape0; ++i) {
        dmod = modf(fc[i * fc_shape1 + k], &integer);
        cart_fc[i * lat_shape1 + j] += dmod * lat[k * lat_shape1 + j];
      }
    }
  }
}

// Given fractional coordinates in the lattice basis, returns corresponding
// fractional coordinates in the lll basis.
// lll_inverse is 3x3 fc is n*3
void get_lll_frac_coords(const double fc[][3], const double lll_inverse[3][3],
                         const int fc_dim0, double fc_lll[][3]) {
  // np.dot(fc, lll_inverse)
  int i, j, k;
  for (j = 0; j < 3; ++j) {
    for (i = 0; i < fc_dim0; ++i) {
      fc_lll[i][j] = 0.0;
    }
    for (k = 0; k < 3; ++k) {
      for (i = 0; i < fc_dim0; ++i) {
        fc_lll[i][j] += fc[i][k] * lll_inverse[k][j];
      }
    }
  }
}


int invert3x3(const double *src, double *dst) {
  double det;
  // Compute adjoint:

  dst[0] = +src[4] * src[8] - src[5] * src[7];
  dst[1] = -src[1] * src[8] + src[2] * src[7];
  dst[2] = +src[1] * src[5] - src[2] * src[4];
  dst[3] = -src[3] * src[8] + src[5] * src[6];
  dst[4] = +src[0] * src[8] - src[2] * src[6];
  dst[5] = -src[0] * src[5] + src[2] * src[3];
  dst[6] = +src[3] * src[7] - src[4] * src[6];
  dst[7] = -src[0] * src[7] + src[1] * src[6];
  dst[8] = +src[0] * src[4] - src[1] * src[3];

  // Compute determinant

  det = src[0] * dst[0] + src[1] * dst[3] + src[2] * dst[6];

  if(ISZERO(det)) return -1;
  // Multiply adjoint with reciprocal of determinant:
  det = 1.0f / det;

  dst[0] *= det;
  dst[1] *= det;
  dst[2] *= det;
  dst[3] *= det;
  dst[4] *= det;
  dst[5] *= det;
  dst[6] *= det;
  dst[7] *= det;
  dst[8] *= det;
  return 0;
}

int invert2x2( const double *matrix, double *result) {
  double det = matrix[0]*matrix[3] - matrix[1]*matrix[2];

  if (ISZERO(det))
      return -1;

  result[0] = matrix[3] / det;
  result[1] = -matrix[1] / det;
  result[2] = -matrix[2] / det;
  result[3] = matrix[0] / det;
  return 0;
}

int invert4x4( const double *matrix, double *result) {
  double *M = matrix;
  double t[12];
  double det;
  int i;

  t[0] = M[10] * M[15];
  t[1] = M[14] * M[11];
  t[2] = M[6] * M[15];
  t[3] = M[14] * M[7];
  t[4] = M[6] * M[11];
  t[5] = M[10] * M[7];
  t[6] = M[2] * M[15];
  t[7] = M[14] * M[3];
  t[8] = M[2] * M[11];
  t[9] = M[10] * M[3];
  t[10] = M[2] * M[7];
  t[11] = M[6] * M[3];

  result[0] = t[0]*M[5] + t[3]*M[9] + t[4]*M[13];
  result[0] -= t[1]*M[5] + t[2]*M[9] + t[5]*M[13];
  result[1] = t[1]*M[1] + t[6]*M[9] + t[9]*M[13];
  result[1] -= t[0]*M[1] + t[7]*M[9] + t[8]*M[13];
  result[2] = t[2]*M[1] + t[7]*M[5] + t[10]*M[13];
  result[2] -= t[3]*M[1] + t[6]*M[5] + t[11]*M[13];
  result[3] = t[5]*M[1] + t[8]*M[5] + t[11]*M[9];
  result[3] -= t[4]*M[1] + t[9]*M[5] + t[10]*M[9];
  result[4] = t[1]*M[4] + t[2]*M[8] + t[5]*M[12];
  result[4] -= t[0]*M[4] + t[3]*M[8] + t[4]*M[12];
  result[5] = t[0]*M[0] + t[7]*M[8] + t[8]*M[12];
  result[5] -= t[1]*M[0] + t[6]*M[8] + t[9]*M[12];
  result[6] = t[3]*M[0] + t[6]*M[4] + t[11]*M[12];
  result[6] -= t[2]*M[0] + t[7]*M[4] + t[10]*M[12];
  result[7] = t[4]*M[0] + t[9]*M[4] + t[10]*M[8];
  result[7] -= t[5]*M[0] + t[8]*M[4] + t[11]*M[8];

  t[0] = M[8]*M[13];
  t[1] = M[12]*M[9];
  t[2] = M[4]*M[13];
  t[3] = M[12]*M[5];
  t[4] = M[4]*M[9];
  t[5] = M[8]*M[5];
  t[6] = M[0]*M[13];
  t[7] = M[12]*M[1];
  t[8] = M[0]*M[9];
  t[9] = M[8]*M[1];
  t[10] = M[0]*M[5];
  t[11] = M[4]*M[1];

  result[8] = t[0]*M[7] + t[3]*M[11] + t[4]*M[15];
  result[8] -= t[1]*M[7] + t[2]*M[11] + t[5]*M[15];
  result[9] = t[1]*M[3] + t[6]*M[11] + t[9]*M[15];
  result[9] -= t[0]*M[3] + t[7]*M[11] + t[8]*M[15];
  result[10] = t[2]*M[3] + t[7]*M[7] + t[10]*M[15];
  result[10]-= t[3]*M[3] + t[6]*M[7] + t[11]*M[15];
  result[11] = t[5]*M[3] + t[8]*M[7] + t[11]*M[11];
  result[11]-= t[4]*M[3] + t[9]*M[7] + t[10]*M[11];
  result[12] = t[2]*M[10] + t[5]*M[14] + t[1]*M[6];
  result[12]-= t[4]*M[14] + t[0]*M[6] + t[3]*M[10];
  result[13] = t[8]*M[14] + t[0]*M[2] + t[7]*M[10];
  result[13]-= t[6]*M[10] + t[9]*M[14] + t[1]*M[2];
  result[14] = t[6]*M[6] + t[11]*M[14] + t[3]*M[2];
  result[14]-= t[10]*M[14] + t[2]*M[2] + t[7]*M[6];
  result[15] = t[10]*M[10] + t[4]*M[2] + t[9]*M[6];
  result[15]-= t[8]*M[6] + t[11]*M[10] + t[5]*M[2];

  det = M[0]*result[0] + M[4]*result[1] + M[8]*result[2] + M[12]*result[3];

  if (ISZERO(det)) return -1;

  det = 1.0 / det;
  for (i = 0; i < 16; i++)
      result[i] *= det;

  return 0;
}



static void _calc_distance_array( const double ref[][3],
                                  int numref,
                                  const double conf[][3],
                                  int numconf,
                                  double* distances) {
  int i, j;
  double dx[3];
  double rsq;
#pragma omp parallel for private(i, j, dx, rsq) shared(distances)
  for (i=0; i<numref; i++) {
    for (j=0; j<numconf; j++) {
      dx[0] = conf[j][0] - ref[i][0];
      dx[1] = conf[j][1] - ref[i][1];
      dx[2] = conf[j][2] - ref[i][2];

      rsq = (dx[0]*dx[0]) + (dx[1]*dx[1]) + (dx[2]*dx[2]);
      *(distances+i*numconf+j) = sqrt(rsq);
    }
  }
}


int find_minimum_image(double *cent_a,
                  const double latt_matrix[9],
                  const double inv_latt_matrix[9]) {
  double rel_frac_a[3];
  double images[3];
  double delta_frac[3];

  rel_frac_a[0] = cent_a[0] * inv_latt_matrix[0] + cent_a[1] * inv_latt_matrix[3] + cent_a[2] * inv_latt_matrix[6];
  rel_frac_a[1] = cent_a[0] * inv_latt_matrix[1] + cent_a[1] * inv_latt_matrix[4] + cent_a[2] * inv_latt_matrix[7];
  rel_frac_a[2] = cent_a[0] * inv_latt_matrix[2] + cent_a[1] * inv_latt_matrix[5] + cent_a[2] * inv_latt_matrix[8];

  images[0] = round(rel_frac_a[0]);
  images[1] = round(rel_frac_a[1]);
  images[2] = round(rel_frac_a[2]);

  delta_frac[0] = rel_frac_a[0] - images[0];
  delta_frac[1] = rel_frac_a[1] - images[1];
  delta_frac[2] = rel_frac_a[2] - images[2];

  if(fabs(delta_frac[0])<0.5 && fabs(delta_frac[1])<0.5 && fabs(delta_frac[2])<0.5)
  {
    return 0;
  }

  if(delta_frac[0]>0.5) { delta_frac[0] =  delta_frac[0] - 1; }
  if(delta_frac[0]<-0.5) { delta_frac[0] = delta_frac[0]+1; }
  if(delta_frac[1]>0.5) { delta_frac[1] = delta_frac[1]-1; }
  if(delta_frac[1]<-0.5) { delta_frac[1] = delta_frac[1]+1; }
  if(delta_frac[2]>0.5) { delta_frac[2] = delta_frac[2]-1; }
  if(delta_frac[2]<-0.5) { delta_frac[2] = delta_frac[2]+1; }
  cent_a[0] = delta_frac[0] * latt_matrix[0] + delta_frac[1] * latt_matrix[3] + delta_frac[2] * latt_matrix[6];
  cent_a[1] = delta_frac[0] * latt_matrix[1] + delta_frac[1] * latt_matrix[4] + delta_frac[2] * latt_matrix[7];
  cent_a[2] = delta_frac[0] * latt_matrix[2] + delta_frac[1] * latt_matrix[5] + delta_frac[2] * latt_matrix[8];
  return 1;
}

void get_images(const double tau[3], double a[27][3]) {
  int i, j, k;
  int index = 0;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        a[index][0] = tau[i];
        a[index][1] = tau[j];
        a[index][2] = tau[k];
        index += 1;
      }
    }
  }
}


// not perfact, best to use reduced coords and reduced lattice
int is_contact_pbc_frac(
        const double *frac_coord, // reduced_coord
        const int num_atoms,
        const double lattice[9], // reduced_lattice
        const double *diameter,
        const double cutoff_scalar) {
  int i, j;
  double dist=0.0;

  int k=0;

  double cutoff=5.0;
  // move atom to the box
  double rr1[3];
  double rr2[3];


  double images_view[78] = {
    -1.00000,  -1.00000,  -1.00000,
    -1.00000,  -1.00000,   0.00000,
    -1.00000,  -1.00000,   1.00000,
    -1.00000,   0.00000,  -1.00000,
    -1.00000,   0.00000,   0.00000,
    -1.00000,   0.00000,   1.00000,
    -1.00000,   1.00000,  -1.00000,
    -1.00000,   1.00000,   0.00000,
    -1.00000,   1.00000,   1.00000,
    0.00000,  -1.00000,  -1.00000,
    0.00000,  -1.00000,   0.00000,
    0.00000,  -1.00000,   1.00000,
    0.00000,   0.00000,  -1.00000,
    0.00000,   0.00000,   1.00000,
    0.00000,   1.00000,  -1.00000,
    0.00000,   1.00000,   0.00000,
    0.00000,   1.00000,   1.00000,
    1.00000,  -1.00000,  -1.00000,
    1.00000,  -1.00000,   0.00000,
    1.00000,  -1.00000,   1.00000,
    1.00000,   0.00000,  -1.00000,
    1.00000,   0.00000,   0.00000,
    1.00000,   0.00000,   1.00000,
    1.00000,   1.00000,  -1.00000,
    1.00000,   1.00000,   0.00000,
    1.00000,   1.00000,   1.00000
  };


  double cart_image[3];
  double delta_image[3];
  double delta_rr[3];
  for(i=0; i<num_atoms; ++i)
  {
    rr1[0] = frac_coord[i*3];
    rr1[1] = frac_coord[i*3+1];
    rr1[2] = frac_coord[i*3+2];

    for(j=i+1; j<num_atoms; ++j)
    {
      rr2[0] = frac_coord[j*3];
      rr2[1] = frac_coord[j*3+1];
      rr2[2] = frac_coord[j*3+2];

      cutoff =(diameter[i]+diameter[j])*cutoff_scalar;
      cutoff = cutoff*cutoff;
      delta_rr[0]= rr1[0]-rr2[0];
      delta_rr[1]= rr1[1]-rr2[1];
      delta_rr[2]= rr1[2]-rr2[2];
#pragma unroll
      for(k=0; k<26; ++k)
      {
        delta_image[0] = delta_rr[0] -images_view[k*3];
        delta_image[1] = delta_rr[1] -images_view[k*3+1];
        delta_image[2] = delta_rr[2] - images_view[k*3+2];

        cart_image[0] = delta_image[0] * lattice[0] + delta_image[1] * lattice[3] + delta_image[2] * lattice[6];
        cart_image[1] = delta_image[0] * lattice[1] + delta_image[1] * lattice[4] + delta_image[2] * lattice[7];
        cart_image[2] = delta_image[0] * lattice[2] + delta_image[1] * lattice[5] + delta_image[2] * lattice[8];

        dist = cart_image[0]*cart_image[0] + cart_image[1]*cart_image[1] + cart_image[2]*cart_image[2];
        if(dist<4 && dist<cutoff) {
          return 1;
        }
      }
    }
  }
  return 0;
}

int is_contact_pbc(const double cart_coord[][3],
                   const int num_atoms,
                   const double lattice[9],
                   const double inverse_lattice[9],
                   const double cutoff)
{
  int i, j;
  double dist;
  double cutoff2 = cutoff*cutoff; // covlent_radis*1,4
  double rr[3];
  int k=0;
  // move atom to the box

  for(i=0; i<num_atoms; ++i)
  {
    for(j=i+1; j<num_atoms; ++j)
    {
      rr[0] = cart_coord[i][0] - cart_coord[j][0];
      rr[1] = cart_coord[i][1] - cart_coord[j][1];
      rr[2] = cart_coord[i][2] - cart_coord[j][2];
      k=find_minimum_image(rr, lattice, inverse_lattice);
      dist = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2];
      if(dist<cutoff2) { return 1; }
    }
  }
  return 0;
}

void pbc_shortest_vectors(const double fcoords1[][3],
                          const double fcoords2[][3],
                          const double lll_matrix[3][3],
                          const double lll_inverse[3][3],
                          const int fc1_dim0,
                          const int fc2_dim0,
                          const int mask[fc1_dim0][fc2_dim0],
                          const int has_mask,
                          const double lll_frac_tol[3],
                          const int has_ftol,
                          double *vectors,
                          double *ds) {
  int i, j, k, l, bestk;
  double images_view[27][3];
  double tau[3] = {-1, 0, 1};
  get_images(tau, images_view);
  double(*fc1_lll)[3] = (double(*)[3])malloc(3 * fc1_dim0 * sizeof(double));
  double(*fc2_lll)[3] = (double(*)[3])malloc(3 * fc2_dim0 * sizeof(double));
  // double lll_inverse[3][3];
  //  invert3x3(lll_mapping, lll_inverse);
  get_lll_frac_coords(fcoords1, lll_inverse, fc1_dim0, fc1_lll);
  get_lll_frac_coords(fcoords2, lll_inverse, fc2_dim0, fc2_lll);

  double(*cart_f1)[3] = (double(*)[3])malloc(3 * fc1_dim0 * sizeof(double));
  double(*cart_f2)[3] = (double(*)[3])malloc(3 * fc2_dim0 * sizeof(double));
  double(*cart_im)[3] = (double(*)[3])malloc(27 * 3 * sizeof(double));
  //int has_ftol = 0;

  dot_2d_mod(fc1_lll, lll_matrix, fc1_dim0, 3, 3, cart_f1);
  dot_2d_mod(fc2_lll, lll_matrix, fc2_dim0, 3, 3, cart_f2);
  get_cart_im(images_view, lll_matrix, cart_im);
  //printf("cart_f1 is %9.5f %9.5f %9.5f\n", cart_f1[0][0], cart_f1[0][1], cart_f1[0][2]);
  //printf("cart_f2 is %9.5f %9.5f %9.5f\n", cart_f2[0][0], cart_f2[0][1], cart_f2[0][2]);

  // double (*ds)[fc2_dim0] =
  // (double(*)[fc2_dim0])malloc(fc1_dim0*fc2_dim0*sizeof(double));
  double d, best,  da, db, dc, fdist;
  d=0.0;
  int within_frac = 1;
  best = 1000000;
  double pre_im[3];
  int index = 0;

  for (i = 0; i < fc1_dim0; ++i)
    for (j = 0; j < fc2_dim0; ++j) {
      within_frac = 0;
      if (!has_mask || mask[i][j] == 0) {
        within_frac = 1;
        if (has_ftol) {
          for (l = 0; l < 3; ++l) {
            fdist = fc2_lll[j][l] - fc1_lll[i][l];
            if (fabs(fdist - round(fdist)) > lll_frac_tol[l]) {
              within_frac = 0;
              break;
            }
          }
        }
        if (within_frac) {
          pre_im[0] = cart_f2[j][0] - cart_f1[i][0];
          pre_im[1] = cart_f2[j][1] - cart_f1[i][1];
          pre_im[2] = cart_f2[j][2] - cart_f1[i][2];
          best = 1000000;
#pragma unroll
          for (k = 0; k < 27; ++k) {
            da = pre_im[0] + cart_im[k][0];
            db = pre_im[1] + cart_im[k][1];
            dc = pre_im[2] + cart_im[k][2];
            d = da * da + db * db + dc * dc;
            if (d < best) {
              best = d;
              bestk = k;
            }
          }
          ds[i*fc2_dim0 +j] = best;
          index = (i * fc2_dim0 + j) * 3;
          vectors[index]     = pre_im[0] + cart_im[bestk][0];
          vectors[index + 1] = pre_im[1] + cart_im[bestk][1];
          vectors[index + 2] = pre_im[2] + cart_im[bestk][2];
        }
      }
      if (!within_frac) {
        ds[i*fc2_dim0 +j] = best;
        index = (i * fc2_dim0 + j) * 3;
        vectors[index] = 100000;
        vectors[index + 1] = 100000;
        vectors[index + 2] = 100000;
      }
    }
  free(cart_f1);
  free(cart_f2);
  free(cart_im);
  free(fc1_lll);
  free(fc2_lll);
}

double dot(const double A[3], const double B[3])
{
  return (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]);
}

double determinant(const double arr[3][3])
{
  return
    arr[0][0] * (arr[1][1]*arr[2][2] - arr[1][2]*arr[2][1]) +
    arr[0][1] * (arr[1][2]*arr[2][0] - arr[1][0]*arr[2][2]) +
    arr[0][2] * (arr[1][0]*arr[2][1] - arr[1][1]*arr[2][0]);
}

double abs_cap(double val, double max_abs_val)
{
   // returns the value with its absolute value capped at max_abs_val.
   // Particularly useful in passing values to trignometrix functions where
   // may result in argument >1 being passed in
   return MAX(MIN(val, max_abs_val), -max_abs_val);
}

void lattice_from_parameter(const double a,
                            const double b,
                            const double c,
                            const double alpha,
                            const double beta,
                            const double gamma,
                            double lattice[3][3])
{
  double alpha_r = degreesToRadians(alpha);
  double beta_r = degreesToRadians(beta);
  double gamma_r = degreesToRadians(gamma);
  double val = (cos(alpha_r)*cos(beta_r) - cos(gamma_r))/(sin(alpha_r)*sin(beta_r));
  double gamma_star = acos(val);
  lattice[0][0] = a*sin(beta_r);
  lattice[0][1] = 0.0;
  lattice[0][2] = a*cos(beta_r);
  lattice[1][0] = -b*sin(alpha_r)*cos(gamma_star);
  lattice[1][1] = b*sin(alpha_r)*sin(gamma_star);
  lattice[1][2] = b*cos(alpha_r);
  lattice[2][0] = 0.0;
  lattice[2][1] = 0.0;
  lattice[2][2] = c;
  int i=0, j=0;
  for(i=0;i <3; i++)
  {
  printf("\n");
    for(j=0; j<3; j++)
        printf("%9.5f \t", lattice[i][j]);
  }
}

int get_niggli_reduced_lattice(const double arr[3][3], double out[3][3])
{
  int result =0;
  const int iterations =100;
  double a[3];
  double b[3];
  double c[3];
  out[0][0]=0;
  out[0][1]=0;
  out[0][2]=0;
  out[1][0]=0;
  out[1][1]=0;
  out[1][2]=0;
  out[2][0]=0;
  out[2][1]=0;
  out[2][2]=0;
  a[0]=arr[0][0]; a[1]=arr[0][1]; a[2]=arr[0][2];
  b[0]=arr[1][0]; b[1]=arr[1][1]; b[2]=arr[1][2];
  c[0]=arr[2][0]; c[1]=arr[2][1]; c[2]=arr[2][2];

  const double origVolume = determinant(arr);
  const double tol = 0.00001* pow(origVolume, 1/3);


  int count=0;
  double G[3][3];
  G[0][0] = dot(a, a); G[0][1] = dot(a, b); G[0][2] = dot(a, c);
  G[1][0] = dot(a, b); G[1][1] = dot(b, b); G[1][2] = dot(b, c);
  G[2][0] = dot(a, c); G[2][1] = dot(c, b); G[2][2] = dot(c, c);
  int i,j;
  for(i=0;i <3; i++)
  {
  printf("\n");
    for(j=0; j<3; j++)
        printf("%9.5f \t", G[i][j]);
  }

  double A, B, C, E, N, Y;
  const double M1[3][3] ={{0.0, -1.0, 0.0},
                        {-1.0, 0.0, 0.0},
                        {0.0, 0.0, -1.0}};

  const double M1T[3][3]={{0.0, -1.0, 0.0},
                        {-1.0, 0.0, 0.0},
                        {0.0, 0.0, -1.0}};

  const double M2[3][3] ={{-1.0, 0.0, 0.0},
                        {0.0, 0.0, -1.0},
                        {0.0, -1.0, 0.0}};

  const double M2T[3][3] ={{-1.0, 0.0, 0.0},
                         {0.0, 0.0, -1.0},
                         {0.0, -1.0, 0.0}};
  double M[3][3];
  double tM[3][3];
  double MT[3][3];

  for(count=0; count<1000; ++count)
  {
    A=G[0][0];
    B=G[1][1];
    C=G[2][2];
    E=2*G[1][2];
    N=2*G[0][2];
    Y=2*G[0][1];
    if((A > (B+tol)) || ((fabs(A-B)<tol) && (fabs(E)> fabs(N)+tol)))
    {
      //A1
      matxmat33(G, M1, M);
      matxmat33(M1T, M, G);
    }
    if(B>(C+tol) || ((fabs(B-C)<tol) && (fabs(N)> fabs(Y)+tol)))
    {
       matxmat33(G, M2, M);
       matxmat33(M2T, M, G);
       continue;
    }

    int i, j, k;
    i=(E>0?-1:1);
    j=(N>0?-1:1);
    k=(Y>0?-1:1);
    if(i*j*k == -1)
    {
      if(fabs(Y)<tol)
      {
        k=-1;
      }
      else if(fabs(N)<tol)
      {
        j=-1;
      }
      else if(fabs(E)<tol)
      {
        i=-1;
      }
      else{
      printf("Error: unassigned and i*j*k<0!\n");
      return result;
      }
    }

    M[0][0] = i; M[0][1] = 0; M[0][2] = 0;
    M[1][0] = 0; M[1][1] = j; M[1][2] = 0;
    M[2][0] = 0; M[2][1] = 0; M[2][2] = k;
    matxmat33(G, M, tM);
    matxmat33(M, tM, G);

    A=G[0][0]; B=G[1][1]; C=G[2][2];
    E=2*G[1][2]; N=2*G[0][2]; Y=2*G[0][1];
    // A5
    if((fabs(E)> B+tol) || ((fabs(E-B) <tol) && (2*N<(Y-tol))) || ((fabs(E+B)<tol)
&& Y+tol <0))
    {
      M[0][0] = 1; M[0][1] = 0; M[0][2] = 0;
      M[1][0] = 0; M[1][1] = 1; M[1][2] =(-1*E/fabs(E));
      M[2][0] = 0; M[2][1] = 0; M[2][2] = 1;

      MT[0][0] = 1; MT[0][1] = 0; MT[0][2] = 0;
      MT[1][0] = 0; MT[1][1] = 1; MT[1][2] = 0;
      MT[2][0] = 0; MT[2][1] = (-1*E/fabs(E)); MT[2][2] = 1;

      matxmat33(G, M, tM);
      matxmat33(MT, tM, G);
    }
    // A6
    else if((fabs(N)> A+tol) || ((fabs(A-N) <tol) && (2*E<(Y-tol))) || (fabs(A+N)<tol
&& Y < -tol))
    {
      M[0][0] = 1; M[0][1] = 0; M[0][2] = (-1*N/fabs(N));
      M[1][0] = 0; M[1][1] = 1; M[1][2] = 0;
      M[2][0] = 0; M[2][1] = 0; M[2][2] = 1;

      MT[0][0] = 1; MT[0][1] = 0; MT[0][2] = 0;
      MT[1][0] = 0; MT[1][1] = 1; MT[1][2] = 0;
      MT[2][0] = (-1*N/fabs(N)); MT[2][1] = 0; MT[2][2] = 1;

      matxmat33(G, M, tM);
      matxmat33(MT, tM, G);
    }
    // A7
    else if((fabs(Y)> A+tol) || ((fabs(A-Y) <tol) && (2*E<(N-tol))) || (fabs(A+Y)<tol
&& N < -tol))
    {
      M[0][0] = 1; M[0][1] =(-1*Y/fabs(Y)) ; M[0][2]=0;
      M[1][0] = 0; M[1][1] = 1; M[1][2] = 0;
      M[2][0] = 0; M[2][1] = 0; M[2][2] = 1;

      MT[0][0] = 1; MT[0][1] =0 ; MT[0][2]=0;
      MT[1][0] = -1*Y/fabs(Y); MT[1][1] = 1; MT[1][2] = 0;
      MT[2][0] = 0; MT[2][1] = 0; MT[2][2] = 1;
      matxmat33(G, M, tM);
      matxmat33(MT, tM, G);
    }
    // A8
    else if((E+N+Y+A+B<0) || ((fabs(E+N+Y+A+B)< tol) && tol<(Y+(A+N)*2)))
    {
      M[0][0] = 1; M[0][1] = 0; M[0][2] = 1;
      M[1][0] = 0; M[1][1] = 1; M[1][2] = 1;
      M[2][0] = 0; M[2][1] = 0; M[2][2] = 1;

      MT[0][0] = 1; MT[0][1] = 0; MT[0][2] = 0;
      MT[1][0] = 0; MT[1][1] = 1; MT[1][2] = 0;
      MT[2][0] = 1; MT[2][1] = 1; MT[2][2] = 1;

      matxmat33(G, M, tM);
      matxmat33(MT, tM, G);
    }
    else
    {
      break;
    }
    }
    if(count>999)
    {
    return result;
    }

  printf("G\n");
    A=G[0][0]; B=G[1][1]; C=G[2][2];
  E=2*G[1][2]; N=2*G[0][2]; Y=2*G[0][1];
  double la = sqrt(A);
  double lb = sqrt(B);
  double lc = sqrt(C);
  double alpha = acos(E/2/lb/lc)/PI*180;
  double beta= acos(N/2/la/lc)/PI*180;
  double gamma= acos(Y/2/la/lb)/PI*180;
  printf("--- %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",la,lb, lc, alpha, beta, gamma );
  lattice_from_parameter(la, lb, lc, alpha, beta, gamma, out);

  result=1;
  return result;
  //
}

/*
检测周期性晶体结构是否合理
  算法分2个部分：
  （1）不同非对称单元在空间的是否碰撞；
    大体参展pbc_shortest_vector的算法过程得到距离值来判断（这里不需要计算那个大矢量）；
    计算两组坐标在pbc条件下的最短vector的距离来判断是否；
    考虑晶体周围的最近邻的image（27个）
    将初始分数坐标转换到lll_matrix基矢下（含分数和笛卡尔）；
    27个image的坐标也跟随lll变换;

    对不同结构单元的两个原子 需要遍历周围image，并记录下最短距离值（笛卡尔距离）
    判断距离值是否大于共价键长（需要考虑原子类型）

    (2)相同非对称单元内的原子与image中的同一结构单元内的原子距离是否<共价键*scale
    这一过程比较简单，需要注意在lll_matrix处理后的坐标上计算，和考虑空间更多方向，
    测试发现考虑到{210}方向就很好了

*/
int detectDistOfSymmetricalAndPeriodic( const double matrix[3][3],
                                        const double lll_matrix[3][3],
                                        const double lll_inverse[3][3],
                                        const double cutoff_scalar,
                                        double *all_coords,
                                        const int dim0,
                                        const int dim1,
                                        const double *diameter,
                                        const int size_dia
                                        ) {

  int i, j, k, l;
  int ii, jj, kk, ll;
  int lll_index;
  double d, best, da, db, dc;

  double cutoff=2.0;

  double dmod;
  double integer;
  double temp_dist2;
  double pre_im[3];

  double *dis2_mtx=malloc(dim1*dim1*sizeof(double));

  double *fc_all_lll=malloc(dim0*dim1*3*sizeof(double));
  double *cart_all_lll=malloc(dim0*dim1*3*sizeof(double));

  double cart_tmp[3*dim1];
  double fc_tmp[3*dim1];
  double *cart_f1=NULL;
  double *cart_f2=NULL;


  double images_view[81] = {
    -1.00000,   0.00000,   0.00000,
     0.00000,  -1.00000,   0.00000,
     0.00000,   0.00000,  -1.00000,
    -1.00000,  -1.00000,   0.00000,
    -1.00000,   0.00000,  -1.00000,
     0.00000,  -1.00000,  -1.00000,
    -1.00000,   0.00000,   1.00000,
    -1.00000,   1.00000,   0.00000,
     1.00000,   0.00000,  -1.00000,
    -1.00000,  -1.00000,   1.00000,
    -1.00000,   1.00000,  -1.00000,
     1.00000,  -1.00000,  -1.00000,
    -1.00000,  -1.00000,  -1.00000,
    0.00000,   0.00000,   0.00000
    -1.00000,   1.00000,   1.00000,
    1.00000,  -1.00000,   0.00000,
    0.00000,   0.00000,   1.00000,
    1.00000,   0.00000,   0.00000,
    0.00000,   1.00000,   0.00000,
    0.00000,  -1.00000,   1.00000,
    0.00000,   1.00000,  -1.00000,
    0.00000,   1.00000,   1.00000,
    1.00000,   0.00000,   1.00000,
    1.00000,   1.00000,   0.00000,
    1.00000,   1.00000,   1.00000,
    1.00000,  -1.00000,   1.00000,
    1.00000,   1.00000,  -1.00000,
  };

  // condition one
  //double *cart_im = malloc(81*sizeof(double));
  double cart_im[42];


  for (jj = 0; jj < 3; jj++) {
    for (ii = 0; ii < 14; ii++) {
      cart_im[ii * 3 + jj] = 0.0;
    }
    for (kk = 0; kk < 3; kk++) {
      for (ii = 0; ii < 14; ii++) {
        cart_im[ii * 3 + jj] += images_view[ii * 3 + kk] * lll_matrix[kk][jj];
      }
    }
  }


  for(i=0; i<dim0;i++){
    double *coord1=&all_coords[i*dim1*3];

    memset(fc_tmp, 0, 3*dim1*sizeof(double));
    memset(cart_tmp, 0, 3*dim1*sizeof(double));

    for(jj=0;jj<3;jj++){
      for (kk = 0; kk < 3; kk++) {
        for (ii = 0; ii < dim1; ii++) {
          fc_tmp[ii*3+jj] += coord1[ii*3 + kk] * lll_inverse[kk][jj];
        }
      }
    }

    for (jj = 0; jj < 3; jj++) {
      for (kk = 0; kk < 3; kk++) {
        for (ii = 0; ii < dim1; ii++) {
          dmod = modf(fc_tmp[ii * 3+ kk], &integer);
          cart_tmp[ii * 3+ jj] += dmod * lll_matrix[kk][jj];
        }
      }
    }

    for(ii=0; ii<dim1;ii++){
      lll_index=i*dim1*3+ii*3;
      fc_all_lll[lll_index] = fc_tmp[ii*3];
      fc_all_lll[lll_index+1]=fc_tmp[ii*3+1];
      fc_all_lll[lll_index+2]=fc_tmp[ii*3+2];

      cart_all_lll[lll_index] = cart_tmp[ii*3];
      cart_all_lll[lll_index+1]=cart_tmp[ii*3+1];
      cart_all_lll[lll_index+2]=cart_tmp[ii*3+2];
    }
  }

  for(i=0; i<dim0-1; i++)
  {
    cart_f1 = &cart_all_lll[i*dim1*3];

    for(j=i+1; j<dim0; ++j)
    {
      cart_f2 =&cart_all_lll[j*dim1*3];
      d=0.0;
      best=1000000;

      for (ii = 0; ii < dim1; ++ii)
      {
        for (jj = 0; jj < dim1; ++jj) {
            pre_im[0] = cart_f2[jj*3] - cart_f1[ii*3];
            pre_im[1] = cart_f2[jj*3+1] - cart_f1[ii*3+1];
            pre_im[2] = cart_f2[jj*3+2] - cart_f1[ii*3+2];
            best = 1000000;
#pragma unroll
            for (kk = 0; kk < 14; ++kk) {
                da = pre_im[0] + cart_im[kk*3];
                db = pre_im[1] + cart_im[kk*3+1];
                dc = pre_im[2] + cart_im[kk*3+2];
                d = da * da + db * db + dc * dc;
                if (d < best) {  best = d; }
            }
            dis2_mtx[ii*dim1+jj] = best;
        }
      }

      for(kk=0; kk<dim1; ++kk)
      {
        for(ll=0; ll<dim1; ++ll)
        {
          temp_dist2 = dis2_mtx[kk*dim1+ll];
          if (temp_dist2<4.0)
          {
            cutoff=(diameter[kk]+diameter[ll])*cutoff_scalar;
            cutoff = cutoff*cutoff;
            if(temp_dist2 < cutoff){
              //collpsed
              //printf("collpsed from condition one\n");
              free(fc_all_lll);
              free(cart_all_lll);
              free(dis2_mtx);
              return 1;
            }
          }
        }
      }
    }
  }
  free(dis2_mtx);
  //free(cart_all_lll);


  // condition two
  // 这一部分考虑 同一个非对称单元在周围晶胞image的分布情况
  // 判断碰撞的标准是2个不同image的原子距离 共价键长＊scala值即可
  double images_view2[150] = {
    -1.00000,   0.00000,   0.00000,
     0.00000,  -1.00000,   0.00000,
     0.00000,   0.00000,  -1.00000,
    -1.00000,  -1.00000,   0.00000,
    -1.00000,   0.00000,  -1.00000,
     0.00000,  -1.00000,  -1.00000,
    -1.00000,   0.00000,   1.00000,
    -1.00000,   1.00000,   0.00000,
     1.00000,   0.00000,  -1.00000,
    -1.00000,  -1.00000,   1.00000,
    -1.00000,   1.00000,  -1.00000,
     1.00000,  -1.00000,  -1.00000,
    -1.00000,  -1.00000,  -1.00000,
    -1.00000,   1.00000,   1.00000,
    1.00000,  -1.00000,   0.00000,
    0.00000,   0.00000,   1.00000,
    1.00000,   0.00000,   0.00000,
    0.00000,   1.00000,   0.00000,
    0.00000,  -1.00000,   1.00000,
    0.00000,   1.00000,  -1.00000,
    0.00000,   1.00000,   1.00000,
    1.00000,   0.00000,   1.00000,
    1.00000,   1.00000,   0.00000,
    1.00000,   1.00000,   1.00000,
    1.00000,  -1.00000,   1.00000,
    1.00000,   1.00000,  -1.00000,
    0.0, 1.0, 2.0,
    0.0, 2.0, 1.0,
    1.0, 0.0, 2.0,
    1.0, 2.0, 0.0,
    2.0, 0.0, 1.0,
    2.0, 1.0, 0.0,
    0.0, -1.0, 2.0,
    0.0, 2.0, -1.0,
    -1.0, 0.0, 2.0,
    -1.0, 2.0, 0.0,
    -2.0, -1.0, 0.0,
    2.0, 0.0, -1.0,
    2.0, -1.0, 0.0,
    0.0, -2.0, 1.0,
    0.0, 1.0, -2.0,
    1.0, 0.0, -2.0,
    1.0, -2.0, 0.0,
    -2.0, 0.0, 1.0,
    -2.0, 1.0, 0.0,
    0.0, -1.0, -2.0,
    0.0, -2.0, -1.0,
    -1.0, 0.0, -2.0,
    -1.0, -2.0, 0.0,
    -2.0, 0.0, -1.0
  };


  double rr1[3];
  double delta_rr[3];
  double cart_image[3];
  double delta_image[3];
  double dist;
/*
  double cart_im2[150];

  for (jj = 0; jj < 3; jj++) {
    for (ii = 0; ii < 50; ii++) {
      cart_im2[ii * 3 + jj] = 0.0;
    }
    for (kk = 0; kk < 3; kk++) {
      for (ii = 0; ii < 50; ii++) {
        cart_im2[ii * 3 + jj] += images_view2[ii * 3 + kk] * lll_matrix[kk][jj];
      }
    }
  }

*/
  for(i=0; i<dim0; i++)
  {
    double *new_coord =&fc_all_lll[i*dim1*3];
    for(ii=0; ii<dim1-1; ii++)
    {
      rr1[0] = new_coord[ii*3];
      rr1[1] = new_coord[ii*3+1];
      rr1[2] = new_coord[ii*3+2];
      for(jj=ii+1; jj<dim1; jj++)
      {
        delta_rr[0]= rr1[0]-new_coord[jj*3];
        delta_rr[1]= rr1[1]-new_coord[jj*3+1];
        delta_rr[2]= rr1[2]-new_coord[jj*3+2];
        cutoff =(diameter[ii]+diameter[jj])*cutoff_scalar;
        cutoff = cutoff*cutoff;
#pragma unroll
        for(k=0; k<50; k++)
        {
          delta_image[0] = delta_rr[0] - images_view2[k*3];
          delta_image[1] = delta_rr[1] - images_view2[k*3+1];
          delta_image[2] = delta_rr[2] - images_view2[k*3+2];
          cart_image[0] = delta_image[0] * lll_matrix[0][0] + delta_image[1] * lll_matrix[1][0] + delta_image[2] * lll_matrix[2][0];
          cart_image[1] = delta_image[0] * lll_matrix[0][1] + delta_image[1] * lll_matrix[1][1] + delta_image[2] * lll_matrix[2][1];
          cart_image[2] = delta_image[0] * lll_matrix[0][2] + delta_image[1] * lll_matrix[1][2] + delta_image[2] * lll_matrix[2][2];
          dist = cart_image[0]*cart_image[0] + cart_image[1]*cart_image[1] + cart_image[2]*cart_image[2];
          //dist=delta_image[0]*delta_image[0]+delta_image[1]*delta_image[1]+delta_image[2]*delta_image[2];
          if(dist<cutoff) {
            // printf("dist is %9.5f\n", dist);
            free(fc_all_lll);
             return 1;
           }
        }
      }
    }
  }
  free(fc_all_lll);
  return 0;
}


int collision_detect(const double *coords1, const double *coords2, const int len_a, const int len_b) {
  int i, j, k;
  double cent_a[3]={0.0, 0.0, 0.0};
  double cent_b[3]={0.0, 0.0, 0.0};
  double radius_a=0.0;
  double radius_b=0.0;
  double radius_tmp=0.0;
  double *delta_a=malloc(len_a*3*sizeof(double));
  double *delta_b=malloc(len_b*3*sizeof(double));
  double delta_radius=0.0;

  for (i = 0; i < len_a; ++i) {
    k = i * 3;
    cent_a[0] += coords1[k];
    cent_a[1] += coords1[k+1];
    cent_a[2] += coords1[k+2];
  }
  cent_a[0] /= len_a;
  cent_a[1] /= len_a;
  cent_a[2] /= len_a;

  for(i=0; i< len_a; ++i)
  {
    k= i*3;
    delta_a[k]   = coords1[k] - cent_a[0];
    delta_a[k+1] = coords1[k+1] - cent_a[1];
    delta_a[k+2] = coords1[k+2] - cent_a[2];

    radius_tmp = dist(&delta_a[k]);

    if(radius_tmp > radius_a)
    {
      radius_a =  radius_tmp;
    }
  }

  for (i = 0; i < len_b; ++i) {
    k = i * 3;
    cent_b[0] += coords2[k];
    cent_b[1] += coords2[k+1];
    cent_b[2] += coords2[k+2];
  }
  cent_b[0] /= len_b;
  cent_b[1] /= len_b;
  cent_b[2] /= len_b;

  for(i=0; i< len_b; ++i)
  {
    k= i*3;
    delta_b[k]   = coords2[k] - cent_b[0];
    delta_b[k+1] = coords2[k+1] - cent_b[1];
    delta_b[k+2] = coords2[k+2] - cent_b[2];

    radius_tmp = dist(&delta_b[k]);
    if(radius_tmp > radius_b)
    {
      radius_b =  radius_tmp;
    }
  }

  free(delta_a);
  free(delta_b);

  radius_tmp = getDist2(cent_a, cent_b);
  delta_radius = sqrt(radius_tmp) - sqrt(radius_a)-sqrt(radius_b);
  if(delta_radius)
    return 0;

  return 1;
}

int find_empiric_bonds(const double *coords, const int num_atoms, const double *diameter, const double cutoff_scalar, int *is_bond)
{
  int i, j;
  double ij_dist=0.0;
  double eps=0.0001;
  double cutoff;

  for(i=0; i< num_atoms-1; ++i) {
    for(j=i+1; j<num_atoms; ++j) {
      ij_dist =getDist2(&coords[i*3], &coords[j*3]);
      cutoff =(diameter[i]+diameter[j])*cutoff_scalar; // 原子间距
      cutoff = cutoff*cutoff;
      if((ij_dist-cutoff)<eps) { is_bond[i*num_atoms+j] = 1; }
      else { is_bond[i*num_atoms+j] = 0; }
    }
  }
  return 1;
}

void get_neighbor_list(const double *protein_coords, const int num_protein_atoms, const double *ligand_coords, const int num_ligand_atoms, const double cutoff_distance, int *segment_id_atom, int *a2a_neighbors, double *a2a_distance)
{
  int i, j;
  double ij_dist=0.0;
  double cutoff=cutoff_distance * cutoff_distance;
  int index=0;
  for(i=0; i< num_protein_atoms; ++i) {
    for(j=0; j<num_ligand_atoms; ++j) {
      ij_dist =getDist2(&protein_coords[i*3], &ligand_coords[j*3]);
      if(ij_dist < cutoff)
      {
         segment_id_atom[index] = i;
         a2a_neighbors[index] = j;
         a2a_distance[index] = ij_dist;
         index+=1;
      }
    }
  }
}
