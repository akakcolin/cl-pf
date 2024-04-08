#include <stdlib.h>
#include <math.h>
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

double innerProduct(double *A, // rotate matrix,
                   const double *coords_a_czero,
                   const double *coords_b_czero,
                   const int num_atoms) {
  double x1, x2, y1, y2, z1, z2;
  double G1 = 0.0f, G2 = 0.0f;
  A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0f;

  int total_number_of_coordinates = 3 * num_atoms;
  for (int i = 0; i < total_number_of_coordinates; i += 3) {
    x1 = coords_a_czero[i];
    y1 = coords_a_czero[i + 1];
    z1 = coords_a_czero[i + 2];

    x2 = coords_b_czero[i];
    y2 = coords_b_czero[i + 1];
    z2 = coords_b_czero[i + 2];

    G1 += x1 * x1 + y1 * y1 + z1 * z1;
    G2 += x2 * x2 + y2 * y2 + z2 * z2;

    A[0] += (x1 * x2);
    A[1] += (x1 * y2);
    A[2] += (x1 * z2);

    A[3] += (y1 * x2);
    A[4] += (y1 * y2);
    A[5] += (y1 * z2);

    A[6] += (z1 * x2);
    A[7] += (z1 * y2);
    A[8] += (z1 * z2);
  }
  return (G1 + G2) * 0.5;
}

double calcRMSDOfTwoConformationsWithQCPMethod(double *A, const double E0,
                                              const int num_atoms,
                                              double *rot_matrix) {
  double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
  double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2, SyzSzymSyySzz2,
      Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2, SxzpSzx, SyzpSzy, SxypSyx,
      SyzmSzy, SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
  double C[4];
  double mxEigenV;
  double b, a, delta, x2;
  double oldg = 0.0f;
  double evalprec = 1e-11f;

  Sxx = A[0];
  Sxy = A[1];
  Sxz = A[2];
  Syx = A[3];
  Syy = A[4];
  Syz = A[5];
  Szx = A[6];
  Szy = A[7];
  Szz = A[8];

  Sxx2 = Sxx * Sxx;
  Syy2 = Syy * Syy;
  Szz2 = Szz * Szz;

  Sxy2 = Sxy * Sxy;
  Syz2 = Syz * Syz;
  Sxz2 = Sxz * Sxz;

  Syx2 = Syx * Syx;
  Szy2 = Szy * Szy;
  Szx2 = Szx * Szx;

  SyzSzymSyySzz2 = 2.0f * (Syz * Szy - Syy * Szz);
  Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

  C[2] = -2.0f * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
  C[1] = 8.0f * (Sxx * Syz * Szy + Syy * Szx * Sxz + Szz * Sxy * Syx -
                 Sxx * Syy * Szz - Syz * Szx * Sxy - Szy * Syx * Sxz);

  SxzpSzx = Sxz + Szx;
  SyzpSzy = Syz + Szy;
  SxypSyx = Sxy + Syx;
  SyzmSzy = Syz - Szy;
  SxzmSzx = Sxz - Szx;
  SxymSyx = Sxy - Syx;
  SxxpSyy = Sxx + Syy;
  SxxmSyy = Sxx - Syy;
  Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

  C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2 +
         (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) *
             (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2) +
         (-(SxzpSzx) * (SyzmSzy) + (SxymSyx) * (SxxmSyy - Szz)) *
             (-(SxzmSzx) * (SyzpSzy) + (SxymSyx) * (SxxmSyy + Szz)) +
         (-(SxzpSzx) * (SyzpSzy) - (SxypSyx) * (SxxpSyy - Szz)) *
             (-(SxzmSzx) * (SyzmSzy) - (SxypSyx) * (SxxpSyy + Szz)) +
         (+(SxypSyx) * (SyzpSzy) + (SxzpSzx) * (SxxmSyy + Szz)) *
             (-(SxymSyx) * (SyzmSzy) + (SxzpSzx) * (SxxpSyy + Szz)) +
         (+(SxypSyx) * (SyzmSzy) + (SxzmSzx) * (SxxmSyy - Szz)) *
             (-(SxymSyx) * (SyzpSzy) + (SxzmSzx) * (SxxpSyy - Szz));
  mxEigenV = E0;
  for (int i = 0; i < 50; ++i) {
    oldg = mxEigenV;
    x2 = mxEigenV * mxEigenV;
    b = (x2 + C[2]) * mxEigenV;
    a = b + C[1];
    delta = ((a * mxEigenV + C[0]) / (2.0f * x2 * mxEigenV + b + a));
    mxEigenV -= delta;
    if (fabs(mxEigenV - oldg) < fabs(evalprec * mxEigenV))
      break;
  }

  if (rot_matrix != NULL) {
    double a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, a41, a42,
        a43, a44;
    double a3344_4334, a3244_4234, a3243_4233, a3143_4133, a3144_4134,
        a3142_4132;
    double q1, q2, q3, q4, normq;
    double evecprec = 1e-6;
    double a2, x2, y2, z2;
    double xy, az, zx, ay, yz, ax;
    a11 = SxxpSyy + Szz - mxEigenV;
    a12 = SyzmSzy;
    a13 = -SxzmSzx;
    a14 = SxymSyx;
    a21 = SyzmSzy;
    a22 = SxxmSyy - Szz - mxEigenV;
    a23 = SxypSyx;
    a24 = SxzpSzx;
    a31 = a13;
    a32 = a23;
    a33 = Syy - Sxx - Szz - mxEigenV;
    a34 = SyzpSzy;
    a41 = a14;
    a42 = a24;
    a43 = a34;
    a44 = Szz - SxxpSyy - mxEigenV;
    a3344_4334 = a33 * a44 - a43 * a34;
    a3244_4234 = a32 * a44 - a42 * a34;
    a3243_4233 = a32 * a43 - a42 * a33;
    a3143_4133 = a31 * a43 - a41 * a33;
    a3144_4134 = a31 * a44 - a41 * a34;
    a3142_4132 = a31 * a42 - a41 * a32;
    q1 = a22 * a3344_4334 - a23 * a3244_4234 + a24 * a3243_4233;
    q2 = -a21 * a3344_4334 + a23 * a3144_4134 - a24 * a3143_4133;
    q3 = a21 * a3244_4234 - a22 * a3144_4134 + a24 * a3142_4132;
    q4 = -a21 * a3243_4233 + a22 * a3143_4133 - a23 * a3142_4132;
    double qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;
    if (qsqr < evecprec) {
      q1 = a12 * a3344_4334 - a13 * a3244_4234 + a14 * a3243_4233;
      q2 = -a11 * a3344_4334 + a13 * a3144_4134 - a14 * a3143_4133;
      q3 = a11 * a3244_4234 - a12 * a3144_4134 + a14 * a3142_4132;
      q4 = -a11 * a3243_4233 + a12 * a3143_4133 - a13 * a3142_4132;
      qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;
      if (qsqr < evecprec) {
        double a1324_1423 = a13 * a24 - a14 * a23,
               a1224_1422 = a12 * a24 - a14 * a22;
        double a1223_1322 = a12 * a23 - a13 * a22,
               a1124_1421 = a11 * a24 - a14 * a21;
        double a1123_1321 = a11 * a23 - a13 * a21,
               a1122_1221 = a11 * a22 - a12 * a21;
        q1 = a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
        q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
        q3 = a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
        q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
        qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;
        if (qsqr < evecprec) {
          q1 = a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
          q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
          q3 = a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
          q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
          qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;
        }
      }
    }
    normq = sqrt(qsqr);
    q1 /= normq;
    q2 /= normq;
    q3 /= normq;
    q4 /= normq;

    a2 = q1 * q1;
    x2 = q2 * q2;
    y2 = q3 * q3;
    z2 = q4 * q4;

    xy = q2 * q3;
    az = q1 * q4;
    zx = q4 * q2;
    ay = q1 * q3;
    yz = q3 * q4;
    ax = q1 * q2;

    // zhuanzhi
    rot_matrix[0] = a2 + x2 - y2 - z2;
    rot_matrix[1] = 2 * (xy - az);
    rot_matrix[2] = 2 * (zx + ay);
    rot_matrix[3] = 2 * (xy + az);
    rot_matrix[4] = a2 - x2 + y2 - z2;
    rot_matrix[5] = 2 * (yz - ax);
    rot_matrix[6] = 2 * (zx - ay);
    rot_matrix[7] = 2 * (yz + ax);
    rot_matrix[8] = a2 - x2 - y2 + z2;
  }
  if(num_atoms == 1){
      return 0.0;
  }
  return sqrt(fabs(2.0f * (E0 - mxEigenV) / num_atoms));
}


double calcRMSDOfTwoConformations(const double *coords_a_czero,
                                 const double *coords_b_czero,
                                 const int num_atoms, double *rot_matrix) {

  double rmsd=1000;
  double A[9];
  double E0 = innerProduct(A, coords_a_czero, coords_b_czero, num_atoms);
  rmsd =calcRMSDOfTwoConformationsWithQCPMethod(A, E0, num_atoms, rot_matrix);
  return rmsd;
}

double mol_rmsd(const double *coords_a, const double *coords_b, const int num_atoms)
{
  int tn =0;
  int tj =0;
  double coords_a_czero[num_atoms * 3];
  double coords_b_czero[num_atoms * 3];
  double cent_a[3]= {0.0, 0.0, 0.0};
  double cent_b[3]= {0.0, 0.0, 0.0};
  double res_rmsd=100.0;
  double rot_matrix[9];

  for (int i=0; i< num_atoms; i++){
    tn = 3*i;
    cent_a[0] += coords_a[tn];
    cent_a[1] += coords_a[tn+1];
    cent_a[2] += coords_a[tn+2];

    cent_b[0] += coords_b[tn];
    cent_b[1] += coords_b[tn+1];
    cent_b[2] += coords_b[tn+2];
  }
  cent_a[0] /= num_atoms;
  cent_a[1] /= num_atoms;
  cent_a[2] /= num_atoms;

  cent_b[0] /= num_atoms;
  cent_b[1] /= num_atoms;
  cent_b[2] /= num_atoms;

  for (int i = 0; i < num_atoms; ++i) {
    tn = 3 * i;
    coords_a_czero[tn] = coords_a[tn] - cent_a[0];
    coords_a_czero[tn + 1] = coords_a[tn + 1] - cent_a[1];
    coords_a_czero[tn + 2] = coords_a[tn + 2] - cent_a[2];
  }

  for (int j = 0; j < num_atoms; ++j) {
    tj = 3 * j;
    coords_b_czero[tj] = coords_b[tj] - cent_b[0];
    coords_b_czero[tj + 1] = coords_b[tj + 1] - cent_b[1];
    coords_b_czero[tj + 2] = coords_b[tj + 2] - cent_b[2];
  }

  res_rmsd=calcRMSDOfTwoConformations(coords_a_czero, coords_b_czero, num_atoms, rot_matrix);
  return res_rmsd;
}

