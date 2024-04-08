

void Kronecker_Product(double *C,  const double *A, int nrows, int ncols,
                                                const double *B, int mrows, int mcols)
{
   int ccols, i, j, k, l;
   int block_increment;
   const double *pB;
   double *pC, *p_C;

   ccols = ncols * mcols;
   block_increment = mrows * ccols;
   for (i = 0; i < nrows; C += block_increment, i++)
      for (p_C = C, j = 0; j < ncols; p_C += mcols, A++, j++)
         for (pC = p_C, pB = B, k = 0; k < mrows; pC += ccols, k++)
            for (l = 0; l < mcols; pB++, l++) *(pC+l) = *A * *pB;

}


void Kronecker_Product_complex(complex *C, const complex *A, int nrows, int ncols,
                                               const complex *B, int mrows, int mcols)
{
   int ccols, i, j, k, l;
   int block_increment;
   const complex *pB;
   complex *pC, *p_C;

   ccols = ncols * mcols;
   block_increment = mrows * ccols;
   for (i = 0; i < nrows; C += block_increment, i++)
      for (p_C = C, j = 0; j < ncols; p_C += mcols, A++, j++)
         for (pC = p_C, pB = B, k = 0; k < mrows; pC += ccols, k++)
            for (l = 0; l < mcols; pB++, l++) *(pC+l) = *A * *pB;

}

