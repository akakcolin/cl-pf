#ifndef SK_DATA_H
/* ========================================================================
   $File: $
   $Date: $
   $Revision: $
   $Creator: Zenghui Liu $
   $Notice: (C) Copyright 2018. $
   ========================================================================
*/
#define SK_DATA_H
#include"common.h"
#include"gtb_io.h"

void alloc_OSlakoCont(OSlakoCont *a);

// parse fortran formated filebuffer
int parse_sk_data_line(char fbuffer[], double *result);

int OldSKData_readFromFile(TOldSKData *skData,
                           TRepSplineIn *repSplineIn,
                           TRepPolyIn *repPolyIn,
                           const char *filename,
                           const int homo,
                           const int iSp1,
                           const int iSp2);


void Alloc_TSKData(TSKData *sk);
void Free_TSKData(TSKData *sk);

void print_TOldSKData(TOldSKData *a);
void print_TRepPolyIn(TRepPolyIn *a);
void print_TRepSplineIn(TRepSplineIn *a);
void print_OSlakoCont(OSlakoCont *a);

void free_TOldSKData(TOldSKData *a);
void free_TRepPolyIn(TRepPolyIn *a);
void free_TRepSplineIn(TRepSplineIn *a);
void free_OSlakoCont(OSlakoCont *a);


#endif
