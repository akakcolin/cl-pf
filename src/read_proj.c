#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<assert.h>
#include<complex.h>

#define DEBUG 0

typedef struct{
    int ichem;
    int lmax;
    int num_same_atoms;
    int L;
    int cols;
    int rows;
    double *b;
}Block;

typedef struct{
    int J;   // index for irreps
    int JD;  // index for the diagonal elements of the irrep
    int num_subblock;
    int lmax;
    Block *subblock;
}ProjectMatrix;

typedef struct{
    double k_vector[3];
    int nel;
    int *nat;
    int *lmax;
    int matrix_order;
    int num_proj_matrix;
    ProjectMatrix *proj_matrix;
    double *dftbmat;
    //double complex *dmat;
}KV;

void printDmat(double *dmat, int matrix_order){
    for(int i=0; i<matrix_order; i++){
        printf("%d\n", i);
        for(int j=0;j<matrix_order;j++){
            printf("(%9.5f, %9.5f)\t", dmat[i*matrix_order*2 + 2*j], dmat[2*i*matrix_order+2*j+1]);
        }
    }
}

void printDCmat(double complex *dmat, int matrix_order){
    int ij=0;
    for(int i=0; i<matrix_order; i++){
        printf("%d\n", i);
        for(int j=0;j<matrix_order;j++){
            ij = i*matrix_order+j;
            printf("(%9.5f, %9.5f)\t", creal(dmat[ij]), cimag(dmat[ij]));
        }
    }
}

void double_array_to_complex(double *dmatrix, double complex *dmat, const int matrix_order)
{
    if(dmat==NULL){printf("dmat is null;\n");
    }
    int ij=0;
    int ji=0;
    for(int i=0;i<matrix_order;i++)
        for(int j=0;j<matrix_order;j++){
            ij=i*matrix_order + j;
            //ji=j*matrix_order + i;
            dmat[ij]= dmatrix[2*ij] + dmatrix[2*ij+1]* _Complex_I;
    }
}


KV * parse_symproj(const char *fname, const int num_kpoints,
                   const int dftb_matrix_order){
    int nel;
    int no_wave_vectors;

    int i,j,k;
    int *nat=NULL;
    int *lmax=NULL;
    int *shift_block=NULL;

    KV *kv=NULL;  // return the ptr

    int *nAllOrb=NULL;
    Block *bb=NULL;

    int matrix_order=0;
    int g=0;
    int wv=0;
    
    double k_vector[3];

    int nAllAtoms=0;
    
    int ind=0;
    int nOrbAtom=1;
    int *iSquare=NULL;
    
    int num_each_subblock=0;
    
    FILE *myPtr;
    myPtr = fopen(fname, "r");
    if(!myPtr) {
        printf("\nUnable to open file: > %s", fname);
        exit(-1);
    }

    fread(&no_wave_vectors, sizeof(int), 1, myPtr);
    fread(&nel, sizeof(int),1,myPtr);

    assert(no_wave_vectors ==1);
    //printf("no_wave_vectors is %d\n", no_wave_vectors);
    //printf("nel is %d\n", nel);

    nat=malloc(nel*sizeof(int));
    lmax=malloc(nel*sizeof(int));

    for(i=0; i<nel; i++){
        fread(&nat[i], sizeof(int),1,myPtr);
        //printf("nat of %d, is %d\n", i, nat[i]);
        nAllAtoms += nat[i];
    }

    nAllOrb = malloc(sizeof(int)*nAllAtoms);
    for(i=0; i<nel; i++){
        fread(&lmax[i], sizeof(int),1,myPtr);
        //printf("lmax of %d, is %d\n", i, lmax[i]);
    }

    shift_block =malloc((nel+1)*sizeof(int));

    for(i=0;i<nel;i++){
        g = 0;
        for(j=0;j<lmax[i]+1;j++){
            g += 2*j +1;
        }
        shift_block[i+1]=g*nat[i];
        matrix_order += nat[i]*g;
    }

    shift_block[0]=0;
    for(i=0;i<nel;i++){
        g=0;
        for(j=0;j<i;j++){
            shift_block[i]+= shift_block[j];
        }

    }
    
    kv = malloc(sizeof(KV)*no_wave_vectors);

    // init kv 
    for(i=0;i<no_wave_vectors; i++){
        // read k_vector
        kv[i].k_vector[0] =0.0;
        kv[i].k_vector[1] =0.0;
        kv[i].k_vector[2] =0.0;

        kv[i].nel =nel;
        kv[i].nat = malloc(sizeof(int)*nel);
        kv[i].lmax = malloc(sizeof(int)*nel);

        for(j=0;j<nel;j++){
            kv[i].nat[j] = nat[j];
            kv[i].lmax[j] = lmax[j];
        }
        
        kv[i].matrix_order = matrix_order;
        kv[i].dftbmat =malloc(sizeof(double)*matrix_order*matrix_order*2);
        //kv[i].dmat =malloc(sizeof(double complex)*matrix_order*matrix_order);

        kv[i].num_proj_matrix = 1;
        for(j=0;j<matrix_order*matrix_order*2;j++){
            kv[i].dftbmat[j] = 0.0;
        }
    }

    iSquare = malloc(nAllAtoms*sizeof(int));
    
    k=0;
    for(i=0; i<nel; i++){
        if(lmax[i]==0) nOrbAtom=1;
        if(lmax[i]==1) nOrbAtom=1+3;
        if(lmax[i]==2) nOrbAtom=1+3+5;
        nOrbAtom *=2;  // careful
        for(j=0; j<nat[i]; j++){
            iSquare[k] = ind;
            ind += nOrbAtom;
            k++;
        }
    }

    while(wv < no_wave_vectors){
        num_each_subblock=0;
        fread(&(kv[wv].k_vector[0]), sizeof(double), 1, myPtr);
        fread(&(kv[wv].k_vector[1]), sizeof(double), 1, myPtr);
        fread(&(kv[wv].k_vector[2]), sizeof(double), 1, myPtr);

        fread(&(kv[wv].num_proj_matrix), sizeof(int), 1, myPtr);
        fread(&num_each_subblock, sizeof(int), 1, myPtr);

        if(kv[wv].num_proj_matrix == 1){
            printf("dmat must be E matix, so just return an matrix_order * matrix_order\n");
            for(j=0;j<matrix_order;j++){
                kv[i].dftbmat[j*matrix_order*2 + j*2] = 1.0;
            }
            continue;
        }

        kv[wv].proj_matrix = malloc(sizeof(ProjectMatrix)*kv[wv].num_proj_matrix);

        for(i=0; i<kv[wv].num_proj_matrix; i++){
            for(j=0;j<nel;j++){ // diff with mat index, 0 is start
                kv[wv].proj_matrix[i].num_subblock = num_each_subblock;
                kv[wv].proj_matrix[i].subblock = malloc(sizeof(Block)*(num_each_subblock));

            }
        }

        for(i=0; i<kv[wv].num_proj_matrix; i++){
            int J=0;
            int JD=0;
            int num_subblock=0;
            int num_same_atoms=0;
            int ichem=0;
            int L=0;
            int itotal=0;
            int ndi=0;
            int cols=0;
            int rows=0;

            fread(&J, sizeof(int), 1, myPtr);
            fread(&JD, sizeof(int), 1, myPtr);
            
            
            kv[wv].proj_matrix[i].J = J;
            kv[wv].proj_matrix[i].JD = JD;

            num_subblock = kv[wv].proj_matrix[i].num_subblock;
            num_same_atoms = nat[ichem];
#if DEBUG
            printf("T-matrix for J = %d, JD = %d\n", J, JD);
            printf("num subblock is %d\n", num_subblock);
            printf("num same atoms in unit cell is %d\n", num_same_atoms);
#endif
            
            kv[wv].proj_matrix[i].subblock=malloc(sizeof(Block)*num_subblock);


            for(j=0; j<num_subblock; j++){
                // ichem L, num_same_atoms, cols, rows
                //printf("subblock index %d\n", j);
                fread(&ichem, sizeof(int), 1, myPtr);
                fread(&L, sizeof(int), 1, myPtr);
                fread(&cols, sizeof(int), 1, myPtr); // itotal
                fread(&rows, sizeof(int), 1, myPtr); // ndi

                num_same_atoms = nat[ichem-1];
                kv[wv].proj_matrix[i].subblock[j].num_same_atoms = num_same_atoms;
                // kv[wv].proj_matrix[i].subblock[j].subindex = j;
                kv[wv].proj_matrix[i].subblock[j].ichem = ichem;
                kv[wv].proj_matrix[i].subblock[j].lmax = lmax[ichem-1];
#if DEBUG
                printf("L is %d and lmax[ichem] is %d\n", L, lmax[ichem-1]);
                printf("ichem is %d\n", ichem);
                printf(" cols is %d rows is %d\n", cols, rows);
                printf("nat[ichem] is %d \n", nat[ichem-1]);
                printf(" cols*rows*2 is %d\n", cols*rows*2);
#endif
                assert(L <= lmax[ichem-1]);

                // point bb just for convenient
                bb=&(kv[wv].proj_matrix[i].subblock[j]);
                bb->b =malloc(sizeof(double)*2*cols*rows);

                bb->L = L;
                bb->cols = cols;
                bb->rows = rows;
                //bb->subindex = j;

                for(k=0; k<cols*rows*2; k++){
                    fread(&(bb->b[k]), sizeof(double), 1, myPtr);
                }
                // change to dftbmat
                // printf subblock
#if DEBUG
                for(k=0; k<cols; k++){
                    printf("Column %d \n", k);
                    for(int kj=0; kj<rows; kj++){
                        printf("(%9.5f %9.5f) ", bb->b[k*rows*2+kj*2], bb->b[k*rows*2+kj*2+1]);
                    }
                    printf("\n");
                }
#endif
            }
        }
    }
    fclose(myPtr);
    
    // free malloc
    
    free(nat);
    free(shift_block);
    free(lmax);
    return kv;
}


void parse_symproj_(const char *fname, const int num_kpoints, const int dftb_matrix_order, double complex *dmatrixs){
    int nel;
    int no_wave_vectors;

    int i,j,k;
    int *nat=NULL;
    int *lmax=NULL;
    int matrix_order=0;
    int g=0;
    int *shift_block=NULL;

    int wv=0;
    KV *kv=NULL;

    double k_vector[3];

    int nAllAtoms=0;
    int *nAllOrb=NULL;
    Block *bb=NULL;


    FILE *myPtr;
    myPtr = fopen(fname, "r");
    if(!myPtr) {
        printf("\nUnable to open file: > %s", fname);
        exit(-1);
    }

    fread(&no_wave_vectors, sizeof(int), 1, myPtr);
    fread(&nel, sizeof(int),1,myPtr);

    assert(no_wave_vectors ==1);
    //printf("no_wave_vectors is %d\n", no_wave_vectors);
    //printf("nel is %d\n", nel);

    nat=malloc(nel*sizeof(int));
    lmax=malloc(nel*sizeof(int));

    for(i=0; i<nel; i++){
        fread(&nat[i], sizeof(int),1,myPtr);
        //printf("nat of %d, is %d\n", i, nat[i]);
        nAllAtoms += nat[i];
    }

    nAllOrb = malloc(sizeof(int)*nAllAtoms);
    for(i=0; i<nel; i++){
        fread(&lmax[i], sizeof(int),1,myPtr);
        //printf("lmax of %d, is %d\n", i, lmax[i]);
    }

    shift_block =malloc((nel+1)*sizeof(int));

    for(i=0;i<nel;i++){
        g = 0;
        for(j=0;j<lmax[i]+1;j++){
            g += 2*j +1;
        }
        shift_block[i+1]=g*nat[i];
        matrix_order += nat[i]*g;
    }

    shift_block[0]=0;
    for(i=0;i<nel;i++){
        g=0;
        for(j=0;j<i;j++){
            shift_block[i]+= shift_block[j];
        }

    }


    //printf("dftb_matrix_order is %d\n", dftb_matrix_order);
    //assert(matrix_order == dftb_matrix_order);

    if(dmatrixs == NULL){
        printf("error: dmatrixs is NULL\n ");
        exit(-1);
    }

    kv = malloc(sizeof(KV)*no_wave_vectors) ;

    for(i=0;i<no_wave_vectors; i++){
        // read k_vector
        kv[i].k_vector[0] =0.0;
        kv[i].k_vector[1] =0.0;
        kv[i].k_vector[2] =0.0;

        kv[i].nel =nel;
        kv[i].nat = malloc(sizeof(int)*nel);
        kv[i].lmax = malloc(sizeof(int)*nel);

        for(j=0;j<nel;j++){
            kv[i].nat[j] = nat[j];
            kv[i].lmax[j] = lmax[j];
        }
        kv[i].matrix_order = matrix_order;
        kv[i].dftbmat =malloc(sizeof(double)*matrix_order*matrix_order*2);
        //kv[i].dmat =malloc(sizeof(double complex)*matrix_order*matrix_order);

        kv[i].num_proj_matrix = 1;
        for(j=0;j<matrix_order*matrix_order*2;j++){
            kv[i].dftbmat[j] = 0.0;
        }
    }

    int ind=0;
    int nOrbAtom=1;
    int *iSquare=malloc(nAllAtoms*sizeof(int));
    k=0;
    for(i=0; i<nel; i++){
        if(lmax[i]==0) nOrbAtom=1;
        if(lmax[i]==1) nOrbAtom=1+3;
        if(lmax[i]==2) nOrbAtom=1+3+5;
        nOrbAtom *=2;  // careful
        for(j=0; j<nat[i]; j++){
            iSquare[k] = ind;
            ind += nOrbAtom;
            k++;
        }
    }

    while(wv < no_wave_vectors){
        int num_each_subblock=0;
        fread(&(kv[wv].k_vector[0]), sizeof(double), 1, myPtr);
        fread(&(kv[wv].k_vector[1]), sizeof(double), 1, myPtr);
        fread(&(kv[wv].k_vector[2]), sizeof(double), 1, myPtr);

        //printf("k_vector is %9.5f %9.5f %9.5f\n", kv[wv].k_vector[0], kv[wv].k_vector[1], kv[wv].k_vector[2]);

        fread(&(kv[wv].num_proj_matrix), sizeof(int), 1, myPtr);
        fread(&num_each_subblock, sizeof(int), 1, myPtr);

        //printf("Project matrixs number is %d\n", kv[wv].num_proj_matrix);
        //printf("each project matrix consists of %d subblocks\n", num_each_subblock);

        if(kv[wv].num_proj_matrix == 1){
            printf("dmat must be E matix, so just return an matrix_order * matrix_order\n");
            for(j=0;j<matrix_order;j++){
                kv[i].dftbmat[j*matrix_order*2 + j*2] = 1.0;
            }
            continue;
        }


        kv[wv].proj_matrix = malloc(sizeof(ProjectMatrix)*kv[wv].num_proj_matrix);

        for(i=0; i<kv[wv].num_proj_matrix; i++){
            for(j=0;j<nel;j++){ // diff with mat index, 0 is start
                kv[wv].proj_matrix[i].num_subblock = num_each_subblock;
                kv[wv].proj_matrix[i].subblock = malloc(sizeof(Block)*(num_each_subblock));

            }
        }

        int index_dftb=0;

        for(i=0; i<kv[wv].num_proj_matrix; i++){
            int J=0;
            int JD=0;
            int num_subblock=0;
            int num_same_atoms=0;
            int ichem=0;
            int L=0;
            int itotal=0;
            int ndi=0;
            int cols=0;
            int rows=0;

            fread(&J, sizeof(int), 1, myPtr);
            fread(&JD, sizeof(int), 1, myPtr);
            
            
            kv[wv].proj_matrix[i].J = J;
            kv[wv].proj_matrix[i].JD = JD;

            num_subblock = kv[wv].proj_matrix[i].num_subblock;
            num_same_atoms = nat[ichem];
#if DEBUG
            printf("T-matrix for J = %d, JD = %d\n", J, JD);
            printf("num subblock is %d\n", num_subblock);
            printf("num same atoms in unit cell is %d\n", num_same_atoms);
#endif
            
            kv[wv].proj_matrix[i].subblock=malloc(sizeof(Block)*num_subblock);


            for(j=0; j<num_subblock; j++){
                // ichem L, num_same_atoms, cols, rows
                //printf("subblock index %d\n", j);
                fread(&ichem, sizeof(int), 1, myPtr);
                fread(&L, sizeof(int), 1, myPtr);
                fread(&cols, sizeof(int), 1, myPtr); // itotal
                fread(&rows, sizeof(int), 1, myPtr); // ndi

                num_same_atoms = nat[ichem-1];
                kv[wv].proj_matrix[i].subblock[j].num_same_atoms = num_same_atoms;
                // kv[wv].proj_matrix[i].subblock[j].subindex = j;
                kv[wv].proj_matrix[i].subblock[j].ichem = ichem;
                kv[wv].proj_matrix[i].subblock[j].lmax = lmax[ichem-1];
#if DEBUG
                printf("L is %d and lmax[ichem] is %d\n", L, lmax[ichem-1]);
                printf("ichem is %d\n", ichem);
                printf(" cols is %d rows is %d\n", cols, rows);
                printf("nat[ichem] is %d \n", nat[ichem-1]);
                printf(" cols*rows*2 is %d\n", cols*rows*2);
#endif
                assert(L <= lmax[ichem-1]);

                // point bb just for convenient
                bb=&(kv[wv].proj_matrix[i].subblock[j]);
                bb->b =malloc(sizeof(double)*2*cols*rows);

                bb->L = L;
                bb->cols = cols;
                bb->rows = rows;
                //bb->subindex = j;

                for(k=0; k<cols*rows*2; k++){
                    fread(&(bb->b[k]), sizeof(double), 1, myPtr);
                }
                // change to dftbmat

                // printf subblock
                /*
                for(k=0; k<cols; k++){
                    printf("Column %d \n", k);
                    for(int kj=0; kj<rows; kj++){
                        printf("(%9.5f %9.5f) ", bb->b[k*rows*2+kj*2], bb->b[k*rows*2+kj*2+1]);
                    }
                    printf("\n");
                }
                */
            }
        }

        int index=0;
        int l3[3]={2,8,18};
        int col_index=0;
        int row_index=0;
        int shift_single=0;

        int num_subblock=0;
        int subblock_index=0;
        int num_s_atoms=0;
        int ichem =0;
        int temp_L=0;
        int ichem_index =0;
        int first_row_step =0;
        int current_block_cols=0;
        int current_block_rows=0;


        int l0_col_index=0;
        int l1_col_index=0;
        int l2_col_index=0;
        int tmp_col_row =0;
        int b_index=0;
        int spd[3];
        while(ichem_index < nel+1){
            subblock_index=0;
            
            for(i=0; i<kv[wv].num_proj_matrix; i++)
                //for(i=kv[wv].num_proj_matrix-1; i>=0; i--)
            {
                
                for(j=0; j<kv[wv].proj_matrix[i].num_subblock; j++){
                    ichem = kv[wv].proj_matrix[i].subblock[j].ichem;
                    bb=&(kv[wv].proj_matrix[i].subblock[j]);
                    if(lmax[ichem-1] == 0) shift_single = l3[0];
                    else if(lmax[ichem-1] == 1) shift_single = l3[1];
                    else if(lmax[ichem-1] == 2) shift_single = l3[2];

                    b_index=0;
                    if(ichem == ichem_index){ // ichem == ichem_index
                        num_s_atoms =bb->num_same_atoms;
                        current_block_cols = bb->cols;
                        current_block_rows = bb->rows;
                        //printf("current_block_cols is %d %d ichem is %d shift_single is %d\n", current_block_cols, current_block_rows, ichem, shift_single);

                        temp_L = bb->L;

                        if(temp_L == 0){
                            first_row_step = shift_block[ichem-1]*2;
#if DEBUG
                            printf("first_row_step = %d\n", first_row_step);
#endif
                            for(int c=0; c < current_block_cols; c++){
                                for(k=0; k<current_block_rows; k++){
                                    row_index = first_row_step + shift_single*k;
                                    b_index = c*current_block_rows*2 + k*2;

#if DEBUG
                                    printf(" b_index %d\n", b_index);
                                    printf(" col and row %d %d\n", l0_col_index, row_index);
#endif
                                    tmp_col_row = (l0_col_index*matrix_order*2+row_index);
                                    kv[wv].dftbmat[tmp_col_row]=bb->b[b_index];
                                    kv[wv].dftbmat[tmp_col_row+1]=bb->b[b_index+1];
                                    //printf("%d (%9.5f %9.5f)", tmp_col_row, kv[wv].dftbmat[tmp_col_row], kv[wv].dftbmat[tmp_col_row+1]);
                                }
                                l0_col_index = l0_col_index + shift_single/2;
                            }
                            //l0_col_index++;
                        }

                        else if(temp_L == 1){

                            first_row_step = shift_block[ichem-1]*2 + 2;
                            
                            l1_col_index = l0_col_index - current_block_cols*shift_single/6 + 1 ;
#if DEBUG
                            printf("first_row_step = %d\n", first_row_step);
                            printf("current_block_cols/3  and number_s_atoms is %d  and %d\n", current_block_cols/3, num_s_atoms);
#endif 
                            for(int c=0; c<current_block_cols/3; c++){
                                row_index = first_row_step;
                              
                                for(k=0; k<current_block_rows/3; k++){
                                    
                                    // b_index  k*2*3;
#if DEBUG
                                    printf(" row_index is %d and b_index is %d\n", row_index, b_index);
                                    printf(" col and row %d %d\n", l1_col_index, row_index);
                                    printf(" col and row %d %d\n", l1_col_index, row_index+2);
                                    printf(" col and row %d %d\n", l1_col_index, row_index+4);
#endif
                                    tmp_col_row = l1_col_index * matrix_order*2 + row_index;
                                    // px
                                    kv[wv].dftbmat[tmp_col_row]=bb->b[b_index];
                                    kv[wv].dftbmat[tmp_col_row+1]=bb->b[b_index+1];
                                    // py
                                    kv[wv].dftbmat[tmp_col_row+2]=bb->b[b_index+2];
                                    kv[wv].dftbmat[tmp_col_row+3]=bb->b[b_index+3];
                                    // pz
                                    kv[wv].dftbmat[tmp_col_row+4]=bb->b[b_index+4];
                                    kv[wv].dftbmat[tmp_col_row+5]=bb->b[b_index+5];

                                    b_index += 6;
                                    row_index += shift_single;
                                }
                                row_index = first_row_step;
                                for(k=0; k<current_block_rows/3; k++){
                                    tmp_col_row = (l1_col_index+1) * matrix_order*2 + row_index;
                                    //b_index += k*2*3;
#if DEBUG
                                    printf(" row_index is %d and b_index is %d\n", row_index, b_index);
                                    printf(" col and row %d %d\n", l1_col_index+1, row_index);
                                    printf(" col and row %d %d\n", l1_col_index+1, row_index+2);
                                    printf(" col and row %d %d\n", l1_col_index+1, row_index+4);
#endif
                                    kv[wv].dftbmat[tmp_col_row]=bb->b[b_index];
                                    kv[wv].dftbmat[tmp_col_row+1]=bb->b[b_index+1];
                                    kv[wv].dftbmat[tmp_col_row+2]=bb->b[b_index+2];
                                    kv[wv].dftbmat[tmp_col_row+3]=bb->b[b_index+3];
                                    kv[wv].dftbmat[tmp_col_row+4]=bb->b[b_index+4];
                                    kv[wv].dftbmat[tmp_col_row+5]=bb->b[b_index+5];

                                    b_index += 6;
                                    row_index += shift_single;
                                }
                                row_index = first_row_step;
                                for(k=0; k<current_block_rows/3; k++){
                                    tmp_col_row = (l1_col_index+2) * matrix_order*2 + row_index;
                                    // b_index +=k*2*3;
#if DEBUG
                                    printf(" row_index is %d and b_index is %d\n", row_index, b_index);
                                    printf(" col and row %d %d\n", l1_col_index+2, row_index);
                                    printf(" col and row %d %d\n", l1_col_index+2, row_index+2);
                                    printf(" col and row %d %d\n", l1_col_index+2, row_index+4);
#endif
                                    kv[wv].dftbmat[tmp_col_row]=bb->b[b_index];
                                    kv[wv].dftbmat[tmp_col_row+1]=bb->b[b_index+1];
                                    kv[wv].dftbmat[tmp_col_row+2]=bb->b[b_index+2];
                                    kv[wv].dftbmat[tmp_col_row+3]=bb->b[b_index+3];
                                    kv[wv].dftbmat[tmp_col_row+4]=bb->b[b_index+4];
                                    kv[wv].dftbmat[tmp_col_row+5]=bb->b[b_index+5];

                                    b_index += 6;
                                    row_index += shift_single;
                                }
                            
                                 l1_col_index += shift_single/2;
                            }
                        }
                        else if(temp_L == 2){
                            first_row_step = shift_block[ichem-1]*2 + 8;
                            printf("first_row_step = %d\n", first_row_step);
                            if(current_block_cols*shift_single%10 != 0 || current_block_cols==0){
                                printf("current_block is empty\n");
                            }

                            l2_col_index = l1_col_index - current_block_cols*shift_single/10+ 3 + 1 ;
                            for(int c=0; c < current_block_cols/5; c++){
                                row_index = first_row_step;
                                for(k=0; k<num_s_atoms; k++){
#if DEBUG
                                    printf(" col and row %d %d\n", l2_col_index, row_index);
                                    printf(" col and row %d %d\n", l2_col_index, row_index+1);
                                    printf(" col and row %d %d\n", l2_col_index, row_index+2);
                                    printf(" col and row %d %d\n", l2_col_index, row_index+3);
                                    printf(" col and row %d %d\n", l2_col_index, row_index+4);
#endif
                                    tmp_col_row = 2*l2_col_index* matrix_order + row_index;
                                   
                                    kv[wv].dftbmat[tmp_col_row]=bb->b[b_index];
                                    kv[wv].dftbmat[tmp_col_row+1]=bb->b[b_index+1];
                                    kv[wv].dftbmat[tmp_col_row+2]=bb->b[b_index+2];
                                    kv[wv].dftbmat[tmp_col_row+3]=bb->b[b_index+3];
                                    kv[wv].dftbmat[tmp_col_row+4]=bb->b[b_index+4];
                                    kv[wv].dftbmat[tmp_col_row+5]=bb->b[b_index+5];
                                    kv[wv].dftbmat[tmp_col_row+6]=bb->b[b_index+6];
                                    kv[wv].dftbmat[tmp_col_row+7]=bb->b[b_index+7];
                                    kv[wv].dftbmat[tmp_col_row+8]=bb->b[b_index+8];
                                    kv[wv].dftbmat[tmp_col_row+9]=bb->b[b_index+9];

                                    b_index +=10;
                                    row_index += shift_single;
                                }
                                row_index = first_row_step;
                                for(k=0; k<num_s_atoms; k++){
#if DEBUG
                                    printf(" col and row %d %d\n", l2_col_index+1, row_index);
                                    printf(" col and row %d %d\n", l2_col_index+1, row_index+1);
                                    printf(" col and row %d %d\n", l2_col_index+1, row_index+2);
                                    printf(" col and row %d %d\n", l2_col_index+1, row_index+3);
                                    printf(" col and row %d %d\n", l2_col_index+1, row_index+4);
#endif
                                    tmp_col_row = 2*(l2_col_index+1) * matrix_order + row_index;
                                   
                                    kv[wv].dftbmat[tmp_col_row]=bb->b[b_index];
                                    kv[wv].dftbmat[tmp_col_row+1]=bb->b[b_index+1];
                                    kv[wv].dftbmat[tmp_col_row+2]=bb->b[b_index+2];
                                    kv[wv].dftbmat[tmp_col_row+3]=bb->b[b_index+3];
                                    kv[wv].dftbmat[tmp_col_row+4]=bb->b[b_index+4];
                                    kv[wv].dftbmat[tmp_col_row+5]=bb->b[b_index+5];
                                    kv[wv].dftbmat[tmp_col_row+6]=bb->b[b_index+6];
                                    kv[wv].dftbmat[tmp_col_row+7]=bb->b[b_index+7];
                                    kv[wv].dftbmat[tmp_col_row+8]=bb->b[b_index+8];
                                    kv[wv].dftbmat[tmp_col_row+9]=bb->b[b_index+9];

                                    b_index +=10;
                                    row_index += shift_single;
                                }
                                row_index = first_row_step;
                                for(k=0; k<num_s_atoms; k++){
#if DEBUG
                                    printf(" col and row %d %d\n", l2_col_index+2, row_index);
                                    printf(" col and row %d %d\n", l2_col_index+2, row_index+1);
                                    printf(" col and row %d %d\n", l2_col_index+2, row_index+2);
                                    printf(" col and row %d %d\n", l2_col_index+2, row_index+3);
                                    printf(" col and row %d %d\n", l2_col_index+2, row_index+4);
#endif
                                    tmp_col_row = 2*(l2_col_index+2) * matrix_order + row_index;
                                  
                                    kv[wv].dftbmat[tmp_col_row]=bb->b[b_index];
                                    kv[wv].dftbmat[tmp_col_row+1]=bb->b[b_index+1];
                                    kv[wv].dftbmat[tmp_col_row+2]=bb->b[b_index+2];
                                    kv[wv].dftbmat[tmp_col_row+3]=bb->b[b_index+3];
                                    kv[wv].dftbmat[tmp_col_row+4]=bb->b[b_index+4];
                                    kv[wv].dftbmat[tmp_col_row+5]=bb->b[b_index+5];
                                    kv[wv].dftbmat[tmp_col_row+6]=bb->b[b_index+6];
                                    kv[wv].dftbmat[tmp_col_row+7]=bb->b[b_index+7];
                                    kv[wv].dftbmat[tmp_col_row+8]=bb->b[b_index+8];
                                    kv[wv].dftbmat[tmp_col_row+9]=bb->b[b_index+9];

                                    b_index +=10;
                                    row_index += shift_single;
                                }
                                row_index = first_row_step;
                                for(k=0; k<num_s_atoms; k++){
#if DEBUG
                                    printf(" col and row %d %d\n", l2_col_index+3, row_index);
                                    printf(" col and row %d %d\n", l2_col_index+3, row_index+1);
                                    printf(" col and row %d %d\n", l2_col_index+3, row_index+2);
                                    printf(" col and row %d %d\n", l2_col_index+3, row_index+3);
                                    printf(" col and row %d %d\n", l2_col_index+3, row_index+4);
#endif
                                   tmp_col_row = 2*(l2_col_index+3) * matrix_order + row_index;
                                
                                    kv[wv].dftbmat[tmp_col_row]=bb->b[b_index];
                                    kv[wv].dftbmat[tmp_col_row+1]=bb->b[b_index+1];
                                    kv[wv].dftbmat[tmp_col_row+2]=bb->b[b_index+2];
                                    kv[wv].dftbmat[tmp_col_row+3]=bb->b[b_index+3];
                                    kv[wv].dftbmat[tmp_col_row+4]=bb->b[b_index+4];
                                    kv[wv].dftbmat[tmp_col_row+5]=bb->b[b_index+5];
                                    kv[wv].dftbmat[tmp_col_row+6]=bb->b[b_index+6];
                                    kv[wv].dftbmat[tmp_col_row+7]=bb->b[b_index+7];
                                    kv[wv].dftbmat[tmp_col_row+8]=bb->b[b_index+8];
                                    kv[wv].dftbmat[tmp_col_row+9]=bb->b[b_index+9];

                                    b_index =10;
                                    row_index += shift_single;
                                }
                                row_index = first_row_step;
                                for(k=0; k<num_s_atoms; k++){
                                    tmp_col_row = 2*(l2_col_index+4) * matrix_order + row_index;
                                 
                                    kv[wv].dftbmat[tmp_col_row]=bb->b[b_index];
                                    kv[wv].dftbmat[tmp_col_row+1]=bb->b[b_index+1];
                                    kv[wv].dftbmat[tmp_col_row+2]=bb->b[b_index+2];
                                    kv[wv].dftbmat[tmp_col_row+3]=bb->b[b_index+3];
                                    kv[wv].dftbmat[tmp_col_row+4]=bb->b[b_index+4];
                                    kv[wv].dftbmat[tmp_col_row+5]=bb->b[b_index+5];
                                    kv[wv].dftbmat[tmp_col_row+6]=bb->b[b_index+6];
                                    kv[wv].dftbmat[tmp_col_row+7]=bb->b[b_index+7];
                                    kv[wv].dftbmat[tmp_col_row+8]=bb->b[b_index+8];
                                    kv[wv].dftbmat[tmp_col_row+9]=bb->b[b_index+9];
#if DEBUG
                                    printf(" col and row %d %d\n", l2_col_index+4, row_index);
                                    printf(" col and row %d %d\n", l2_col_index+4, row_index+1);
                                    printf(" col and row %d %d\n", l2_col_index+4, row_index+2);
                                    printf(" col and row %d %d\n", l2_col_index+4, row_index+3);
                                    printf(" col and row %d %d\n", l2_col_index+4, row_index+4);
                            
#endif
                                    b_index=10;
                                    row_index += shift_single;
                                }
                                l2_col_index += shift_single/2;
                            }

                        }
                        else {
                            printf("L > 2 is not supported yet \n");
                            exit(1);
                        }
                        //subblock_index++;
                    }
                }
            }
            ichem_index++;
        }
        //printf("output  D mat for kpoint (%9.5f %9.5f %9.5f)\n", kv[wv].k_vector[0], kv[wv].k_vector[1], kv[wv].k_vector[2]);

        //printDmat(kv[wv].dftbmat, matrix_order);
        double_array_to_complex(kv[wv].dftbmat, dmatrixs, matrix_order);

        //memcpy(dmatrixs[wv], kv[wv].dftbmat, matrix_order*matrix_order*2*sizeof(double));
        wv++;
    }


    // free malloc
    free(nat);
    free(shift_block);
    free(lmax);

    for(i=0; i<no_wave_vectors; i++){
        for(j=0; j<kv[i].num_proj_matrix; j++){
            for(k=0; k<kv[i].proj_matrix[j].num_subblock; k++){
                free(kv[i].proj_matrix[j].subblock[k].b);
            }
            free(kv[i].proj_matrix[j].subblock);
        }
        free(kv[i].proj_matrix);
        free(kv[i].nat);
        free(kv[i].lmax);
        free(kv[i].dftbmat);
        //free(kv[i].dmat);
    }

    free(kv);
    free(nAllOrb);
    free(iSquare);
}
/*
int main(int argc, char *argv[]){
    int no_wave_vectors=1;
    int matrix_order = 96;

    double complex **dmatrixs=NULL;

    if(dmatrixs == NULL){
        printf("WARNING: dmatrixs is NULL\n ");
        dmatrixs = (double complex **)malloc(sizeof(double complex*)*no_wave_vectors);
        for(int i=0;i<no_wave_vectors;i++)
        {
            dmatrixs[i] = (double complex*)malloc(sizeof(double complex)*matrix_order*matrix_order);
        }
    }

    parse_symproj_(argv[1], 1, 96, dmatrixs);
    //printDmat(dmatrixs[0], 96);
    //double_array_to_complex(dmatrixs[0], dmat, matrix_order);
    printDCmat(dmatrixs[0], matrix_order);
    //free(dmatrixs[0]);
    free(dmatrixs);
    //free(dmat);
}
*/
