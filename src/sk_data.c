#include"sk_data.h"


void print_TOldSKData(TOldSKData *a){
    printf("SKData->dist is %9.5f\n", a->dist);
    printf("SKData->nGrid is %d\n", a->nGrid);
    printf("SKData->skSelf is %9.5f %9.5f %9.5f %9.5f\n", a->skSelf[0], a->skSelf[1], a->skSelf[2], a->skSelf[3]);
    printf("SKData->skHubbu is %9.5f %9.5f %9.5f %9.5f\n", a->skHubbU[0], a->skHubbU[1], a->skHubbU[2], a->skHubbU[3]);
    printf("SKData->skOccp is %9.9f %9.9f %9.5f %9.5f\n", a->skOccp[0], a->skOccp[1], a->skOccp[2], a->skOccp[3]);
    printf("SKData->mass is %9.5f\n", a->mass);
    for(int i=0;i<a->nGrid;i++)
    {
        int ii=i*nSKInterOld;
        printf("skHam %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f\n",
               a->skHam[ii], a->skHam[ii+1], a->skHam[ii+2], a->skHam[ii+3], a->skHam[ii+4],
               a->skHam[ii+5], a->skHam[ii+6], a->skHam[ii+7], a->skHam[ii+8], a->skHam[ii+9]);
        printf("skOver %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f %9.9f\n",
               a->skOver[ii], a->skOver[ii+1], a->skOver[ii+2], a->skOver[ii+3], a->skOver[ii+4],
               a->skOver[ii+5], a->skOver[ii+6], a->skOver[ii+7], a->skOver[ii+8], a->skOver[ii+9]);
    }
}

void print_TRepSplineIn(TRepSplineIn *a){
    printf("RepSplinein->nInt is %d\n", a->nInt);
    printf("RepSplineIn->cutoff is %9.5f\n", a->cutoff);
    printf("Repsplinein->expcoeffs is %9.5f %9.5f %9.5f\n", a->expCoeffs[0], a->expCoeffs[1], a->expCoeffs[2]);
    printf("Repsplinein->spLastcoeffs is %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
           a->spLastCoeffs[0], a->spLastCoeffs[1], a->spLastCoeffs[2],
           a->spLastCoeffs[3], a->spLastCoeffs[4], a->spLastCoeffs[5]);
    for(int i=0;i<a->nInt;i++)
    {
        int ii=i*4;
        printf("RepSplinein->xStart is %9.5f\n", a->xStart[i]);
        printf("RepSplinein->spCoeff is %9.5f %9.5f %9.5f %9.5f\n",
               a->spCoeff[ii], a->spCoeff[ii+1], a->spCoeff[ii+2], a->spCoeff[ii+3]);
    }
}

void print_TRepPolyIn(TRepPolyIn *a){
    printf("RepPolyIn->cutoff is %9.5f\n", a->cutoff);
    printf("RepPolyIn->polyCoeffs is\n");
    printf("%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
           a->polyCoeffs[0], a->polyCoeffs[1], a->polyCoeffs[2], a->polyCoeffs[3],
           a->polyCoeffs[4], a->polyCoeffs[5], a->polyCoeffs[6], a->polyCoeffs[7]);
}

void free_TOldSKData(TOldSKData *a){
    free(a->skHam);
    free(a->skOver);
}

void free_TRepSplineIn(TRepSplineIn *a){
    free(a->xStart);
    free(a->spCoeff);
}


int parse_sk_data_line(char fbuffer[], double *result){
    char *p=NULL, *q=NULL, buf[24];
    p=blank_advance(fbuffer);
    int l=0, i;
    int has_star=0;
    int index=0;
    double t_value;
    while(1){
        if ( (*p == EOS) || (*p == '!') ) break;
        q = nonblank_advance(p);
        sscanf(p, "%s", buf);
        has_star=0;
        l = strlen(buf);
        for(i=0; i<l; i++)
        {
            if(buf[i] == '*')
            {
                has_star=1;
                char *fir_res = strtok(buf, "*");
                char *last_res = strtok(NULL, "*");
                int number = atoi(fir_res);
                t_value = atof(last_res);
                for(int k=0;k<number;k++){
                    result[index++]=t_value;
                }
            }
        }

        if(has_star==0 && *p!=EOS){
            t_value = atof(buf);
            result[index++]=t_value;
        }
        p=blank_advance(q);
    }
    return index;
}


int OldSKData_readFromFile(TOldSKData *skData,
                           TRepSplineIn *repSplineIn,
                           TRepPolyIn *repPolyIn,
                           const char *filename,
                           const int homo,
                           const int iSp1,
                           const int iSp2)
{
    int nInt;
    double rTmp, rTmp2, rDummy;
    int k;
    //int iostat;
    int tExtended;
    int nShell;
    double dist;

    FILE *fid;
    fid = fopen(filename, "r");
    if(!fid) exit(0);

    int rs;
    char fbuffer[FGETS_MAXSIZE];

    double coeffs[8]; // 2-9 coeffs
    double polyCutoff;
    double data_line[64];
    if(NULL == fgets(fbuffer, 1024, fid)) return 0;

    int nGrid=0;
    if(fbuffer[0]=='@'){
        tExtended=1;
        nShell=4;
        fprintf(stderr,"not finished all yet \n");
        exit(1);
    }else{
        tExtended=0;
        nShell=3;
    }

    double line0[3];
    rs=parse_sk_data_line(fbuffer, line0);
    dist = line0[0];
    nGrid= (int)line0[1];
    skData->nGrid = nGrid;

    if(homo){
        k=atoi(fgets(fbuffer, 1024, fid));
        // d p s rDummy d p s d p s
        double line1[10];
        rs=parse_sk_data_line(fbuffer, line1);
        skData->skSelf[0]=line1[0];
        skData->skSelf[1]=line1[1];
        skData->skSelf[2]=line1[2];
        rDummy = line1[3];

        skData->skHubbU[0]=line1[4];
        skData->skHubbU[1]=line1[5];
        skData->skHubbU[2]=line1[6];

        skData->skOccp[0]=line1[7];
        skData->skOccp[1]=line1[8];
        skData->skOccp[2]=line1[9];
        //rs = sscanf(fbuffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &skData->skSelf[0], &skData->skSelf[1],&skData->skSelf[2],
        //            &rDummy, &skData->skHubbU[0], &skData->skHubbU[1], &skData->skHubbU[2],
        //            &skData->skOccp[0], &skData->skOccp[1], &skData->skOccp[2]);
/*
        printf("skData-> %9.5f %9.5f %9.5f| %9.5f %9.5f %9.5f| %9.5f %9.5f %9.5f \n",
               skData->skSelf[0], skData->skSelf[1], skData->skSelf[2],
               skData->skHubbU[0], skData->skHubbU[1], skData->skHubbU[2],
               skData->skOccp[0], skData->skOccp[1], skData->skOccp[2]);
*/
        if(rs<10){
            fprintf(stderr, "read sk file line 2 error, number is %d\n", rs);
            exit(0);
        }
        //FGets(fbuffer, fid);
        k = fgets(fbuffer, 1024, fid);
        // mass c2-c9 rcut d1-d10;
        double line3[20];
        int index=parse_sk_data_line(fbuffer, line3);
        if(!index){
            error_exit("parse_sk_data_line error\n");
        }
        skData->mass = line3[0]* amu__au; // convert to atomic units
        for(int i=0; i<8;i++){
            coeffs[i]=line3[i+1];
        }
        polyCutoff= line3[9];
    }else{
        k=fgets(fbuffer, 1024, fid);
        parse_sk_data_line(fbuffer, data_line);
        for(int i=0; i<8;i++){
            coeffs[i]=data_line[i+1];
        }
        polyCutoff=data_line[9];
    }

    // present(repPloyIn) save data to repPloyin
    for(int i=0;i<8;i++)
        repPolyIn->polyCoeffs[i] = coeffs[i];
    repPolyIn->cutoff=polyCutoff;

    // init skdata->skHam skOver
    int size_skHam;
    if(tExtended){
        size_skHam = sizeof(double)*skData->nGrid*nSKInter;
    }else{
        size_skHam = sizeof(double)*skData->nGrid*nSKInterOld;
    }

    skData->skHam =(double *)malloc(size_skHam);
    skData->skOver=(double *)malloc(size_skHam);


    for(int i=0; i<skData->nGrid; i++){
        k = fgets(fbuffer, 1024, fid);
        int num_ele=parse_sk_data_line(fbuffer, data_line);
        //printf(" %d num_ele is %d\n",i, num_ele);
        if(num_ele==20){
            for(int j=0; j<nSKInterOld; j++){
                skData->skHam[i*nSKInterOld+j]=data_line[j];
                skData->skOver[i*nSKInterOld+j]=data_line[nSKInterOld+j];
            }
        }else if(num_ele==40){
           for(int j=0; j<nSKInter; j++){
                skData->skHam[i*nSKInter+j]=data_line[j];
                skData->skOver[i*nSKInter+j]=data_line[nSKInter+j];
            }
        }else {
            fprintf(stderr, "%d in nGrid is %d\n", i, nGrid);
            fprintf(stderr, "error read ham and over data from %s\n", filename);
            return 0;
        }
    }
    // handle spline data

    char line[64];
    while(1){
        k=fgets(fbuffer,1024,fid);
        rs = sscanf(fbuffer, "%s", line);
        if(!strcmp(line, "Spline")){
            break;
        }
    }
    k=fgets(fbuffer, 1024, fid);
    rs = sscanf(fbuffer, "%d %lf", &nInt, &repSplineIn->cutoff);
    k=fgets(fbuffer, 1024, fid);
    rs = sscanf(fbuffer, "%lf %lf %lf", &repSplineIn->expCoeffs[0], &repSplineIn->expCoeffs[1], &repSplineIn->expCoeffs[2]);

    repSplineIn->xStart = (double *)malloc(sizeof(double)*nInt);
    repSplineIn->spCoeff = (double *)malloc(sizeof(double)*nInt*4);
    repSplineIn->nInt = nInt;
    rTmp = 0.0;

    for(int i=0; i<nInt; i++){
        k=fgets(fbuffer, 1024, fid);
        rs = sscanf(fbuffer, "%lf %lf %lf %lf %lf %lf", &repSplineIn->xStart[i], &rTmp2, &repSplineIn->spCoeff[i*4],
                    &repSplineIn->spCoeff[i*4+1], &repSplineIn->spCoeff[i*4+2], &repSplineIn->spCoeff[i*4+3]);
        if( i !=0 && fabs(repSplineIn->xStart[i] - rTmp) > 0.00001){
            fprintf(stderr, "Repulsive not continuous for specie pair number %d and %d", iSp1, iSp2);
            return 0;
        }
        rTmp = rTmp2;
    }
    repSplineIn->cutoff = rTmp2;
    fclose(fid);
    return 1;
};

void Alloc_TSKData(TSKData *sk){
    sk->skData= (TOldSKData *)malloc(sizeof(TOldSKData));
    sk->repSplineIn = (TRepSplineIn *)malloc(sizeof(TRepSplineIn));
    sk->repPolyIn = (TRepPolyIn *) malloc(sizeof(TRepPolyIn));

}

void Free_TSKData(TSKData *sk){
    free(sk->skData);
    free(sk->repSplineIn);
    free(sk->repPolyIn);
}

void init_OSlakoCont(OSlakoCont *a, int nSpecie){
    if(!a->tInit)
    {
        a->nSpecie = nSpecie;
        a->mInt = 0;
        a->tDataOK=0;
        a->tInit=1;
        a->slakos=(PSlaKo_ *)malloc(nSpecie*nSpecie*sizeof(PSlaKo_));
        a->cutoff=0.0;
    }
    return;
}

void free_OSlakoCont(OSlakoCont *a){
    free(a->slakos);
}

int main()
{
    printf("read sk data\n");
    TOldSKData *skData=(TOldSKData *)malloc(sizeof(TOldSKData));
    TRepSplineIn *repSplineIn=(TRepSplineIn *)malloc(sizeof(TRepSplineIn));
    TRepPolyIn *repPolyIn=(TRepPolyIn *)malloc(sizeof(TRepPolyIn));
    int homo=1;
    int iSp1=6;
    int iSp2=6;
    printf("in old\n");
    OldSKData_readFromFile(skData, repSplineIn, repPolyIn, "B-B.skf", homo, iSp1, iSp2);
    printf("end old\n");
    print_TOldSKData(skData);
    print_TRepPolyIn(repPolyIn);
    print_TRepSplineIn(repSplineIn);

    free_TOldSKData(skData);
    free_TRepSplineIn(repSplineIn);
    free(skData);
    free(repSplineIn);
    free(repPolyIn);

}

