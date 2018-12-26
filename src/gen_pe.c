#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

double distance(double a[3], double b[3])
{
    double x,y,z;
    x = a[0]-b[0];
    y = a[1]-b[1];
    z = a[2]-b[2];

    return sqrt(x*x + y*y + z*z);
}

int run_task(const char *task) {
    int aux_system;
    aux_system = system(task);
    return 1;
}


void simple_write_to_psi4(const char *fname, const char *element, const char *basis, const char *scf_type, double b[3])
{
    char *fmt;
    FILE *o;
    o = fopen(fname,"w");
    fprintf(o, "#scan C-Ha\n");
    fprintf(o, "\n");
    fprintf(o, "memory 24000 mb\n");
    fprintf(o, "\n");
    fprintf(o, "molecule abc_ar {\n");
    fprintf(o, "0 1\n");
    fprintf(o, "C                  0.473555413895     0.820128990805    -0.000380416252\n");
    fprintf(o, "C                  1.863816689780     0.816506827686     0.000070474802\n");
    fprintf(o, "C                  2.555558246250     2.025940988460    -0.000487079124\n");
    fprintf(o, "C                  1.865757314295     3.231247121135    -0.001474401526\n");
    fprintf(o, "C                  0.477033715570     3.226102333495    -0.001912343544\n");
    fprintf(o, "C                 -0.224602666951     2.022377324919    -0.001372423211\n");
    fprintf(o, "H                  2.407642774474    -0.121187157469     0.000842384598\n");
    fprintf(o, "H                  3.639667771283     2.009068574457    -0.000134980633\n");
    fprintf(o, "H                  2.407729021001     4.169868300643    -0.001901027610\n");
    fprintf(o, "H                 -0.079555907814     4.156578240828    -0.002682058058\n");
    fprintf(o, "H                 -1.308582549918     2.024600882960    -0.001722146049\n");
    fprintf(o, "%s                %9.5f    %9.5f    %9.5f\n", element, b[0], b[1], b[2]);
    fprintf(o, "}\n");
    fprintf(o, "set g_convergence gau_loose\n");
    fprintf(o, "set basis %s\n", basis);
    fprintf(o, "set guess sad\n");
    fprintf(o, "set freeze_core True\n");
    fprintf(o, "set scf_type %s\n", scf_type);
    fprintf(o, "e=energy('b3lyp')\n");
    fprintf(o, "print(e)\n");


    fclose(o);
}
int main(){
    double    pos_C[3]={ -0.224602666951, 2.022377324919,-0.001372423211};
    double pos_Cl[3]={ -0.719089110597,    -1.245364631768,     0.000569845141};
    double len;
    len = distance(pos_C, pos_Cl);
    printf("%9.5f\n", len);
    double change_range=1.0;
    int num_p=20;
    double change_list[20];
    double step= change_range/num_p;
    for(int i=0; i<num_p; i++)
    {
        change_list[i]= -1*(change_range/2) + i*step;
    }


    char element[4] ="I";
    char basis[20] ="6-31G*";
    char scf_type[20] ="df"; // direct
    char basePsi4[60]="scan1_";
    for(int i=0; i<num_p; i++)
    {
        char psifile[60];
        char psiout[60];
        double b[3];
        double dis;
        char intstr[10];
        sprintf(intstr, "%d", i);
        b[0] = pos_Cl[0];
        b[2] = pos_Cl[2];
        b[1] = pos_Cl[1] + change_list[i];
        dis = distance(pos_C, b);
        printf("distance is %9.5f\n", dis);
        strcpy(psifile, basePsi4);
        strcat(psifile, intstr);
        strcpy(psiout, psifile);
        strcat(psifile, ".in");
        strcat(psiout, ".out");
        simple_write_to_psi4(psifile, element, basis, scf_type, b);
        //run_task("psi4 -i %s -o %s -n 4", psifile, psiout);
        char command[200];
        sprintf(command, "grep \"Total Energy =                     \" %s",psiout);
        //run_task(command);
        //printf("psi4 -i %s -n 4\n", psifile);
    }

}


