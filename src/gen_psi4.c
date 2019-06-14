#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>
#include<getopt.h>
#include<stdint.h>
#include<regex.h>

#define APP_NAME "gen_psi4"
#define VERSION "0.3"

bool opt_debug = false;
static char psi4_filename[128];
static char xyz_filename[128];
static char output_psi4[128];
static int num_scale= 4;
static int density = 1;
static int charge = 0;
static int spin = 1;
static int mem = 24000;
static int inputfile_type =0;
static char psi4_job_type[64];
static int num_threads =4;
static char const usage[]="\
\n\
\n\
\t *** gen_psi4 for generating psi4 in file *** \n\
\n\
\n\
";

static char const short_options[]="hi:t:n:";

static struct option const options[] = {
    {"file", 1, NULL, 'i'},
    {"psi4_job_type", 1, NULL, 't'},
    {"num_threads", 0, NULL, 'n'},
    {"help", 0, NULL, 'h'}
};

static void show_usage_and_exit(int status){
    if(status) fprintf(stderr, "Try" APP_NAME " --help' for more information\n");
    else
        printf(usage);
    exit(status);
}

static void show_version_and_exit(void){
    printf("%s\n", VERSION);
    exit(0);
}

// support xyz or psi4's out
int get_file_type(const char *filename){
    const char ch = '.';
    char *ret;
    ret = strrchr(filename, ch);
    if(strcmp(ret,".xyz") == 0){
        return 1;
    }
    else if(strcmp(ret, ".out") == 0) {
        return 2;
    }
    return 0;
}

static void parse_arg(int key, char *arg)
{
    int enum_type=0;
    switch(key){
        case 'i':
            enum_type = get_file_type(arg);
            if(enum_type == 2){
                strcpy(psi4_filename, arg);
                inputfile_type = enum_type;
            }
            else if(enum_type == 1){
                strcpy(xyz_filename, arg);
                inputfile_type = enum_type;
            }
            break;
        case 't':
            strcpy(psi4_job_type, arg);
            break;
        case 'n':
            num_threads=atoi(arg);
            break;
        case 'h':
            show_usage_and_exit(0);
        default:
            show_usage_and_exit(0);
    }
}


static void parse_cli(int argc, char *argv[])
{
    int key;
    while(1)
    {
      key = getopt_long(argc, argv, short_options, options, NULL);
      if(key<0)
        break;
      parse_arg(key, optarg);
    }
  if(optind <argc){
    fprintf(stderr, "%s: unsupported non-option argument '%s'\n", argv[0], argv[optind]);
    show_usage_and_exit(0);
  }
}

#define FMF 0 // malloc/free debug  model
#define HPD 0
void *my_malloc(size_t size, const char *file, int line, const char *func);
void my_free(void *p, const char *file, int line, const char *func);

#define hp_malloc(X) my_malloc(X, __FILE__, __LINE__, __FUNCTION__)
#define hp_free(X) my_free(X, __FILE__, __LINE__, __FUNCTION__)


#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))


#define PI 3.141592653589793
#define EPSILON 8.8817841970012523e-016 /* 4.0 * DBL_EPSILON */
#define B2A  0.529177249  // 1Bohr = 0.529177249 Angstorm

#define TORAD(A)     ((A)*0.017453293)
#define DP_TOL     0.001

#define N_REFINE_STEPS 25

void *my_malloc(size_t size, const char *file, int line, const char *func) {
  void *p = malloc(size);
  if (!p) {
    printf("Error: couldn't allocate memory .Allocated = %s, %i, %s, %p[%li]\n",
           file, line, func, p, size);
    exit(EXIT_FAILURE);
  }
#if FMF
  printf("Allocated = %s, %i, %s, %p[%li]\n", file, line, func, p, size);
#endif
  return p;
}

void my_free(void *p, const char *file, int line, const char *func) {
#if FMF
  printf("Free = %s, %i, %s, %p[%li]\n", file, line, func, p, sizeof(p));
#endif
  free(p);
}


// map element type to int number
struct entry {
  char *str;
  int n;
};

static const double an2masses[119] = {
0.0,1.00782503207,4.00260325415,7.016004548,9.012182201,11.009305406,
12.0,14.00307400478,15.99491461956,18.998403224,19.99244017542,
22.98976928087,23.985041699,26.981538627,27.97692653246,30.973761629,
31.972070999,34.968852682,39.96238312251,38.963706679,39.962590983,
44.955911909,47.947946281,50.943959507,51.940507472,54.938045141,
55.934937475,58.933195048,57.935342907,62.929597474,63.929142222,
68.925573587,73.921177767,74.921596478,79.916521271,78.918337087,
85.910610729,84.911789737,87.905612124,88.905848295,89.904704416,
92.906378058,97.905408169,98.906254747,101.904349312,102.905504292,
105.903485715,106.90509682,113.90335854,114.903878484,119.902194676,
120.903815686,129.906224399,126.904472681,131.904153457,132.905451932,
137.905247237,138.906353267,139.905438706,140.907652769,141.907723297,
144.912749023,151.919732425,152.921230339,157.924103912,158.925346757,
163.929174751,164.93032207,165.930293061,168.93421325,173.938862089,
174.940771819,179.946549953,180.947995763,183.950931188,186.955753109,
191.96148069,192.96292643,194.964791134,196.966568662,201.970643011,
204.974427541,207.976652071,208.980398734,208.982430435,210.987496271,
222.017577738,222.01755173,228.031070292,227.027752127,232.038055325,
231.03588399,238.050788247,237.048173444,242.058742611,243.06138108,
247.07035354,247.07030708,251.079586788,252.082978512,257.095104724,
258.098431319,255.093241131,260.105504,263.112547,255.107398,259.114500,
262.122892,263.128558,265.136151,281.162061,272.153615,283.171792,283.176451,
285.183698,287.191186,292.199786,291.206564,293.214670};

struct entry dict[] = {
  "h", 1,
  "H", 1,
  "HE", 2,
  "He", 2,
  "Li", 3,
  "Be", 4,
  "B", 5,
  "c", 6,
  "C", 6,
  "n", 7,
  "N", 7,
  "o", 8,
  "O", 8,
  "f", 9,
  "F", 9,
  "Ne", 10,
  "Na", 11,
  "Mg", 12,
  "Al", 13,
  "SI", 14,
  "Si", 14,
  "p", 15,
  "P", 15,
  "s", 16,
  "S", 16,
  "CL", 17,
  "Cl", 17,
  "Ar", 18,
  "K", 19,
  "Ca", 20,
  "Sc", 21,
  "SC", 21,
  "TI", 22,
  "Ti", 22,
  "V", 23,
  "Cr", 24,
  "Mn", 25,
  "FE", 26,
  "Fe", 26,
  "Co", 27,
  "Ni", 28,
  "CU", 29,
  "Cu", 29,
  "ZN", 30,
  "Zn", 30,
  "GA", 31,
  "Ga", 31,
  "GE", 32,
  "Ge", 32,
  "As", 33,
  "Se", 34,
  "Br", 35,
  "Kr", 36,
};

int get_ele_number(char *key){
  int num=1;
  int i=0;
  char *name = dict[i].str;
  while(name) {
     if(strcmp(name, key) == 0)
       return dict[i].n;
     name = dict[++i].str;
  }
  return 0;
}

char * get_ele_name(const int num){
   char *name;
   int n=0;
   for(int i=0; i<53; i++)
   {
      n = dict[i].n;
      if (n==num)
        name = dict[i].str;
   }
   return name;
}

int get_file_length(const char *filename){
  FILE *fp;
  char buf[1000];
  int len =0;
  fp = fopen(filename, "r");
  if(!fp){
    fprintf(stderr, "unable to open %s\n", filename);
    return 0;
  }
  while(fgets(buf, 1000, fp)){
    len++;
  }
  fclose(fp);
  return len;
}

void read_xyz_file(const char *filepath, const int num_ele, int *ele_num, double *coord){
  int x;
  char c[40]; //sometimes the atom names are given which can be longer than 4 chars
  FILE *fid;
  fid= fopen(filepath, "r");
  int num_xyz_len=0;
  char s;
  char buf[100];
  fscanf(fid, "%d ", &num_xyz_len);
  fgets(buf, 100, fid);
  //fscanf(fid, "%s ", &s);
  // have some error, need to fix
  for(int i=0; i< num_ele; i++){
    fscanf(fid, " %s %lf %lf %lf ", c, &coord[i*3], &coord[i*3+1], &coord[i*3+2]);
    ele_num[i] = get_ele_number(c);
  }
  fclose(fid);
}

void read_xyz_out(const char *filepath, const int num_len, int *ele_num, double *coord){
  int x;
  char c[40]; //sometimes the atom names are given which can be longer than 4 chars
  FILE *fid;
  fid= fopen(filepath, "r");
  int number;
  char s[64];
  fscanf(fid,"%d", &number);
  fscanf(fid, "%s", s);
  for(int i=0; i< num_len; i++){
    fscanf(fid, " %s %lf %lf %lf ", c, &coord[i*3], &coord[i*3+1], &coord[i*3+2]);
    ele_num[i] = get_ele_number(c);
  }
  fclose(fid);
}


// generate psi4 inputfile
// m: default is 0;
// n: default is 1
//
void generate_psi4_inputfile(const char *filepath,
                             const int num_ele,
                             const double *coord,
                             const int *ele_num,
                             const int m,
                             const int n,
                             const int memory,
                             const int num_threads,
                             const char *name,
                             const char *basis,
                             const char *scf_type,
                             const char *job_type
                             )
{
  FILE *fid;
  fid = fopen(filepath, "w");
  fprintf(fid, "#psi4 generated by gen_pe\n");
  fprintf(fid, "\n");
  fprintf(fid, "memory %d mb", memory);
  fprintf(fid, "\n");
  fprintf(fid, "molecule abc_ar {\n");
  fprintf(fid, "%d %d\n", m, n);
  char ele[4];
  for(int i=0; i< num_ele; i++){
    strcpy(ele,get_ele_name(ele_num[i]));
    fprintf(fid, "%s %9.6f %9.6f %9.6f\n", ele, coord[i*3], coord[i*3+1], coord[i*3+2]);
  }
  fprintf(fid, "}\n");
  fprintf(fid, "\n");
  fprintf(fid, "set_num_threads(%d)\n", num_threads),
  fprintf(fid, "set basis %s\n", basis);
  fprintf(fid, "set guess sad\n");
  fprintf(fid, "set_max_iter=100\n");
  fprintf(fid, "set scf_type %s\n", scf_type);
  fprintf(fid, "set freeze_core True\n");
  fprintf(fid, "%s(\'%s\')\n", job_type, name);

  fclose(fid);
}

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


int main(int argc, char** argv){
    parse_cli(argc, argv);
    int num_ele=0;
    double *orig_coords=NULL;
    double *coords=NULL;
    int *ele_num= NULL;
    char output_filename[128], xyzshifted_filename[128];
    int need_to_shift_xyz_center=0;
    int output_len =0;
    char c[128];
    if(inputfile_type == 1){
        num_ele=get_file_length(xyz_filename);
        output_len=strlen(xyz_filename);
        output_len -=3;
        snprintf(psi4_filename, output_len, "%s", xyz_filename);
        strcat(psi4_filename, "_psi4.in");
        num_ele = num_ele-2;
        orig_coords=hp_malloc(num_ele*3*sizeof(double));
        coords=hp_malloc(num_ele*3*sizeof(double));
        ele_num = hp_malloc(num_ele*sizeof(int));
        read_xyz_file(xyz_filename, num_ele, ele_num, orig_coords);
        need_to_shift_xyz_center=0;
    }
    printf("num_ele is %d\n", num_ele);
    if(need_to_shift_xyz_center){
        double center_x =0.0;
        double center_y =0.0;
        double center_z =0.0;
        double sum_mass=0.0;
        double single_mass=0.0;

        for(int i=0; i<num_ele;i++){
            single_mass = an2masses[ele_num[i]];
            center_x += orig_coords[i*3]*single_mass;
            center_y += orig_coords[i*3+1]*single_mass;
            center_z += orig_coords[i*3+2]*single_mass;
            sum_mass += single_mass;
        }
        center_x /= sum_mass;
        center_y /= sum_mass;
        center_z /= sum_mass;
        for(int i=0; i< num_ele; i++){
            orig_coords[i*3] -= center_x;
            orig_coords[i*3+1] -= center_y;
            orig_coords[i*3+2] -= center_z;
        }
    }
    for(int i=0; i<num_ele;i++){
        coords[i*3] =   orig_coords[i*3];
        coords[i*3+1] = orig_coords[i*3+1];
        coords[i*3+2] = orig_coords[i*3+2];
    }
    generate_psi4_inputfile(psi4_filename,
                            num_ele,
                            coords,
                            ele_num,
                            charge,
                            spin,
                            mem,
                            num_threads,
                            "b3lyp",
                            "6-31G*",
                            "df",
                            psi4_job_type);
    hp_free(ele_num);
    hp_free(coords);
    hp_free(orig_coords);

}
