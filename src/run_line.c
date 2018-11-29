#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <mpi.h>
#include <omp.h>

#define APP_NAME "run_line"
#define VERSION "0.2"

bool opt_debug = false;
static char input_filename[128];
static int num_lines = 0;
static int max_line_len = 258;
static int lines_allocated = 128;

#define MAX_LINE_WIDTH 512

static char const usage[] = "\
\n\
\n\
\t *** run_line for multi runing comand from comands.txt  *** \n\
\t    build :  mpicc -fopenmp run_line.c -o run_line\n\
\n\
Usage: " APP_NAME " [OPTIONS]\n\
\n\
Options:\n\
-a , --add=4   input file\n\
-h, --help\n\
Example: \n\
mpirun -np 12 run_line -a comands.txt\n\
";

static char const short_options[] = "ha:";

static struct option const options[] = {{"add", 1, NULL, 'a'},
                                        {"help", 0, NULL, 'h'}};

static void show_usage_and_exit(int status) {
  if (status)
    fprintf(stderr, "Try" APP_NAME " --help' for more information\n");
  else
    printf(usage);
  exit(status);
}

static void show_version_and_exit(void) {
  printf("%s\n", VERSION);
  exit(0);
}

int get_file_length(const char *filename) {
  FILE *fp;
  char buf[1000];
  int len = 0;
  fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "unable to open %s\n", filename);
    return 0;
  }
  while (fgets(buf, 1000, fp)) {
    len++;
  }
  fclose(fp);
  return len;
}

static void parse_arg(int key, char *arg) {
  int enum_type = 0;
  switch (key) {
  case 'a':
    num_lines = get_file_length(arg);
    strcpy(input_filename, arg);
    break;
  case 'h':
    show_usage_and_exit(0);
  default:
    show_usage_and_exit(0);
  }
}

static void parse_cli(int argc, char *argv[]) {
  int key;
  while (1) {
    key = getopt_long(argc, argv, short_options, options, NULL);
    if (key < 0)
      break;
    parse_arg(key, optarg);
  }
  if (optind < argc) {
    fprintf(stderr, "%s: unsupported non-option argument '%s'\n", argv[0],
            argv[optind]);
    show_usage_and_exit(0);
  }
}

void readline_from_file(const char *file_name, const int num_lines,
                        char **commands) {

  char line[num_lines];
  FILE *fp;
  char buf[1000];
  fp = fopen(file_name, "r");

  if (fp == NULL) {
    fprintf(stderr, "Error opening file.\n");
    exit(2);
  }
  int i;
  for (i = 0; i < num_lines; i++) {
    int j;
    if (fgets(commands[i], max_line_len - 1, fp) == NULL)
      break;
    /* Get rid of CR or LF at end of line */
    for (j = strlen(commands[i]) - 1;
         j >= 0 && (commands[i][j] == '\n' || commands[i][j] == '\r'); j--)
      ;
    commands[i][j + 1] = '\0';
  }
  fclose(fp);
}

void send_task(const char *task, const int rank) {
  int curLength = strlen(task) + 1;
  MPI_Send(&curLength, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
  MPI_Send(task, curLength, MPI_CHAR, rank, 0, MPI_COMM_WORLD);
}

int recv_task(char *single_task) {
  int str_size;
  MPI_Status status;
  MPI_Recv(&str_size, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  if (status.MPI_TAG == 1)
    return 1;
  MPI_Recv(single_task, str_size, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
           &status);
  return 0;
}

int run_task(const char *task) {
  int aux_system;
  aux_system = system(task);
  return 1;
}

int main(int argc, char **argv) {
  int omp_rank;
  int provided, required=MPI_THREAD_FUNNELED;

  MPI_Init_thread(&argc, &argv, required, &provided);

  //MPI_Init(NULL, NULL);
  int world_size;
  int world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int number;
  char **commands = NULL;
  int counter = 0;
  int totalCount = 0;
  int i, j;
  if (world_rank == 0) {
    parse_cli(argc, argv);
    if(num_lines==0){return 0;}

    /* Allocate lines of text */
    commands = (char **)malloc(sizeof(char *) * num_lines);
    for (i = 0; i < num_lines; i++) {
      commands[i] = malloc(max_line_len);
    }
    readline_from_file(input_filename, num_lines, commands);
  }
  if (world_size == 1) {
    for (i = 0; i < num_lines; i++) {
      run_task(commands[i]);
    }
    return 0;
  }

  if (world_rank == 0) {
    for (i = 1; i < world_size; i++) {
      send_task(commands[counter], i);
      counter++;
      if (counter == num_lines)
        break;
    }

    while (counter != num_lines) {
      MPI_Status status;
      int temp = 0;
      MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      send_task(commands[counter], status.MPI_SOURCE);
      counter++;
    }

    for (i = 1; i < world_size; ++i) {
      MPI_Send(0, 0, MPI_INT, i, 1, MPI_COMM_WORLD);
    }
  }
  else {
    while (true) {
      char single_task[MAX_LINE_WIDTH];
      int poison;
      poison = recv_task(single_task);
      if (poison)
        return 0;
      run_task(single_task);
      int givemework = 1;
      MPI_Send(&givemework, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
  }

  if(world_rank==0){
  for(i= num_lines-1; i>0;i--)
      free(commands[i]);
    free(commands);
  }


  MPI_Finalize();
  return 0;
}
