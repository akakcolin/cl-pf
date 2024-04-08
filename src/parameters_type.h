#ifndef PARAMETERS_TYPE_H
#define PARAMETERS_TYPE_H

#include "defines.h"

typedef struct s_input_parameters {
  char local_execute[MAX_PATH];
  char config_vina[MAX_PATH_FILE_NAME];
  char log_vina[MAX_PATH];
  char out_vina[MAX_PATH];
  char script_ligand4[MAX_PATH_FILE_NAME];
  char script_receptor4[MAX_PATH_FILE_NAME];
  char in_receptor[MAX_PATH];
  char out_receptor[MAX_PATH];
  char in_compounds[MAX_PATH];
  char out_ligands[MAX_PATH];
  char options_vina[MAX_PATH];
  char options_ligand4[MAX_PATH];
  char options_receptor4[MAX_PATH];
} input_parameters_t;

#endif
