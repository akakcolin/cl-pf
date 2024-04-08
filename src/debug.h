#ifndef DEBUG_H
#define DEBUG_H

#include"common.h"


#ifdef __DEBUG_ON
#define debug_printf(fmt, ...) do {                            \
    printf(__CYAN "  %s:%d:%s  " fmt __RESET,                  \
      __FILE__, __LINE__, __func__ ,__VA_ARGS__);              \
  } while(0)

#define debug_puts(str) do {                                   \
    printf(__CYAN "  %s:%d:%s  " str __RESET,                  \
      __FILE__, __LINE__, __func__);                           \
  } while (0)

#define debug_print_matd33(mat) do {                           \
  debug_printf("content of %s is:\n\t%10.5f  %10.5f  %10.5f\n\t%10.5f  %10.5f  %10.5f\n\t%10.5f  %10.5f  %10.5f\n",\
    __STR(mat),                                                \
    mat[0][0], mat[0][1], mat[0][2],                           \
    mat[1][0], mat[1][1], mat[1][2],                           \
    mat[2][0], mat[2][1], mat[2][2]);                          \
  } while(0)

#define debug_print_mati33(mat) do {                           \
  debug_printf("content of %s is:\n\t%6d  %6d  %6d\n\t%6d  %6d  %6d\n\t%6d  %6d  %6d\n",\
    __STR(mat),                                                \
    mat[0][0], mat[0][1], mat[0][2],                           \
    mat[1][0], mat[1][1], mat[1][2],                           \
    mat[2][0], mat[2][1], mat[2][2]);                          \
  } while(0)

#define debug_print_vecd3(vec) do {                            \
  debug_printf("content of %s is:\n\t%10.5f  %10.5f  %10.5f\n",\
    __STR(vec),                                                \
    vec[0], vec[1], vec[2]);                                   \
  } while(0)

#define debug_print_vecdx(vec, size) do {                      \
    debug_printf("len of %s is %d , content:\n",               \
        __STR(vec), size);                                     \
    printf(__CYAN);                                            \
    for (int i=0; i!=size; ++i) {                              \
      if (i % 5 == 0) {                                        \
        printf("\t");                                          \
      }                                                        \
      printf("%8.5f", vec[i]);                                 \
      if ((i + 1) % 5 == 0) {                                  \
        printf("\n");                                          \
      }                                                        \
    }                                                          \
    printf("\n" __RESET);                                      \
  } while(0)

#define debug_print_veci3(vec) do {                            \
  debug_printf("content of %s is:\n\t%6d  %6d  %6d\n",         \
    __STR(vec),                                                \
    vec[0], vec[1], vec[2]);                                   \
  } while(0)

#define debug_print_vecix(vec, size) do {                      \
    debug_printf("len of %s is %d , content:\n",               \
        __STR(vec), size);                                     \
    printf(__CYAN);                                            \
    for (int i=0; i!=size; ++i) {                              \
      if (i % 5 == 0) {                                        \
        printf("\t");                                          \
      }                                                        \
      printf("%8d", vec[i]);                                   \
      if ((i + 1) % 5 == 0) {                                  \
        printf("\n");\
      }                                                        \
    }                                                          \
    printf("\n" __RESET);                                      \
  } while(0)

#define debug_print_vecstr(vec, size) do {                     \
  debug_printf("len of %s is %d , content:\n",                 \
      __STR(vec), size);                                       \
  printf(__CYAN);                                              \
  for(int i=0; i!=size; ++i) {                                 \
    if (i % 5 == 0) {                                          \
      printf("\t");                                            \
    }                                                          \
    printf("%20s", vec[i]);                                    \
    if ((i + 1) % 5 == 0) {                                    \
      printf("\n");                                            \
    }                                                          \
  }                                                            \
  printf("\n" __RESET);                                        \
} while(0)

#define debug_print_matdx3(mat, rows) do {                     \
    debug_printf("shape of %s is %d x 3 , content:\n",         \
        __STR(mat), rows);                                     \
    printf(__CYAN);                                            \
    for (int i=0; i!=rows; ++i) {                              \
      printf("\t%10.5f  %10.5f  %10.5f\n",                     \
          mat[i][0], mat[i][1], mat[i][2]);                    \
    }                                                          \
    printf(__RESET);                                           \
  } while(0)

#else

#define debug_puts(str)
#define debug_printf(...)
#define debug_print_matd33(mat)
#define debug_print_mati33(mat)
#define debug_print_vecd3(vec)
#define debug_print_veci3(vec)
#define debug_print_vecds(mat)
#define debug_print_matdx3(mat, rows)
#define debug_print_vecix(vec, size)
#define debug_print_vecdx(vec, size)
#define debug_print_vecstr(vec, size)

#endif
#endif
