#include <stdio.h>

void int_to_string(int *p_ichar1, char *p_ichar2, int max_position )
{
    int i;
    for (i=0; i < max_position; i++){
        *p_ichar2= *(p_ichar1+i);
        p_ichar2++;
    }
    *p_ichar2= '\0';
    return;
}

void string_to_int(char *p_ichar2, int *p_ichar1, int max_position )
{
    int i;
    while ( *p_ichar2 != '\0' && i < max_position ){
        *p_ichar1= *p_ichar2;
        p_ichar1++;
        p_ichar2++;
        i++;
    }
    return;
}

/* subroutine for searhing the read in arrays for a string of defined length */

int locate_string( char *p_key, int *p_char, int num_of_chars )
{
    int i;
    char *p_start;
    p_start= p_key;
    for (i=0; i <= num_of_chars; ++i){
        if ( *(p_char+i) == *p_key ) p_key++;
        else p_key = p_start;
        if (*(p_key) == '\0') return 1;
    }
    return 0;
}
