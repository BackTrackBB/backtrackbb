/*------------------------------------------------------------*
 * IA_Ealloc
 *
 * Wrapper for malloc() with error handling
 *------------------------------------------------------------*/
#include <sys/types.h>
#include <string.h>

#include "IA_Kdiag.h"


/** @name IA_Ealloc
 * malloc() with error handling and initialization to zero
 * Dies and logs fail_msg to syslog in case of failure

 @param size       size of memory block to alloacate
 @param fail_msg   message to log in case of failure
 */
void *IA_Ealloc(size_t size, const char *fail_msg)
{
    void *p;

    if((p = (void *)calloc(size,1)) == NULL) {
        fprintf(stderr, "%s", fail_msg);
        exit(1);
    }

    return(p);
}

/* allocate a 2-d array */
void **alloc2 (size_t n1, size_t n2, size_t size)
{
    size_t i2;
    void **p;

    if ((p=(void**)malloc(n2*sizeof(void*)))==NULL)
        return NULL;
    if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
        free(p);
        return NULL;
    }
    for (i2=0; i2<n2; i2++)
        p[i2] = (char*)p[0]+size*n1*i2;
    return p;
}

/* free a 2-d array */
void free2 (void **p)
{
    free(p[0]);
    free(p);
}

/* allocate a 2-d array of floats */
float **alloc2float(size_t n1, size_t n2)
{
    return (float**)alloc2(n1,n2,sizeof(float));
}

/* free a 2-d array of floats */
void free2float(float **p)
{
    free2((void**)p);
}

/* allocate a 2-d array of doubles */
double **alloc2double(size_t n1, size_t n2)
{
    return (double**)alloc2(n1,n2,sizeof(double));
}

/* free a 2-d array of doubles */
void free2double(double **p)
{
    free2((void**)p);
}

/* allocate a 2-d array of floats, zero memory */
float **ealloc2float(size_t n1, size_t n2)
{
    float **p;

    if ((p=alloc2float(n1, n2))==NULL) {
        fprintf(stderr, "ealloc2float: malloc failed");
        exit(1);
    }
    memset(*p,0,n1*n2*sizeof(float));
    return p;
}

/* allocate a 2-d array of doubles, zero memory */
double **ealloc2double(size_t n1, size_t n2)
{
    double **p;

    if ((p=alloc2double(n1, n2))==NULL) {
        fprintf(stderr, "ealloc2double: malloc failed");
        exit(1);
    }

    memset(*p,0,n1*n2*sizeof(double));
    return p;
}
