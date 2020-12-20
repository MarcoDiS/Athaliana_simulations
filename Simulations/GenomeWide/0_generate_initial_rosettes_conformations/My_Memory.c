#include<stdio.h>
#include<stdlib.h>

/* Integer */
/* Memory allocation and initialization to zero for an int vector of size dim */
int *vector1_i(int dim)
{
  /* MDS 2011 */
  int i;
  int *v;
    
  if((v = (int *)malloc((size_t) dim * sizeof(int))) == NULL) 
    {    
      fprintf(stderr, "Memory allocation problems");
      exit(1);
    } 
  
  for(i = 0; i < dim; ++i)
    {
      v[i] = 0; 
    }     
  return(v);
}
/*******************************/

/* Memory allocation for an int matrix of size dim1 * dim2 */
int **matrix2_i(int dim1, int dim2)
{
    /* MDS 2011 */
    int i;
    int **m;

    m = (int **)malloc(dim1 * sizeof(int *));
    if(m == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
        m[i] = vector1_i(dim2);
    }

    return(m);
}
/*******************************/

/* Memory allocation for an int tensor of size dim1 * dim2 * dim3 */
int ***tensor3_i(int dim1, int dim2, int dim3)
{
    /* MDS 2011 */
    int i;
    int ***t;

    t = (int ***)malloc(dim1 * sizeof(int **));
    if(t == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
      t[i] = matrix2_i(dim2, dim3);
    }

    return(t);
}
/*******************************/

/* Memory allocation for a double tensor of size dim1 * dim2 * dim3 * dim4 */
int ****tensor4_i(int dim1, int dim2, int dim3, int dim4)
{
    /* MDS 2011 */
    int i;
    int ****t;

    t = (int ****)malloc(dim1 * sizeof(int ***));
    if(t == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
      t[i] = tensor3_i(dim2, dim3, dim4);
    }

    return(t);
}
/*******************************/

/* Float */
/* Memory allocation and initialization to zero for an int vector of size dim */
float *vector1_f(int dim)
{
  /* MDS 2011 */
  int i;
  float *v;
    
  if((v = (float *)malloc((size_t) dim * sizeof(float))) == NULL) 
    {    
      fprintf(stderr, "Memory allocation problems");
      exit(1);
    } 
  
  for(i = 0; i < dim; ++i)
    {
      v[i] = 0.0; 
    }     
  return(v);
}
/*******************************/

/* Memory allocation for an int matrix of size dim1 * dim2 */
float **matrix2_f(int dim1, int dim2)
{
    /* MDS 2011 */
    int i;
    float **m;

    m = (float **)malloc(dim1 * sizeof(float *));
    if(m == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
        m[i] = vector1_f(dim2);
    }

    return(m);
}
/*******************************/

/* Memory allocation for an int tensor of size dim1 * dim2 * dim3 */
float ***tensor3_f(int dim1, int dim2, int dim3)
{
    /* MDS 2011 */
    int i;
    float ***t;

    t = (float ***)malloc(dim1 * sizeof(float **));
    if(t == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
      t[i] = matrix2_f(dim2, dim3);
    }

    return(t);
}
/*******************************/

/* Memory allocation for a double tensor of size dim1 * dim2 * dim3 * dim4 */
float ****tensor4_f(int dim1, int dim2, int dim3, int dim4)
{
    /* MDS 2011 */
    int i;
    float ****t;

    t = (float ****)malloc(dim1 * sizeof(float ***));
    if(t == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
      t[i] = tensor3_f(dim2, dim3, dim4);
    }

    return(t);
}
/*******************************/

/* Double*/
/* Memory allocation and initialization to zero for a double vector of size dim */
double *vector1_d(int dim)
{
  /* MDS 2011 */
  int i;
  double *v;
    
  if((v = (double *)malloc((size_t) dim * sizeof(double))) == NULL) 
    {    
      fprintf(stderr, "Memory allocation problems");
      exit(1);
    } 
  
  for(i = 0; i < dim; ++i)
    {
      v[i] = 0.0; 
    } 
    
  return(v);
}
/*******************************/

/* Memory allocation for a double matrix of size dim1 * dim2 */
double **matrix2_d(int dim1, int dim2)
{
    /* MDS 2011 */
    int i;
    double **m;

    m = (double **)malloc(dim1 * sizeof(double *));
    if(m == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
        m[i] = vector1_d(dim2);
    }

    return(m);
}
/*******************************/

/* Memory allocation for a double tensor of size dim1 * dim2 * dim3 */
double ***tensor3_d(int dim1, int dim2, int dim3)
{
    /* MDS 2011 */
    int i;
    double ***t;

    t = (double ***)malloc(dim1 * sizeof(double **));
    if(t == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
      t[i] = matrix2_d(dim2, dim3);
    }

    return(t);
}
/*******************************/

/* Memory allocation for a double tensor of size dim1 * dim2 * dim3 * dim4 */
double ****tensor4_d(int dim1, int dim2, int dim3, int dim4)
{
    /* MDS 2011 */
    int i;
    double ****t;

    t = (double ****)malloc(dim1 * sizeof(double ***));
    if(t == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
      t[i] = tensor3_d(dim2, dim3, dim4);
    }

    return(t);
}
/*******************************/

/* Char */
/* Memory allocation and initialization to zero for an char vector of size dim */
char *vector1_c(int dim)
{
  /* MDS 2011 */
  int i;
  char *v;
    
  if((v = (char *)malloc((size_t) dim * sizeof(char))) == NULL) 
    {    
      fprintf(stderr, "Memory allocation problems");
      exit(1);
    } 
  
  return(v);
}
/*******************************/

/* Memory allocation for a char matrix of size dim1 * dim2 */
char **matrix2_c(int dim1, int dim2)
{
    /* MDS 2011 */
    int i;
    char **m;

    m = (char **)malloc(dim1 * sizeof(char *));
    if(m == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
        m[i] = vector1_c(dim2);
    }

    return(m);
}
/*******************************/

/* Memory allocation for a char tensor of size dim1 * dim2 * dim3 */
char ***tensor3_c(int dim1, int dim2, int dim3)
{
    /* MDS 2011 */
    int i;
    char ***t;

    t = (char ***)malloc(dim1 * sizeof(char **));
    if(t == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
      t[i] = matrix2_c(dim2, dim3);
    }

    return(t);
}
/*******************************/

/* Memory allocation for a char tensor of size dim1 * dim2 * dim3 * dim4 */
char ****tensor4_c(int dim1, int dim2, int dim3, int dim4)
{
    /* MDS 2011 */
    int i;
    char ****t;

    t = (char ****)malloc(dim1 * sizeof(char ***));
    if(t == NULL) 
    {    
        fprintf(stderr, "Memory allocation problems");
        exit(1);
    }     

    for(i = 0; i < dim1; ++i)
    {
      t[i] = tensor3_c(dim2, dim3, dim4);
    }

    return(t);
}
/*******************************/

/* Free memory */

/* Integer */
void free1_i(int *v)
{
  free(v);
}
/*******************************/

void free2_i(int **m, int dim1)
{
  int i; 

  for(i = 0; i < dim1; ++i)
    {
      free1_i(m[i]);
    }

  free(m);
}
/*******************************/

void free3_i(int ***t, int dim1, int dim2)
{
  int i,j;
 
  for(i = 0; i < dim1; ++i)
    {
      free2_i(t[i], dim2);
    }
  free(t);
}
/*******************************/

void free4_i(int ****t, int dim1, int dim2, int dim3)
{
  int i,j;
 
  for(i = 0; i < dim1; ++i)
    {
      free3_i(t[i], dim2, dim3);
    }
  free(t);
}
/*******************************/

/* Float */
void free1_f(float *v)
{
  free(v);
}
/*******************************/

void free2_f(float **m, int dim1)
{
  int i; 

  for(i = 0; i < dim1; ++i)
    {
      free1_f(m[i]);
    }

  free(m);
}
/*******************************/

void free3_f(float ***t, int dim1, int dim2)
{
  int i,j;
 
  for(i = 0; i < dim1; ++i)
    {
      free2_f(t[i], dim2);
    }
  free(t);
}
/*******************************/

void free4_f(float ****t, int dim1, int dim2, int dim3)
{
  int i,j;
 
  for(i = 0; i < dim1; ++i)
    {
      free3_f(t[i], dim2, dim3);
    }
  free(t);
}
/*******************************/

/* Double */
void free1_d(double *v)
{
  free(v);
}
/*******************************/

void free2_d(double **m, int dim1)
{
  int i; 

  for(i = 0; i < dim1; ++i)
    {
      free1_d(m[i]);
    }

  free(m);
}
/*******************************/

void free3_d(double ***t, int dim1, int dim2)
{
  int i,j;
 
  for(i = 0; i < dim1; ++i)
    {
      free2_d(t[i], dim2);
    }
  free(t);
}
/*******************************/

void free4_d(double ****t, int dim1, int dim2, int dim3)
{
  int i,j;
 
  for(i = 0; i < dim1; ++i)
    {
      free3_d(t[i], dim2, dim3);
    }
  free(t);
}
/*******************************/

/* Char */
void free1_c(char *v)
{
  free(v);
}
/*******************************/

void free2_c(char **m, int dim1)
{
  int i; 

  for(i = 0; i < dim1; ++i)
    {
      free1_c(m[i]);
    }

  free(m);
}
/*******************************/

void free3_c(char ***t, int dim1, int dim2)
{
  int i,j;
 
  for(i = 0; i < dim1; ++i)
    {
      free2_c(t[i], dim2);
    }
  free(t);
}
/*******************************/

void free4_c(char ****t, int dim1, int dim2, int dim3)
{
  int i,j;
 
  for(i = 0; i < dim1; ++i)
    {
      free3_c(t[i], dim2, dim3);
    }
  free(t);
}
/*******************************/
