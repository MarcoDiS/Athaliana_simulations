#include<stdio.h>
#include<stdlib.h>

int *vector1_i(int dim);
int **matrix2_i(int dim1, int dim2);
int ***tensor3_i(int dim1, int dim2, int dim3);
int ****tensor4_i(int dim1, int dim2, int dim3, int dim4);
float *vector1_f(int dim);
float **matrix2_f(int dim1, int dim2);
float ***tensor3_f(int dim1, int dim2, int dim3);
float ****tensor4_f(int dim1, int dim2, int dim3, int dim4);
double *vector1_d(int dim);
double **matrix2_d(int dim1, int dim2);
double ***tensor3_d(int dim1, int dim2, int dim3);
double ****tensor4_d(int dim1, int dim2, int dim3, int dim4);
void free1_i(int *v);
void free2_i(int **m, int dim1);
void free3_i(int ***t, int dim1, int dim2);
void free4_i(int ****t, int dim1, int dim2, int dim3);
void free1_f(float *v);
void free2_f(float **m, int dim1);
void free3_f(float ***t, int dim1, int dim2);
void free4_f(float ****t, int dim1, int dim2, int dim3);
void free1_d(double *v);
void free2_d(double **m, int dim1);
void free3_d(double ***t, int dim1, int dim2);
void free4_d(double ****t, int dim1, int dim2, int dim3);

