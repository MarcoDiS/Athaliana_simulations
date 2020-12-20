#include<stdio.h>
#include<stdlib.h>
#include"My_IO.h"
#include"My_Memory.h"

#define MAXLINE 200

double ***lammps2xyz(char *filename_input, int n_chrs, int n_beads)
{
  int i, j, n_atoms, chr, bead; 
  double min, max;
  double *xdim, *x, *ix;
  double ***coordinates;
  char line[200];
  FILE *fp_input;

  //Apertura del file di input
  fp_input = open_r(filename_input);
  
  fgets(line, MAXLINE, fp_input);
  fgets(line, MAXLINE, fp_input);
  fgets(line, MAXLINE, fp_input);

  fgets(line, MAXLINE, fp_input); 
  sscanf(line, "%d", &n_atoms); 

  /* Allocating Memory*/  
  coordinates = tensor3_d(n_chrs, n_beads, 3);
  x    = vector1_d(3);
  xdim = vector1_d(3);
  ix   = vector1_d(3);

  fgets(line, MAXLINE, fp_input);

  fgets(line, MAXLINE, fp_input);
  sscanf(line, "%lf %lf", &min, &max); 
  xdim[0] = max - min;

  fgets(line, MAXLINE, fp_input);
  sscanf(line, "%lf %lf", &min, &max); 
  xdim[1] = max - min;
  
  fgets(line, MAXLINE, fp_input);
  sscanf(line, "%lf %lf", &min, &max); 
  xdim[2] = max - min;

  fgets(line, MAXLINE, fp_input);
    
  while( fgets(line, MAXLINE, fp_input) != NULL)
    {
      sscanf(line, "%i %lf %lf %lf %lf %lf %lf", &i, &x[0], &x[1], &x[2], &ix[0], &ix[1], &ix[2]);
      chr  = (i - 1) / n_beads;
      bead = (i - 1) % n_beads; 
      for(j = 0; j < 3; j++)  coordinates[chr][bead][j] = x[j];         
    }
  fclose(fp_input);
  
  /* Free Memory */
  free1_d(x);
  free1_d(xdim);
  free1_d(ix);

  return(coordinates);   
}


