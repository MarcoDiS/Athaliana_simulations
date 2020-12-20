#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include"My_Memory.h"
#include"My_Memory.c"
#include"My_IO.h"
#include"My_IO.c" 

#define RADIUS 12.0
#define PARTICLE_RADIUS 0.5
#define MAXLINE 50

long int idum= (-1);

/* Functions */
int check_point_inside_the_confining_environment(double *P, double a, double b, double c);
int check_point_inside_the_nucleolus_environment(double *P, double a, double b, double c);
int check_clashes(double **Seg1_P1, double **Seg1_P0, int N_CHR);
float ran2 (long *idum);
void rods(double **Seg1_P1, double **Seg1_P0, double *length, int N_CHR, double a, double b, double c, double iradius, double eradius);
void rosette(double ***a, double *length, int n_chr, int n_beads);
void rosettes_rototranslation(double ***a, double **Seg1_P1, double **Seg1_P0, int *n_beads, int N_CHR);
double scal_d (double *a, double *b, int dim);
double norm_d (double *a, int dim);
double distance_between_segments_3d(double *Seg1_P1, double *Seg1_P0, double *Seg2_P1, double *Seg2_P0);
double distance(double *P0, double *P1, int dim);
int check_overlap_with_beads(double *P, double ***Chr, int n_chr, int *n_beads, int dim, double size);
int check_bead_inside(double ***Chr, int n_chr, int *n_beads, int dim, double a, double b, double c);
int check_bead_overlap(double ***Chr, int n_chr, int *n_beads, int dim, double size);
void write_input_file(double ***Chr, int n_atoms, int *n_beads, int *offset, int N_CHR, double A, double B, double C, double **P, int N_NUCL);
void print_help();

int main(int argc,char** argv)
{
  int c, i, j, k, flag, cnt, n_atoms, N_CHR, N_NUCL, MAX_BEADS, SEED;
  double A, B, C, CNT, I, E;
  int *n_beads, *offset;
  long int t;
  double l, *length, r, dist;
  double **Seg1_P1, **Seg1_P0, ***Chr, **P;
  char chr[50][10], line[MAXLINE];
  FILE *fp_input;
  
  if(argc < 1) 
    { 
      print_help(); 
    }
  /* Get options executing the file */
  while ((c = getopt (argc, argv, "h:a:b:c:n:s:m:i:e:")) != -1)
    switch (c)
      {
      case 'e':
	E = atof(optarg);
	fprintf(stderr, "External radius: %lf\n", A);
	break;
      case 'i':
	I = atof(optarg);
	fprintf(stderr, "Internal radius: %lf\n", A);
	break;
      case 'a':
	A = atof(optarg);
	fprintf(stderr, "x-axis of the ellipse of the confining environment: %lf\n", A);
	break;
      case 'b':
	B = atof(optarg);
	fprintf(stderr, "y-axis of the ellipse of the confining environment: %lf\n", B);
	break;
      case 'c':
	C = atof(optarg);
	fprintf(stderr, "z-axis of the ellipse of the confining environment: %lf\n", C);
	break;
      case 'n':
	N_NUCL = atof(optarg);
	fprintf(stderr, "Number of nucleoli: %d\n", N_NUCL);
	break;
      case 's':
	SEED = atoi(optarg);
	fprintf(stderr, "Number of random numbers to discard: %d\n", SEED);
	break;
      case 'm':
	CNT = atof(optarg);
	fprintf(stderr, "Number of desired conformations: %lf\n", CNT);
	break;
      case '?':
	if (optopt == 'c')
	  fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	else if (isprint (optopt))
	  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	else
	  fprintf (stderr,
		   "Unknown option character `\\x%x'.\n",
		   optopt);
	return 1;
      default:
	abort ();
      }
  
  /* Convert axes in semi-axes */
  A = A * 0.5 - 1.0;
  B = B * 0.5 - 1.0;
  C = C * 0.5 - 1.0; 
  
  fprintf(stderr, "Radius of the chromosomes: %lf\n\n", RADIUS);
  
  /* Opening the input files: its format should be %s (string of the chr name) %d (number of beads) */
  fp_input = open_r("INPUT_PARAMETERS.DAT");
  
  /* Getting the number of chromosomes in the system */
  N_CHR = 0;
  fprintf(stderr, "Getting the number of chromosomes in the system\n");
  while( fgets(line, MAXLINE, fp_input) != NULL ) N_CHR++;
  fprintf(stderr, "Number of chromosomes: %d\n", N_CHR);
  
  /* Memory Allocation 1 */ 
  n_beads   = vector1_i(N_CHR);     /* This array  of int    contains the number of beads for each chromosome */
  offset    = vector1_i(N_CHR);     /* This array  of int    contains the offset for each chromosome */
  length    = vector1_d(N_CHR);     /* This array  of double contains the length of the Model rosette chromosomes */
  Seg1_P1   = matrix2_d(N_CHR, 3);  /* This matrix of double contains the coordinates of initial point of the segments mimicking the Model rosette chromosomes */
  Seg1_P0   = matrix2_d(N_CHR, 3);  /* This matrix of double contains the coordinates of final   point of the segments mimicking the Model rosette chromosomes */
  P         = matrix2_d(N_NUCL, 3); /* This matrix of double contains the coordinates of the nucleoli */
  
  /* Getting the maximum number of beads per chromosome in the system */
  fprintf(stderr, "Getting the maximum number of beads per chromosome in the system\n");
  rewind(fp_input);                /* This is to read from the beginning the input file */
  N_CHR = MAX_BEADS = n_atoms = 0; /* Initialising the variables containing the chr counter, the max number of beads and the total number of atoms respectively */
  while( fgets(line, MAXLINE, fp_input) != NULL)
    {
      sscanf(line, "%s %d", &chr[N_CHR], &n_beads[N_CHR]);       /* reading each line of the input file */
      if(n_beads[N_CHR] > MAX_BEADS) MAX_BEADS = n_beads[N_CHR]; /* checking the maximum length of the chromosomes */
      
      offset[N_CHR] =  n_atoms;                                  /* storing the offset of the chr offset */
      n_atoms       += n_beads[N_CHR];                           /* updating the number of beads in the system */
      
      N_CHR++;                                                   /* Increasing the chr counter */
    }
  fprintf(stderr, "Maximum number of beads: %d\n", MAX_BEADS); 
  
  /* Memory Allocation 2 */ 
  Chr = tensor3_d(N_CHR, MAX_BEADS, 3);                          /* This array of double contains the coordinates of initial Model rosette chromosomes */
  
  /* Closing the input file */
  fclose(fp_input);
  
  /* Construction of the model rosette chromosomal structure one by one to know the length of the segments */
  for(i = 0; i < N_CHR; i++)
    {
      rosette(Chr, &l, i, n_beads[i]);    
      length[i] = l;
      fprintf(stderr, "%d %lf\n",i,length[i]);
    } 
  fprintf(stderr, "\n\n");
  
  /* Constructing the initial chromosome random conformations */    
  //  fprintf(stderr, "Constructing the initial chromosome random conformations\n");
  flag = cnt = t = 0;
  for(i=0;i<SEED;i++) j=ran2(&idum);
  while(cnt < CNT)
    {
    begin:
      t++;
      /* Guess of the initial rod conformation: */
      /* 1 - each rod is choosen in at random and placed inside the confining evironment in a random position and with random orientation */
      /* 2 - possible clashes between generated rods are checked */
      fprintf(stderr, "Generating the initial conformation of the rods\n");
      rods(Seg1_P1, Seg1_P0, length, N_CHR, A, B, C, I, E);
      
      /* Building the Model rosette chromosomal initial conformation */
      fprintf(stderr, "Building the Model rosette chromosomal initial conformation\n");
      for(i = 0; i < N_CHR; i++) rosette(Chr, &l, i, n_beads[i]);    
      
      /* Roto-translation of the Model rosette chromosomal initial conformation according to the segment position and orientation */
      fprintf(stderr, "Roto-translation of the Model rosette chromosomal initial conformation according to the segment position and orientation\n");
      rosettes_rototranslation(Chr, Seg1_P1, Seg1_P0, n_beads, N_CHR);
      
      /* Checking that the beads are all inside the confining environment */
      fprintf(stderr, "Checking if the beads are inside the confining environment\n");
      if(check_bead_inside(Chr, N_CHR, n_beads, 3, A, B, C) == 1) goto begin;      
      
      /* Checking that chromosome not bearing NORs are outside the central sphere */
      for(i = 0; i < N_CHR; i++) 
	{
	  if ( i==2 || i==3 || i==6 || i==7 ) continue ;
	  //chr1A 30427671 10143
	  //chr1B 30427671 10143
	  //chr2A 23298289  7767
	  //chr2B 23298289  7767
	  //chr3A 23459830  7820
	  //chr3B 23459830  7820
	  //chr4A 22585056  7529
	  //chr4B 22585056  7529
	  //chr5A 26975502  8992
	  //chr5B 26975502  8992
	  for (j = 0 ; j < n_beads[i]; j++)
	    {
	      if(check_point_inside_the_confining_environment(Chr[i][j], E, E, E) == 1) goto begin;
	    }
	}
      
      /* Checking that the beads are not overlapping: 1 -> overlap -> go to begin; 0 -> no-overlap -> store the conformation */
      fprintf(stderr, "Checking overlaps between the beads\n");
      if(check_bead_overlap(Chr, N_CHR, n_beads, 3, PARTICLE_RADIUS*2.0) == 1) goto begin;      
      
      /* increase the number of accepted initial conformations */
      cnt++;
      
      /* Writing the final rod conformation */
      fprintf(stderr, "Succesfully generated configuration number %d after %ld attempts\n",cnt,t);
      for(i = 0; i < N_CHR; ++i)
	{
	  fprintf(stderr, "%lf %lf %lf\n", Seg1_P0[i][0],Seg1_P0[i][1],Seg1_P0[i][2]);
	  fprintf(stderr, "%lf %lf %lf\n\n\n", Seg1_P1[i][0],Seg1_P1[i][1],Seg1_P1[i][2]);
	}
            
      /* OUTPUT on stdout */
      write_input_file(Chr, n_atoms, n_beads, offset, N_CHR, A, B, C, P, N_NUCL);
      fflush(stdout);
    }  
  
  /* Free Memory */
  free1_i(n_beads);
  free1_i(offset);    
  free1_d(length);  
  free2_d(Seg1_P1, N_CHR);
  free2_d(Seg1_P0, N_CHR);
  free3_d(Chr, N_CHR, MAX_BEADS); 
  
  return(0);
}


/* Functions */

void print_help()
{
  fprintf(stderr," \n");
  fprintf(stderr,"OPTIONS: \n");
  fprintf(stderr,"\t-h:\t\t print this help and exit.\n");
  fprintf(stderr,"\t-a:\t\t x-axis of the basal ellipse of the confining environment");
  fprintf(stderr,"\t-b:\t\t y-axis of the basal ellipse of the confining environment");
  fprintf(stderr,"\t-c:\t\t z-axis of the confining environment");
  fprintf(stderr,"\t-n:\t\t Number of nucleoli\n");
  fprintf(stderr,"\t-m:\t\t Number of desired conformations");
  fprintf(stderr," \n");
  exit(1);
}

/*************************/

/* Random number generator */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2 (long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  float temp;


  if (*idum <= 0)
    {
      if (-(*idum) < 1)
	*idum = 1;
      else
	*idum = -(*idum);
      idum2 = (*idum);
      for (j = NTAB + 7; j >= 0; j--)
	{
	  k = (*idum) / IQ1;
	  *idum = IA1 * (*idum - k * IQ1) - k * IR1;
	  if (*idum < 0)
	    *idum += IM1;
	  if (j < NTAB)
	    iv[j] = *idum;
	}
      iy = iv[0];
    }
  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k * IQ1) - k * IR1;
  if (*idum < 0)
    *idum += IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if (idum2 < 0)
    idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1)
    iy += IMM1;
  if ((temp = AM * iy) > RNMX)
    return RNMX;
  else
    return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/*******************************/

/* Construction of the rods initial conformation
void rods(double **Seg1_P1, double **Seg1_P0, double *length, int N_CHR, double a, double b, double c)
{
  int    i, j, k, l, m, flag, cnt, tentative;
  int    order[N_CHR];
  double x, y, z, temp_theta, temp_phi, rescaled_dist_orig;
  double dist0, dist1;

  //  fprintf(stderr, "%lf %lf %lf\n", a,b,c); 


 begin:
      
   for(i = 0; i < N_CHR; i++)
    {     
      j = i;

      tentative = 0;
    position:      
      if((tentative++) > 100000)
	{
	  fprintf(stderr,"Too many tentative: I will start over!\n");
	  goto begin;
	}
      fprintf(stderr, "Trying to position terminus 0\n");
      fprintf(stderr, "Length = %lf\n", length[j]);
       do
      {
      //Seg1_P0[j][0] = (2.0*ran2(&idum)-1.0)*a*0.5;
      //#Seg1_P0[j][1] = 140. + ran2(&idum) * b    ;
      //#Seg1_P0[j][2] = (2.0*ran2(&idum)-1.0)*c*0.5;
	Seg1_P0[j][0] = (2.0*ran2(&idum)-1.0)*(a-RADIUS);
	Seg1_P0[j][1] = (2.0*ran2(&idum)-1.0)*(b-RADIUS);
	Seg1_P0[j][2] = (2.0*ran2(&idum)-1.0)*(c-RADIUS);
      } while (check_point_inside_the_confining_environment(Seg1_P0[j], a, b, c) == 1);
      fprintf(stderr,"%lf %lf %lf\n", Seg1_P0[j][0], Seg1_P0[j][1], Seg1_P0[j][2]);
      
      //      fprintf(stderr, "Successfully positioned terminus 0 of chromosome %d\n", j);
      //      fprintf(stderr, "Trying to position terminus 1\n");
      //
      temp_theta    = acos(2.0*ran2(&idum)-1.0);
      temp_phi      = 2*acos(-1.0)*ran2(&idum);
      Seg1_P1[j][0] = Seg1_P0[j][0] + length[j] * cos(temp_phi) * sin(temp_theta);
      Seg1_P1[j][1] = Seg1_P0[j][1] + length[j] * sin(temp_phi) * sin(temp_theta);
      Seg1_P1[j][2] = Seg1_P0[j][2] + length[j] * cos(temp_theta);      
      //
      if(check_point_inside_the_confining_environment(Seg1_P1[j], a, b, c) == 1) goto position;
      fprintf(stderr, "Successfully positioned terminus 1 of chromosome %d\n", j);
      fprintf(stderr,"%lf %lf %lf\n", Seg1_P1[j][0], Seg1_P1[j][1], Seg1_P1[j][2]);
      //
      if(check_clashes(Seg1_P1, Seg1_P0, i+1) == 1) goto position;

      //      fprintf(stderr, "Successfully positioned chromosome %d\n", j);
    }   

  return;
}
*/

/* Construction of the rods initial conformation */
void rods(double **Seg1_P1, double **Seg1_P0, double *length, int N_CHR, double a, double b, double c, double iradius, double eradius)
{
  int    i, j, k, l, m, flag, cnt, tentative;
  int    order[N_CHR];
  double x, y, z, temp_theta, temp_phi, rescaled_dist_orig;
  double dist0, dist1;

  //  fprintf(stderr, "%lf %lf %lf\n", a,b,c); 


 begin:
  
  for(i = 0; i < N_CHR; i++)
    {     
      j = i;
      
      tentative = 0;
    position:      
      if((tentative++) > 100000)
	{
	  fprintf(stderr,"Too many tentative: I will start over!\n");
	  goto begin;
	}
      fprintf(stderr, "Trying to position terminus 0\n");
      fprintf(stderr, "Length = %lf\n", length[j]);
      do
	{
	  if ( j==2 || j==3 || j==6 || j==7 )
	    {
	      Seg1_P0[j][0] = (2.0*ran2(&idum)-1.0)*iradius;
	      Seg1_P0[j][1] = (2.0*ran2(&idum)-1.0)*iradius;
	      Seg1_P0[j][2] = (2.0*ran2(&idum)-1.0)*iradius;
	    } else {
	    do 
	      {
		Seg1_P0[j][0] = (2.0*ran2(&idum)-1.0)*(a-RADIUS);
		Seg1_P0[j][1] = (2.0*ran2(&idum)-1.0)*(b-RADIUS);
		Seg1_P0[j][2] = (2.0*ran2(&idum)-1.0)*(c-RADIUS);	      
	      }  while (check_point_inside_the_confining_environment(Seg1_P0[j], eradius, eradius, eradius) == 0);
	  }
	} while (check_point_inside_the_confining_environment(Seg1_P0[j], a, b, c) == 1);
      fprintf(stderr,"%lf %lf %lf\n", Seg1_P0[j][0], Seg1_P0[j][1], Seg1_P0[j][2]);
      
      //      fprintf(stderr, "Successfully positioned terminus 0 of chromosome %d\n", j);
      //      fprintf(stderr, "Trying to position terminus 1\n");
      
      temp_theta    = acos(2.0*ran2(&idum)-1.0);
      temp_phi      = 2*acos(-1.0)*ran2(&idum);
      Seg1_P1[j][0] = Seg1_P0[j][0] + length[j] * cos(temp_phi) * sin(temp_theta);
      Seg1_P1[j][1] = Seg1_P0[j][1] + length[j] * sin(temp_phi) * sin(temp_theta);
      Seg1_P1[j][2] = Seg1_P0[j][2] + length[j] * cos(temp_theta);      
      //
      if(check_point_inside_the_confining_environment(Seg1_P1[j], a, b, c) == 1) goto position;
      fprintf(stderr, "Successfully positioned terminus 1 of chromosome %d\n", j);
      fprintf(stderr,"%lf %lf %lf\n", Seg1_P1[j][0], Seg1_P1[j][1], Seg1_P1[j][2]);
      //
      if(check_clashes(Seg1_P1, Seg1_P0, i+1) == 1) goto position;
      
      //      fprintf(stderr, "Successfully positioned chromosome %d\n", j);
    }   
  
  return;
}

 /*
void rods(double **Seg1_P1, double **Seg1_P0, double *length, int N_CHR, double a, double b, double c)
{
  int    i, j, k, l, m, flag, cnt, tentative;
  int  order[N_CHR];
  double x, y, z, temp_theta, temp_phi, rescaled_dist_orig;
  double dist0, dist1;

  for (i = 0; i < N_CHR; i++)
    {
      do
	{
	  cnt = 0;
	  m = round((ran2(&idum) * (N_CHR-1))) ;
	  order[i] = m;
	  for(j=0;j<i;j++) if(order[j]==m) cnt++;
	} while (cnt == 1);
      fprintf(stderr, "%d %d %lf\n",i,order[i], (float) order[i] * 2. * RADIUS - (float) N_CHR * RADIUS + ((float)RADIUS));
    }

  for(i = 0; i < N_CHR; i++)
    {     

      j = i;

      tentative = 0;
      
      Seg1_P0[j][0] = (2.0*ran2(&idum)-1.0)*a*0.5;
      Seg1_P0[j][1] = 140. + ran2(&idum) * b    ;
      Seg1_P0[j][2] = (2.0*ran2(&idum)-1.0)*c*0.5;
      
      if ( tentative == 1)
	{
	  fprintf(stderr,"%lf %lf %lf\n", Seg1_P0[j][0], Seg1_P0[j][1], Seg1_P0[j][2]);
	}
      temp_theta    = acos(2.0*ran2(&idum)-1.0);
      temp_phi      = 2*acos(-1.0)*ran2(&idum);
      Seg1_P1[j][0] = Seg1_P0[j][0] + length[j] * cos(temp_phi) * sin(temp_theta);
      Seg1_P1[j][1] = Seg1_P0[j][1] + length[j] * sin(temp_phi) * sin(temp_theta);
      Seg1_P1[j][2] = Seg1_P0[j][2] + length[j] * cos(temp_theta);      

      fprintf(stderr,"%lf %lf %lf\n", Seg1_P1[j][0], Seg1_P1[j][1], Seg1_P1[j][2]);

    }   

  return;
}
*/

/*******************************/

/* Distance inside the ellipsoid 
double distance_inside_the_ellipsoid(double *P, double a, double b, double c)
{
  return(sqrt((P[0]*P[0])/(a*a)+(P[1]*P[1])/(b*b)+(P[2]*P[2])/(c*c)));
  } */

/*******************************/

/* Check if the point is inside the confining environment */
int check_point_inside_the_confining_environment(double *P, double a, double b, double c)
{
  /*  if( (((P[0]*P[0])/(a*a) + (P[1]*P[1])/(b*b))) > 1.0) return(1);
      if( (P[2]*P[2]) > (c*c) ) return(1); */

  /* For a sphere or ellipsoid */
  if( (((P[0]*P[0])/(a*a) + (P[1]*P[1])/(b*b) + (P[2]*P[2])/(c*c))) > 1.0) return(1); 
  
  return(0);
} 

/* Check if the point is inside the confining environment */
int check_point_inside_the_nucleolus_environment(double *P, double a, double b, double c)
{
  /*  if( (((P[0]*P[0])/(a*a) + (P[1]*P[1])/(b*b))) > 1.0) return(1);
      if( (P[2]*P[2]) > (c*c) ) return(1); */

  /* For a sphere or ellipsoid */
  if( (((P[0]*P[0])/(a*a) + (P[2]*P[2])/(c*c))) > 1.0) return(1); 
  
  return(0);
} 

/*******************************/

/* Check for PBC clashes within a sphere */
int check_clashes(double **Seg1_P1, double **Seg1_P0, int N_CHR)
{
  int i, j, k, flag;
  double dist;
  double AuxSeg_P1[3], AuxSeg_P0[3];

  flag = 1;

  for(i = 0; i < N_CHR; i++)
    {
      for(k = 0; k < 3; k++) AuxSeg_P0[k] = Seg1_P0[i][k];
      for(k = 0; k < 3; k++) AuxSeg_P1[k] = Seg1_P1[i][k];

      /* Check for steric slashes */
      for(j = i+1; j < N_CHR; j++)
	{    
	  /* Comparison between a Segment and All the Segments in the system */
	  if((dist = distance_between_segments_3d(Seg1_P1[j], Seg1_P0[j], AuxSeg_P1, AuxSeg_P0)) < (2.0 * RADIUS))
	    {
	      fprintf(stderr, "Clash between segments %3d and %3d at distance: %lf\n", i, j, dist); 
	      return(flag);
	    }
	}
    }

  flag = 0;
  return(flag);
}

/*******************************/

/* Construction of the model rosette chromosomal initial conformation */
void rosette(double ***a, double *length, int chr, int n_beads)
{
  int j;
  double phi, dist, x, y, z, x_0, y_0, z_0;

  phi = 0.0;
  j   = 0;
  /* First bead position */
  x_0 = RADIUS * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * cos(phi);
  y_0 = RADIUS * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * sin(phi);
  z_0 = phi / (2.0 * 3.141592654);   
  a[chr][j][0] = x_0; 
  a[chr][j][1] = y_0; 
  a[chr][j][2] = z_0;
  //  if(j == 0) fprintf(stderr,"First bead for Chr %d in position: %lf %lf %lf\n", chr, a[chr][j][0], a[chr][j][1], a[chr][j][2]);
  /* Growing the chain */
  for(j = 1; j < n_beads; j++)
    {
      dist = 0.0;
      do
	{   
	  phi += 0.001;
	  x = RADIUS * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * cos(phi);
	  y = RADIUS * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * sin(phi);
	  z = phi / (2.0 * 3.141592654);                                  
	  dist = sqrt((x - x_0)*(x - x_0) + (y - y_0)*(y - y_0) + (z - z_0)*(z - z_0));
	} while(dist < (PARTICLE_RADIUS*2.0));
      //      fprintf(stderr,"particle %d value of phi %lf\n",j,phi);
      x_0 = x;
      y_0 = y;                    
      z_0 = z;
      a[chr][j][0] = x; 
      a[chr][j][1] = y; 
      a[chr][j][2] = z;
      if(dist > (PARTICLE_RADIUS*2.0)*1.2) fprintf(stderr, "%lf %d %d %d\n", dist, j-1, j, chr);

      //      if(j == (n_beads-1)) fprintf(stderr,"Last bead for Chr %d in position: %lf %lf %lf\n", chr, a[chr][j][0], a[chr][j][1], a[chr][j][2]);
    }
  *length = a[chr][n_beads-1][2]-a[chr][0][2];

  return;
}  

/*******************************/

void rosettes_rototranslation(double ***Chr, double **Seg1_P1, double **Seg1_P0, int *n_beads, int N_CHR)
{
  int i, j;
  double vector[3], theta[3];
  double x_temp_2, y_temp_2, z_temp_2; 
  double x_temp_1, y_temp_1, z_temp_1; 
  double x, y, z; 

  for(i = 0; i < N_CHR; i++)
    {
      for(j = 0; j < 3; j++) vector[j] = (Seg1_P1[i][j] - Seg1_P0[i][j]);
      
      /* Rotation Angles*/       
      theta[0] = atan2(vector[1],vector[2]);
      
      x_temp_2 =  vector[0];
      y_temp_2 =  cos(theta[0]) * vector[1] - sin(theta[0]) * vector[2];                
      z_temp_2 =  sin(theta[0]) * vector[1] + cos(theta[0]) * vector[2];
      
      theta[1] = atan2(x_temp_2,z_temp_2);
      
      x_temp_1 =  cos(theta[1]) * x_temp_2 - sin(theta[1]) * z_temp_2;
      y_temp_1 =  y_temp_2;
      z_temp_1 =  sin(theta[1]) * x_temp_2 + cos(theta[1]) * z_temp_2;
      
      if(z_temp_1 < 0.0)
        {
	  z_temp_1 = -z_temp_1;
	  theta[2] = 3.141592654;
        }
      else
        {
	  theta[2] = 0.0;    
        }
      
      /* Chromosome roto-translations */
      for(j = 0; j < n_beads[i]; ++j)    
        {
	  
	  x_temp_2 =  Chr[i][j][0];
	  y_temp_2 =  cos(theta[2]) * Chr[i][j][1] + sin(theta[2]) * Chr[i][j][2];
	  z_temp_2 = -sin(theta[2]) * Chr[i][j][1] + cos(theta[2]) * Chr[i][j][2];
          
	  x_temp_1 =  cos(theta[1]) * x_temp_2 + sin(theta[1]) * z_temp_2;
	  y_temp_1 =  y_temp_2;
	  z_temp_1 = -sin(theta[1]) * x_temp_2 + cos(theta[1]) * z_temp_2;
	  
	  x =  x_temp_1;
	  y =  cos(theta[0]) * y_temp_1 + sin(theta[0]) * z_temp_1;
	  z = -sin(theta[0]) * y_temp_1 + cos(theta[0]) * z_temp_1;
	  
	  /* Chromosome translations */
	  Chr[i][j][0] = Seg1_P0[i][0] + x;
	  Chr[i][j][1] = Seg1_P0[i][1] + y;
	  Chr[i][j][2] = Seg1_P0[i][2] + z;
        }
    }    

  return;
}

/*******************************/

/* Funzione per il prodotto scalare di due vettori di lunghezza dim*/
double scal_d (double *a, double *b, int dim)
{
  /* CM 2006 */
  int i;
  double temp;
  
  temp = 0.0;
  for (i = 0; i < dim; i++)
    {
      temp += a[i] * b[i];
    }
  return (temp);
}

/*******************************/

/* Funzione per il calcolo della norma di un vettore di lunghezza dim*/
double norm_d (double *a, int dim)
{
    return (sqrt (scal_d (a, a, dim)));
}

/*******************************/
/* finds the minimum distance between two segments */
/* Seg1_P1 etc. are arrays with three components
   e.g. double Seg_P1[3]; etc. */
   
double distance_between_segments_3d(double *Seg1_P1, double *Seg1_P0, double *Seg2_P1, double *Seg2_P0)
{
    /*Inspiration: http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm 
    Copyright 2001, softSurfer (www.softsurfer.com)
    This code may be freely used and modified for any purpose
    providing that this copyright notice is included with it.
    SoftSurfer makes no warranty for this code, and cannot be held
    liable for any real or imagined damage resulting from its use.
    Users of this code must verify correctness for their application.
    */
    double u[3], v[3], w[3], dP[3], a, b, c, d, e, D, sc, sN, sD, tc, tN, tD;
    int i;

    for(i = 0; i < 3; i++)
    {
        /* vectors for each segment */  
        u[i]=   Seg1_P1[i]-Seg1_P0[i];
        v[i]=   Seg2_P1[i]-Seg2_P0[i];
        w[i]=   Seg1_P0[i]-Seg2_P0[i];
    }       

    a=scal_d(u,u,3);        
    b=scal_d(u,v,3);        
    c=scal_d(v,v,3);        
    d=scal_d(u,w,3);        
    e=scal_d(v,w,3);        
    D  = a*c - b*b;
    sD = tD = D;

    if (D< (1.0e-7)) 
    { /* segments almost parallel */
        sN = 0.0;
        sD = 1.0;
        tN = e;
        tD = c;
    }  
    else 
    { /* get the closest points on the infinite lines */
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) 
        {   /* sc < 0 => the s=0 edge is visible */
            sN = 0.0;
            tN = e;
            tD = c;
        }

        else if (sN > sD) {  /* sc > 1 => the s=1 edge is visible */
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

  if (tN < 0.0) {           /* tc < 0 => the t=0 edge is visible */
        tN = 0.0;
        /* recompute sc for this edge */
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      /* tc > 1 => the t=1 edge is visible */
        tN = tD;
        /* recompute sc for this edge */
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);    
            sD = a;
        }
    }
    /* finally do the division to get sc and tc */
    sc = (fabs(sN) < (1.0e-7) ? 0.0 : sN / sD);
    tc = (fabs(tN) < (1.0e-7) ? 0.0 : tN / tD);
     
    /* get the difference of the two closest points */
    
    for(i=0; i <3; i++)
    {
        dP[i]=w[i]+(sc*u[i])-(tc*v[i]); /* = S1(sc) - S2(tc) */
    }       

    return (norm_d(dP,3));   /* return the closest distance */
}

/*******************************/

void write_input_file(double ***Chr, int n_atoms, int *n_beads, int *offset, int N_CHR, double A, double B, double C, double **P, int N_NUCL)
{
  int i, j, k;

  fprintf(stdout, "LAMMPS input data file \n\n");
  fprintf(stdout, "%9d atoms \n" , n_atoms+N_NUCL );
  fprintf(stdout, "%9d bonds \n" , n_atoms - N_CHR );
  fprintf(stdout, "%9d angles \n\n" , 79416); //n_atoms - 2*N_CHR );
  fprintf(stdout, "%9s atom types \n" , "50" );
  fprintf(stdout, "%9s bond types \n" , "3" );
  fprintf(stdout, "%9s angle types \n\n" , "1" );
  fprintf(stdout, "%6.3lf    %6.3lf     %s %s \n", -(A+3.0), A+3.0, "xlo", "xhi");
  fprintf(stdout, "%6.3lf    %6.3lf     %s %s \n", -(B+3.0), B+3.0, "ylo", "yhi");
  fprintf(stdout, "%6.3lf    %6.3lf     %s %s \n", -(C+3.0), C+3.0, "zlo", "zhi");
  
  fprintf(stdout, "\n Atoms \n\n");              
  k=0;
  for(i = 0; i < N_CHR; i++)
    {
      for(j = 0; j < n_beads[i]; j++)    
        {
	  fprintf(stdout, "%5d %s %s %7.4lf %7.4lf %7.4lf \n", offset[i] + j + 1, "1", "1", Chr[i][j][0], Chr[i][j][1], Chr[i][j][2]);
	  k++;
        }
    }
  /*
  for(i = 0; i < N_NUCL; i++)
    {
      k++;
      fprintf(stdout, "%5d %s %s %7.4lf %7.4lf %7.4lf \n", k, "1", "1", P[i][0], P[i][1], P[i][2]);
    }
 
  fprintf(stdout, "\n Bonds \n\n");
  k = 0;
  for(i = 0; i < N_CHR; i++)
    {
      for(j = 0; j < (n_beads[i] - 1) ; j++)        
        {
	  k++; 
	  fprintf(stdout, " %4d %s %4d %4d \n", k, "1", k + i, k + i + 1);           
        }    
    }
  
  fprintf(stdout, "\n Angles \n\n");
  k = 0;
  for(i = 0; i < N_CHR; i++)
    {
      for(j = 0; j < (n_beads[i] - 2); j++)        
        {
	  k++;
	  fprintf(stdout, "%4d %s %5d %5d %5d \n", k, "1", k + (i * 2), k + (i * 2) + 1, k + (i * 2) + 2);           
        }    
    }
  */

  return;
}

/*******************************/

double distance(double *P0, double *P1, int dim)
{
  int k;
  double distance;
  
  distance = 0.0;
  for(k = 0; k < dim; k++)
    {
      distance += (P0[k] - P1[k])*(P0[k] - P1[k]);
    }   
  distance = sqrt(distance);
  
  return(distance);
} 

/*******************************/

int check_bead_overlap(double ***Chr, int n_chr, int *n_beads, int dim, double size)
{
  int i,j; /* Indeces of the 3D tensor: i is chr1, j is for the beads in chr1 */
  int l,m; /* Indeces of the 3D tensor: l is chr2, m is for the beads in chr2 */
  
  for(i = 0; i < n_chr; i++)
    {
      for(l = i+1; l < n_chr; l++)
	{
	  for(j = 0; j < n_beads[i]; j++)
	    {
	      for(m = 0; m < n_beads[l]; m++)
		{ 
		  if( (distance(Chr[i][j], Chr[l][m], 3)) < size )
		    {
		      return(1);
		      //		      fprintf(stderr,"Bead %d of chr %d in position (%lf %lf %lf) and bead %d of chr %d position (%lf %lf %lf) are overlapping\n", j, i, Chr[i][j][0], Chr[i][j][1], Chr[i][j][2], m, l, Chr[l][m][0], Chr[l][m][1], Chr[l][m][2]);
		    }
		}	      
	    }
	}	  
    }
  return(0);
}

/*******************************/

int check_overlap_with_beads(double *P, double ***Chr, int n_chr, int *n_beads, int dim, double size)
{
  int i,j; /* Indeces of the 3D tensor: i is chr1, j is for the beads in chr1 */

  for(i = 0; i < n_chr; i++)
    {
      for(j = 0; j < n_beads[i]; j++)
        {
	  if( (distance(P, Chr[i][j], 3)) < size )
	    {
	      return(1);
	      //                      fprintf(stderr,"The nucleolus in position (%lf %lf %lf) and bead %d of chr %d position (%lf %lf %lf) are overlapping\n", P[0], P[1], P[2], i, j, Chr[i][j][0], Chr[i][j][1], Chr[i][j][2]);
	    }
        }
    }
  return(0);
}

/*******************************/

int check_bead_inside(double ***Chr, int n_chr, int *n_beads, int dim, double a, double b, double c)
{
  int i,j; /* Indeces of the 3D tensor: i is chr, j is for the beads in chr */


  for(i = 0; i < n_chr; i++)
    {
      //      fprintf(stderr, "---%lf %lf %lf %d %d\n", a,b,c,n_chr,n_beads[i]); 
      for(j = 0; j < n_beads[i]; j++)
	{
	  if(check_point_inside_the_confining_environment(Chr[i][j], a, b, c) == 1)
	    {
	      //	      fprintf(stderr,"Bead %d of chr %d in position: %lf %lf %lf is outside (dimensions %lf %lf %lf)\n", j, i, Chr[i][j][0], Chr[i][j][1], Chr[i][j][2], a, b, c);
	      return(1);
	    }
  	}	  
    }
  return(0);
}
