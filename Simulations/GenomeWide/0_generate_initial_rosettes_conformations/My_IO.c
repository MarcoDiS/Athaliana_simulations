#include<stdio.h>
#include<stdlib.h>

FILE *open_r(char *filename)
{
  FILE *fp;

  if( (fp = fopen(filename, "r")) == NULL )
    {
      fprintf(stderr, "The file %s does not exist \n", filename);
      exit(0);    
    } 

  return(fp);
}
/*************************/

FILE *open_w(char *filename)
{
  FILE *fp;

  if( (fp = fopen(filename, "w")) == NULL )
    {
      fprintf(stderr, "It's impossible to open the file %s \n", filename);
      exit(0);    
    } 

  return(fp);
}
/*************************/

FILE *open_a(char *filename)
{
  FILE *fp;

  if( (fp = fopen(filename, "a")) == NULL )
    {
      fprintf(stderr, "It's impossible to open the file %s \n", filename);
      exit(0);    
    } 

  return(fp);
}
/*************************/
