/*
 * A simple solution to the Case Study exercise from the MPP course.
 * Note that this uses dynamic memory allocation for the 2D arrays.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "arralloc.h"
#include "pgmio.h"

#define MAXITER   1500
#define PRINTFREQ  200

int main (int argc, char **argv)
{
  double **old, **new, **edge;
  double **buf;

  int M, N;
  int i, j, iter, maxiter;
  char *filename;

  filename = "edge192x128.pgm";

  pgmsize(filename, &M, &N);

  printf("Processing %d x %d image\n", M, N);
  printf("Number of iterations = %d\n", MAXITER);

  /*
   * Allocate arrays dynamically using special malloc routine.
   */

  buf  = (double **) arralloc(sizeof(double), 2, M, N);

  old  = (double **) arralloc(sizeof(double), 2, M+2, N+2);
  new  = (double **) arralloc(sizeof(double), 2, M+2, N+2);
  edge = (double **) arralloc(sizeof(double), 2, M+2, N+2);

  printf("\nReading <%s>\n", filename);
  pgmread(filename, &buf[0][0], M, N);
  printf("\n");


  for (i=1;i<M+1;i++)
    {
      for (j=1;j<N+1;j++)
	{
	  edge[i][j]=buf[i-1][j-1];
	}
    }

  for (i=0; i<M+2;i++)
    {
      for (j=0;j<N+2;j++)
	{
	  old[i][j]=255.0;
	}
    }

  for (iter=1;iter<=MAXITER; iter++)
    {
      if(iter%PRINTFREQ==0)
	{
	  printf("Iteration %d\n", iter);
	}

      for (i=1;i<M+1;i++)
	{
	  for (j=1;j<N+1;j++)
	    {
	      new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]
			      - edge[i][j]);
	    }
	}
	
      for (i=1;i<M+1;i++)
	{
	  for (j=1;j<N+1;j++)
	    {
	      old[i][j]=new[i][j];
	    }
	}
    }

  printf("\nFinished %d iterations\n", iter-1);

  for (i=1;i<M+1;i++)
    {
      for (j=1;j<N+1;j++)
	{
	  buf[i-1][j-1]=old[i][j];
	}
    }

  filename="image192x128.pgm";
  printf("\nWriting <%s>\n", filename); 
  pgmwrite(filename, &buf[0][0], M, N);

  free(buf);
  free(new);
  free(old);
  free(edge);
} 
