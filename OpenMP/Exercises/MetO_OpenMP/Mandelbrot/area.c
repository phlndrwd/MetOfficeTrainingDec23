#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h> 

# define NPOINTS 2000
# define MAXITER 2000

struct complex{
  double real;
  double imag;
};

void serial() {
  int numoutside = 0;
  double area, error;
  double start, finish;
  struct complex z, c;

  start = omp_get_wtime();
  for (int i=0; i<NPOINTS; i++) {
    for (int j=0; j<NPOINTS; j++) {
      c.real = -2.0+2.5*(double)(i)/(double)(NPOINTS)+1.0e-7;
      c.imag = 1.125*(double)(j)/(double)(NPOINTS)+1.0e-7;
      z=c;
      for (int iter = 0; iter < MAXITER; iter++) {
	      double ztemp = (z.real*z.real) - (z.imag*z.imag) + c.real;
	      z.imag = z.real * z.imag * 2 + c.imag; 
	      z.real = ztemp; 
	      if ((z.real * z.real + z.imag * z.imag) > 4.0e0) {
	        numoutside++;
	        break;
	      }
      }
    }
  }
  finish = omp_get_wtime();  
  printf("SERIAL: ");
  area = 2.0 * 2.5 * 1.125 * (double)(NPOINTS*NPOINTS-numoutside) / (double)(NPOINTS*NPOINTS);
  error = area / (double)NPOINTS;
  printf("Area of Mandlebrot set = %12.8f +/- %12.8f - ", area, error);
  printf("Time = %12.8f seconds\n", finish-start);
}

void parallel_man() {
  int numoutside = 0;
  double area, error;
  double start, finish;
  struct complex z, c;

  start = omp_get_wtime();
  int numThreads = 0;

#pragma omp parallel default(none) private(c,z) shared(numThreads) reduction(+:numoutside) 
  {
    numThreads = omp_get_num_threads();
    int thisThread = omp_get_thread_num();
    int numThreadPoints = ceil(NPOINTS / (double)numThreads);
    int numPoints = NPOINTS - (numThreads * numThreadPoints);
    int forStart = thisThread * numThreadPoints;
    int forEnd = forStart + numThreadPoints;

    if (thisThread == numThreads - 1) {
      forEnd += numPoints;
    }
    //printf("thisThread %d, numThreadPoints %d, forStart %d, forEnd %d\n", thisThread, numThreadPoints, forStart, forEnd);  
    for (int i=forStart; i<forEnd; i++) {
      for (int j=0; j<NPOINTS; j++) {
        c.real = -2.0+2.5*(double)(i)/(double)(NPOINTS)+1.0e-7;
        c.imag = 1.125*(double)(j)/(double)(NPOINTS)+1.0e-7;
  
        z=c;
        for (int iter = 0; iter < MAXITER; iter++) {
          double ztemp = (z.real*z.real) - (z.imag*z.imag) + c.real;
          z.imag = z.real * z.imag * 2 + c.imag; 
          z.real = ztemp; 
          if ((z.real * z.real + z.imag * z.imag) > 4.0e0) {
            numoutside++;
            break;
          }
        }
      }
    }
  }
  finish = omp_get_wtime();  
  printf("PARALLEL_MAN %d: ", numThreads);
  area = 2.0 * 2.5 * 1.125 * (double)(NPOINTS*NPOINTS-numoutside) / (double)(NPOINTS*NPOINTS);
  error = area / (double)NPOINTS;
  printf("Area of Mandlebrot set = %12.8f +/- %12.8f - ", area, error);
  printf("Time = %12.8f seconds\n", finish-start);
}

void parallel_for() {
  int numoutside = 0;
  double area, error;
  double start, finish;
  struct complex z, c;

  start = omp_get_wtime();
  int numThreads = 0;

#pragma omp parallel for default(none) private(c,z) shared(numThreads) reduction(+:numoutside) schedule(static, 1)
  for (int i=0; i<NPOINTS; i++) {
    numThreads = omp_get_num_threads();
    for (int j=0; j<NPOINTS; j++) {
      c.real = -2.0+2.5*(double)(i)/(double)(NPOINTS)+1.0e-7;
      c.imag = 1.125*(double)(j)/(double)(NPOINTS)+1.0e-7;

      z=c;
      for (int iter = 0; iter < MAXITER; iter++) {
	      double ztemp = (z.real*z.real) - (z.imag*z.imag) + c.real;
	      z.imag = z.real * z.imag * 2 + c.imag; 
	      z.real = ztemp; 
	      if ((z.real * z.real + z.imag * z.imag) > 4.0e0) {
	        numoutside++;
	        break;
	      }
      }
    }
  }
  finish = omp_get_wtime();  
  printf("PARALLEL_FOR %d: ", numThreads);
  area = 2.0 * 2.5 * 1.125 * (double)(NPOINTS*NPOINTS-numoutside) / (double)(NPOINTS*NPOINTS);
  error = area / (double)NPOINTS;
  printf("Area of Mandlebrot set = %12.8f +/- %12.8f - ", area,error);
  printf("Time = %12.8f seconds\n", finish-start);
}

void parallel_nest() {
  int numoutside = 0;
  double area, error;
  double start, finish;
  struct complex z, c;

  start = omp_get_wtime();
  int numThreads = 0;
  numThreads = omp_get_max_threads();

  omp_set_nested(1);
#pragma omp parallel for default(none) private(c,z) reduction(+:numoutside) schedule(static, 1)
  for (int i=0; i<NPOINTS; i++) {
#pragma omp parallel for default(none) private(c,z) shared(i) reduction(+:numoutside) schedule(static, 1)
    for (int j=0; j<NPOINTS; j++) {
      c.real = -2.0+2.5*(double)(i)/(double)(NPOINTS)+1.0e-7;
      c.imag = 1.125*(double)(j)/(double)(NPOINTS)+1.0e-7;

      z=c;
      for (int iter = 0; iter < MAXITER; iter++) {
	      double ztemp = (z.real*z.real) - (z.imag*z.imag) + c.real;
	      z.imag = z.real * z.imag * 2 + c.imag; 
	      z.real = ztemp; 
	      if ((z.real * z.real + z.imag * z.imag) > 4.0e0) {
	        numoutside++;
	        break;
	      }
      }
    }
  }
  finish = omp_get_wtime();  
  printf("PARALLEL_NEST %d: ", numThreads);
  area = 2.0 * 2.5 * 1.125 * (double)(NPOINTS*NPOINTS-numoutside) / (double)(NPOINTS*NPOINTS);
  error = area / (double)NPOINTS;
  printf("Area of Mandlebrot set = %12.8f +/- %12.8f - ", area, error);
  printf("Time = %12.8f seconds\n", finish-start);
}

int main() {
  //serial();
  parallel_man();
  parallel_for();
  parallel_nest();
}
