#include <cstdlib>
#include <cmath>

#include "mpi.h"

#include <iostream>

#define N 840

int main(void) {
  double pi, exactpi;

  // MPI variables
  MPI_Comm comm;
  MPI_Status status;

  int rank, size, source, tag;

  // Other variables
  int istart, istop;
  double partialpi, recvpi;

  std::cout << "---- Computing approximation to pi using N = " << N << std::endl;

  // Initialise MPI and compute number of processes and local rank
  comm = MPI_COMM_WORLD;

  MPI_Init(NULL, NULL);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  std::cout << "---- Hello from rank " << rank << std::endl;

  // Now make sure output only comes from one process
  if (rank == 0) {
    std::cout << "---- Running on " << size << " process(es)" << std::endl;
  }
  partialpi = 0.0;

  //
  // Compute an approximation to pi using a simple series expansion for pi/4
  // Ensure each process computes a separate section of the summation
  // NOTE: here I assume that N is exactly divisible by the number of processes
  //

  istart = N / size * rank + 1;
  istop  = istart + N / size - 1;

  std::cout << "---- On rank " << rank << " istart = " << istart << ", istop = " << istop << std::endl;
  for (int i = istart; i <= istop; i++) {
    partialpi = partialpi + 1.0 / (1.0 + pow((((double) i) - 0.5) / ((double) N), 2.0));
  }
  std::cout << "---- On rank " << rank << " partialpi = " << partialpi << std::endl;

  //
  // Compute global value of pi by sending partial values to rank 0 
  // NOTE: this would be more efficiently done using MPI_REDUCE 
  //

  if (rank == 0) {
    // Initialise pi to locally computed parial sum
    pi = partialpi;
    // Add in contribution from other processes
    for (source = 1; source < size; source++) {
      // receive partialpi from rank=source and place value in recvpi
      // all messages are tagged as zero
      tag = 0;
  
      MPI_Recv(&recvpi, 1, MPI_DOUBLE, MPI_ANY_SOURCE, tag, comm, &status);
  
      std::cout << "---- rank 0 receiving from rank " << status.MPI_SOURCE << std::endl;
  
      // add to running total
      pi = pi + recvpi;
    }
  } else {
    // all other processes send their partial value to rank 0
    tag = 0;
    std::cout << "---- rank " << rank << " sending to rank 0" << std::endl;
    MPI_Ssend(&partialpi, 1, MPI_DOUBLE, 0, tag, comm);
  }

  pi = pi * 4.0/((double) N);
  exactpi = 4.0 * atan(1.0);
  
  if (rank == 0) {
    std::cout << "pi = " << pi << ", % error = " << fabs(100.0 * (pi-exactpi) / exactpi) << std::endl;
  }

  MPI_Finalize();
  return 0;
}
