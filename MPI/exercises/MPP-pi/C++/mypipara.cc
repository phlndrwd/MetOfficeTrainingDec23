#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

#include "mpi.h"

#define N 84000000

int main() {
  std::chrono::time_point<std::chrono::system_clock> startTime = std::chrono::system_clock::now();
  std::chrono::time_point<std::chrono::system_clock> endTime;

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Status status;

  int rank;
  int size;
  int source;
  int tag;

  long double pi;
  long double recvpi;

  constexpr auto max_precision{std::numeric_limits<long double>::digits10 + 1}; 

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  if (rank == 0) {
    std::cout << "Computing approximation to pi using N = " << N << std::endl;
  }
  
  int istart = N / size * rank + 1;
  int istop  = istart + N / size - 1;

  long double partialpi = 0.0;
  for (int i = istart; i <= istop; i++) {
    partialpi = partialpi + 1.0/( 1.0 + pow((((long double)i) - 0.5) / ((long double)N), 2.0));
  }

  if (rank == 0) {
    pi = partialpi;
    for (source = 1; source < size; source++) {
      tag = 0;
      MPI_Recv(&recvpi, 1, MPI_LONG_DOUBLE, MPI_ANY_SOURCE, tag, comm, &status);
      pi = pi + recvpi;
    }
  } else {
    tag = 0;
    MPI_Ssend(&partialpi, 1, MPI_LONG_DOUBLE, 0, tag, comm);
  }

  pi = pi * 4.0 / ((long double) N);

  long double exactpi = 4.0 * atan(1.0);
  long double error = fabs(100.0 * (pi- exactpi) / exactpi);
  
  if (rank == 0) {
    std::cout << "pi = " << std::setprecision(max_precision) << pi << ", error = " << error << "%" << std::endl;
  }

  MPI_Finalize();
  if (rank == 0) {
    std::cout << "Time taken = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime).count() << "ms" << std::endl;
  }
  return 0;
}
