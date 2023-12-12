#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

#include "mpi.h"

int main() {
  std::chrono::time_point<std::chrono::system_clock> startTime = std::chrono::system_clock::now();
  std::chrono::time_point<std::chrono::system_clock> endTime;

  MPI_Init(NULL, NULL);
  MPI_Comm comm = MPI_COMM_WORLD;

  int rank;
  int sum = 0;

  MPI_Comm_rank(comm, &rank);
  MPI_Allreduce(&rank, &sum, 1, MPI_INT, MPI_SUM, comm);
  MPI_Finalize();

  if (rank == 0) {
    std::cout << "sum = " << sum << std::endl;
    std::cout << "Time taken = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime).count() << "ms" << std::endl;
  }
  return 0;
}
