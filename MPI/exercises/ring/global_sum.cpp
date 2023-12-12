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
  MPI_Request requestSend;
  MPI_Request requestRec;
  MPI_Status status;

  int rank;
  int numPEs;
  int messageIn;
  int messageOut;

  int tag = 0;
  int sum = 0;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &numPEs);

  int dest = rank + 1;
  if (dest >= numPEs) {
    dest = 0;
  }
  int source = rank - 1;
  if (source < 0) {
    source = numPEs - 1;
  }
  
  messageOut = rank;
  for(int i = 0; i < numPEs; ++i) {
    MPI_Issend(&messageOut, 1, MPI_INT, dest, tag, comm, &requestSend);
    MPI_Irecv(&messageIn, 1, MPI_I
    NT, source, tag, comm, &requestRec);
    MPI_Wait(&requestSend, &status);
    MPI_Wait(&requestRec, &status);

    sum += messageIn;
    messageOut = messageIn;
  }

  MPI_Finalize();
  if (rank == 0) {
    std::cout << "sum = " << sum << std::endl;
    std::cout << "Time taken = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime).count() << "ms" << std::endl;
  }
  return 0;
}
