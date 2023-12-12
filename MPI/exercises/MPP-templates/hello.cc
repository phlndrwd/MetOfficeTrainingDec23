#include "mpi.h"
#include <iostream>

int main()
{
  MPI_Init(NULL, NULL);

  int rank;
  int namelen;
  char procname[MPI_MAX_PROCESSOR_NAME];
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(&procname[0], &namelen);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::cout << "Hello from " << procname << " and rank " << rank << " with " << size << " PEs " << std::endl;

  MPI_Finalize();
}
