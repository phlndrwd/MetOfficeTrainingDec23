//
// Prints 'Hello World' from rank 0 and 
// prints what processor it is out of the total number of processors from all ranks
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(){
  int rank, size, ierr;
  MPI_Comm comm;

  comm  = MPI_COMM_WORLD;

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(comm, &rank);            
  MPI_Comm_size(comm, &size); 

  // Each process prints out its rank
  printf("I am %d out of %d\n", rank,size);

  // Only processor 0 prints 
  if(rank == 0) {
   printf("Hello World!\n"); 
  }

  MPI_Finalize();

}
