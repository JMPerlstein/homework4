#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

int my_Allgather(int *sendbuf, int n, int nprocs, int *recvbuf) {

    for (int i=0; i<n*nprocs; i++) {
        recvbuf[i] = 0;
    }

    /* recursive doubling-based code goes here */

    return 0;
}

int main(int argc, char **argv) {

    int rank, nprocs;
    int i;

    /* Initialize MPI Environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* change the input size here */
    int n;
    if (argc == 2) {
        n = atoi(argv[1]);
    } else {
        n = 4;
    }
    if (rank == 0) {
        fprintf(stderr, "n: %d\n", n);
    }

    int *sendbuf;
    int *recvbuf1;
    int *recvbuf2;
    
    sendbuf  = (int *) malloc(n*sizeof(int));
    assert(sendbuf != 0);

    for (i=0; i<n; i++) {
        sendbuf[i] = (rank+1);
    }

    recvbuf1 = (int *) malloc(n*nprocs*sizeof(int));
    assert(recvbuf1 != 0);

    recvbuf2 = (int *) malloc(n*nprocs*sizeof(int));
    assert(recvbuf2 != 0);

    MPI_Allgather(sendbuf, n, MPI_INT, recvbuf1, n, MPI_INT, MPI_COMM_WORLD);
    my_Allgather(sendbuf, n, nprocs, recvbuf2);
    
    /* verify that my_Allgather works correctly */
    for (i=0; i<n*nprocs; i++) {
        assert(recvbuf1[i] == recvbuf2[i]);
    }

    free(sendbuf); free(recvbuf1); free(recvbuf2);

    /* Terminate MPI environment */
    MPI_Finalize();

    return 0;
}
