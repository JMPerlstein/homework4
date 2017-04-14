#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

int main(int argc, char **argv) {

    int rank, nprocs;
    int i;

    /* Initialize MPI Environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Hello world */
    /*
    printf("Hello world from rank %3d of %3d\n", rank, nprocs);
    */
    
    /* Exchange some messages between rank 0 and rank n-1, 
       while other processes wait */

    MPI_Status stat;
    if (rank == 0) {

        int sendval = 42;
        int message1_tag = 1;
        int send_retval = 1;
        send_retval = MPI_Send(&sendval, 1, MPI_INT, nprocs-1, message1_tag, MPI_COMM_WORLD);
        /* MPI_Send will return 0 on successful send */
        assert(send_retval == 0);

        char recvstr[100];
        int message2_tag = 2;
        MPI_Status stat;
        MPI_Recv(&recvstr, 100, MPI_CHAR, nprocs-1, message2_tag, MPI_COMM_WORLD, &stat);
        fprintf(stderr, "Task 0 received string %s from task %d. " 
                "Message tag is %d\n", recvstr, stat.MPI_SOURCE, stat.MPI_TAG);
    }

    if (rank == (nprocs-1)) {

        int recvval = 0;
        /* setting tag to 1 */
        MPI_Recv(&recvval, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);
        fprintf(stderr, "Task %d received value %d from task %d.\n",
                nprocs-1, recvval, stat.MPI_SOURCE, stat.MPI_TAG);
 
        char sendstr[100] = "Hello world";
        int message2_tag = 2;
        MPI_Send(&sendstr, 100, MPI_CHAR, 0, message2_tag, MPI_COMM_WORLD);
    }

    /* A barrier is needed so that other processes wait */
    MPI_Barrier(MPI_COMM_WORLD);



    /* The above example was quite contrived. We usually have all processes
       sending and receiving data */
    /* Here's an example of every process sending some data to the left and 
       receiving from the right */
    int sendval = rank;
    int recvval = -1;

    int my_left_rank  = rank - 1;
    int my_right_rank = rank + 1;

    if (rank == 0)
        my_left_rank = nprocs-1;

    if (rank == (nprocs-1))
        my_right_rank = 0;

    int send_tag = rank;

    /* receiver tag can be MPI_ANY_TAG */
    MPI_Sendrecv(&sendval, 1, MPI_INT, my_left_rank, send_tag,
           &recvval, 1, MPI_INT, my_right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &stat); 

    /* We know the value to be received is my_right_rank */
    assert(recvval == my_right_rank);    
    MPI_Barrier(MPI_COMM_WORLD);


    /* Gather a bunch of integer values to one process */
    sendval = rank;
    int *recvbuf;
   
    /* Let's gather data to process 0 */
    int root = 0; 
    if (rank == root) {
        recvbuf = (int *) malloc(nprocs * sizeof(int));
        assert(recvbuf != NULL);
    }

    /* There are implicit barriers before/after collectives */
    MPI_Gather(&sendval, 1, MPI_INT, recvbuf, 1, MPI_INT, root, MPI_COMM_WORLD);
    
    if (rank == root) {
        /* verify received data */
        for (i=0; i<nprocs; i++) {
            assert(recvbuf[i] == i);
        }
        free(recvbuf);
    }    

    MPI_Barrier(MPI_COMM_WORLD);

    /* Try Allreduce */
    double *A;
    int A_size = 4096;

    double *B;

    A = (double *) malloc(A_size * sizeof(double));
    assert(A != NULL); 

    B = (double *) malloc(A_size * sizeof(double));
    assert(B != NULL); 

    for (i=0; i<A_size; i++) {
        A[i] = ((double) rank);
    }

    MPI_Allreduce(A, B, A_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    /* Verify */
    for (i=0; i<A_size; i++) {
        assert(B[i] == ((nprocs*(nprocs-1))/2.0));
    }

    free(A); free(B);



    /* Terminate MPI environment */
    MPI_Finalize();

    return 0;
}
