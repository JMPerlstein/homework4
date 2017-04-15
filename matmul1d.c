#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>

#define USE_MPI 1

#if USE_MPI
#include <mpi.h>
#endif

static double timer() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
}

int main(int argc, char **argv) {

    int rank, num_tasks;

    /* Initialize MPI */
#if USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // printf("Hello world from rank %3d of %3d\n", rank, num_tasks);
#else
    rank = 0;
    num_tasks = 1;
#endif

    if (argc != 2) {
        if (rank == 0) {
            fprintf(stderr, "%s <n>\n", argv[0]);
            fprintf(stderr, "Program for parallel dense matrix-matrix multiplication\n");
            fprintf(stderr, "with 1D row partitioning\n");
            fprintf(stderr, "<n>: matrix dimension (an nxn dense matrix is created)\n");
#if USE_MPI
            MPI_Abort(MPI_COMM_WORLD, 1);
#else
            exit(1);
#endif
        }
    }

    int n;

    n = atoi(argv[1]);
    assert(n > 0);
    assert(n < 10000);

    /* ensure that n is a multiple of num_tasks */
    n = (n/num_tasks) * num_tasks;
    

    int n_p = (n/num_tasks);

    /* print new n to let user know n has been modified */
    if (rank == 0) {
        fprintf(stderr, "n: %d, n_p: %d\n", n, n_p);
        fprintf(stderr, "Requires %3.6lf MB of memory per task\n", ((3*4.0*n_p)*n/1e6));
    }

    float *A, *B, *C;
    
    A = (float *) malloc(n_p * n * sizeof(float));
    assert(A != 0);

    B = (float *) malloc(n_p * n * sizeof(float));
    assert(A != 0);
    
    C = (float *) malloc(n_p * n * sizeof(float));
    assert(C != 0);

    /* linearized matrices in row-major storage */
    /* A[i][j] would be A[i*n+j] */

    int i, j;

    /* static initalization, so that we can verify output */
    /* using very simple initialization right now */
    /* this isn't a good check for parallel debugging */
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
    for (i=0; i<n_p; i++) {
        for (j=0; j<n; j++) {
            A[i*n+j] = (rank+1);
            B[i*n+j] = 1;
            C[i*n+j] = 0;
        }
    }

#if USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    double elt = 0.0;
    if (rank == 0) 
        elt = timer();

#if USE_MPI
    /* Parallel matmul code goes here, see lecture slides for idea */
    /* The matrix C should be updated correctly */
    MPI_Status stat;
    
    int nextRank = rank+1;
    int prevRank = rank-1;

    int size = n* n_p;
    int offset = (rank*n_p);
    int start = n*offset;
    j = start;

    if (rank==0)
        prevRank = (num_tasks-1);
    if (rank==(num_tasks-1))
        nextRank = 0;

    B = B[start];

    //Circular Shift
    for(int shift == 0; shift < num_tasks-1; task++)
    {
        for (i=start; i<start+size; i++) 
        {
            for(j=0; j<size; j++)
            {
                C[i] = C[i] + A[i] * B[j]
            } 
        }

        int send_tag = rank;
        MPI_Sendrecv(&B[j], size, MPI_FLOAT, nextRank, MPI_ANY_TAG,
                &B, size, MPI_FLOAT, prevRank, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    }

    MPI_Barrier(MPI_COMM_WORLD);

#else
    int k;
    for (i=0; i<n_p; i++) 
    {
        for (j=0; j<n; j++) 
        {
            float c_ij = 0;
            for (k=0; k<n; k++) 
            {
                c_ij += A[i*n+k]*B[k*n+j];
            }
            C[i*n+j] = c_ij;
        }
    }
#endif

    if (rank == 0) 
        elt = timer() - elt;

    /* Verify */
    int verify_failed = 0;
    for (i=0; i<n_p; i++) {
        for (j=0; j<n; j++) {
            if (C[i*n+j] != ((rank+1)*n))
                verify_failed = 1;
        }
    }

    if (verify_failed) {
        fprintf(stderr, "ERROR: rank %d, verification failed, exiting!\n", rank);
#if USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 2);
#else
        exit(2);
#endif
    }

    if (rank == 0) {
        fprintf(stderr, "Time taken: %3.3lf s.\n", elt);
        fprintf(stderr, "Performance: %3.3lf GFlop/s\n", (2.0*n*n)*n/(elt*1e9));
    }

    /* free memory */
    free(A); free(B); free(C);

    /* Shut down MPI */
#if USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
