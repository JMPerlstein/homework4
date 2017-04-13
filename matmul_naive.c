#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>

static double timer() {
    
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);

}

int main(int argc, char **argv) {

    if (argc != 2) {
        fprintf(stderr, "%s <n>\n", argv[0]);
        fprintf(stderr, "<n>: matrix dimension (nxn dense matrices are created)\n");
        exit(1);
    }

    int n;

    n = atoi(argv[1]);
    assert(n > 0);
    assert(n < 10000);

    fprintf(stderr, "n: %d\n", n);
    fprintf(stderr, "Requires %3.6lf MB memory\n", ((3*8.0*n)*n/1e6));

    double *A, *B, *C;
    
    A = (double *) malloc(n * n * sizeof(double));
    assert(A != 0);
    B = (double *) malloc(n * n * sizeof(double));
    assert(B != 0);
    C = (double *) malloc(n * n * sizeof(double));
    assert(C != 0);


    /* linearized matrices in row-major storage */
    /* A[i][j] would be A[i*n+j] */

    int i, j;

    /* static initalization, so that we can verify output */
    /* using very simple initialization right now */
    /* this isn't a good check for parallel debugging */
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            A[i*n+j] = 1;
            B[i*n+j] = 1;
            C[i*n+j] = 0;
        }
    }

    double elt = 0.0;
    elt = timer();

    int k, ii;
    for (ii=0; ii<n; ii+= 4) 
    {
        for (j=0; j<n; j+= 2) 
        {
            for (i = ii; i < ii + 4; i+= 2)
            {    
                int in = i*n;
                int iin = (i+1)*n;
                double c_ij[4] = {0};
                for (k=0; k<n; k+=2) 
                {
                    c_ij[0] += A[in+k]*B[k*n+j];
                    c_ij[1] += A[in+k]*B[k*n+(j+1)];
                    c_ij[2] += A[iin+k]*B[k*n+j];
                    c_ij[3] += A[iin+k]*B[k*n+(j+1)];

                    c_ij[0] += A[in+(k+1)]*B[(k+1)*n+j];
                    c_ij[1] += A[in+(k+1)]*B[(k+1)*n+(j+1)];
                    c_ij[2] += A[iin+(k+1)]*B[(k+1)*n+j];
                    c_ij[3] += A[iin+(k+1)]*B[(k+1)*n+(j+1)];
                }
                C[in+j] = c_ij[0];
                C[in+(j+1)] = c_ij[1];
                C[iin+j] = c_ij[2];
                C[iin+(j+1)] = c_ij[3];
            }
        }
    }

    elt = timer() - elt;

    /* Verify */
    int verify_failed = 0;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            if (C[i*n+j] != n)
                verify_failed = 1;
        }
    }

    if (verify_failed) {
        fprintf(stderr, "ERROR: verification failed, exiting!\n");
        exit(2);
    }

    fprintf(stderr, "Time taken: %3.3lf s.\n", elt);
    fprintf(stderr, "Performance: %3.3lf GFlop/s\n", (2.0*n*n*n)/(elt*1e9));

    /* free memory */
    free(A); free(B); free(C);

    return 0;
}
