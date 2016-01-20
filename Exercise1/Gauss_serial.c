#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <errno.h>

/*double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}*/

double drand (double low, double high) {
    return ((double)rand() * (high - low)) / (double)RAND_MAX + low;
}

static struct timeval start, end, tstart, tend;

int main(int argc, char **argv) {
	int i, j, k, N;
	double **A, l, mtime, total;
	FILE *output;

	gettimeofday(&tstart, NULL);

	if (argc != 2) {
		perror("Wrong arguments");
		return EXIT_FAILURE;
	}

	srand(13);
	N = atoi(argv[1]);

	A = (double**)malloc(N*sizeof(double*));
	
	for ( i = 0 ; i < N ; i++) {
		A[i] = (double*)malloc(N*sizeof(double));
	}
	
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++) { 
				A[i][j] = drand(-10, 10);
		}
	}

	/*output = fopen("result.txt","w");
	
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			fprintf(output, "%lf\t" , A[i][j]);
		}
		fprintf(output, "\n");
	}*/

	gettimeofday(&start, NULL);
	
	for (k = 0 ; k < N-1 ; k++) {
		for (i = k+1; i < N ; i++) {
			l = A[i][k]/A[k][k];
			//printf("l == %lf\n",l);
			for (j = k+1;j<N;j++) {
			    A[i][j] = A[i][j] - l*A[k][j];
			}
		}
	}

	gettimeofday(&end, NULL);
	mtime = (end.tv_sec  - start.tv_sec) + (end.tv_usec - start.tv_usec)*0.000001 ;
	printf("Gauss Serial: Elapsed Calculatio time: %lf seconds\n", mtime);
	
	/*for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			fprintf(output, "%lf\t" , A[i][j]);
		}
		fprintf(output, "\n");
	}*/

	gettimeofday(&tend, NULL);
	total = (tend.tv_sec  - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec)*0.000001;
	printf("Gauss Serial: Elapsed total time: %lf seconds\n", total);
	
	return 0;
}