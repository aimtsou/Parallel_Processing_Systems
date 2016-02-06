#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <errno.h>
#include <mpi.h>

double drand (double low, double high) {
    return ((double)rand() * (high - low)) / (double)RAND_MAX + low;
}

static struct timeval start, end, tstart, tend;

int main(int argc, char **argv) {
	int i, j, k, N, M, n, tmp, tmp2, tmp3, rank, size;
	double *A, *recvbuf, l, *pline, mtime, total;
	FILE *output;

	gettimeofday(&tstart, NULL);
	
	if (argc != 2) {
		perror("Wrong arguments");
		return EXIT_FAILURE;
	}

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	srand(13);
	N = atoi(argv[1]);

	if (N%size != 0)
		M = N + size - (N%size);
	else 
		M = N;
		
	if (rank == 0) {
		// allocate matrix and randomize
		A = (double*) malloc (M*M*sizeof(double));
		for (i = 0; i < M; i++){
			for (j = 0; j < M; j++) {
				if ((i < N) && (j < N)) 
					A[M*i + j] = drand(-10, 10);
				else
					A[M*i + j] = 0;
			}
		}
	}
	
	// allocate recvbuf and scatter
	tmp = M/size;
	recvbuf = (double*) malloc (tmp*M*sizeof(double));
	MPI_Status stat;
	
	for (j=0; j<M; j++) {
		if (rank==0) {
			MPI_Send(&A[M*j], M, MPI_DOUBLE, j%size, 50, MPI_COMM_WORLD);
		}
		if (rank==(j%size)) {
			MPI_Recv(&recvbuf[M*(j/size)], M, MPI_DOUBLE, 0, 50, MPI_COMM_WORLD, &stat);
		}				
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	// start timer
	if (rank == 0) 
		gettimeofday(&start, NULL);
		
	for (k = 0; k < M-1; k++) {
		tmp2 = k % size;
		int dest;
		pline = (double*) malloc (M*sizeof(double));
		if (tmp2 == rank) {
			tmp3 = k/size;
			for (n = 0; n < M-k; n++) {
				pline[n] = recvbuf[tmp3*M + n + k];
			}
			for (dest=0; dest<size; dest++) {	
				if (dest!=rank) {	
					MPI_Send (pline,M-k,MPI_DOUBLE,dest,50,MPI_COMM_WORLD);
				}
			} 
		}
		else {	
			MPI_Recv(pline,M-k,MPI_DOUBLE,tmp2,50,MPI_COMM_WORLD,&stat);
		}

		for (i = k+1; i < N; i++) {
			tmp3 = i/size;
			if ((i % size) == rank) {
				if (pline[0] != 0) {
					l = recvbuf[tmp3*M + k]/pline[0];
					for (j = k+1; j < N; j++) {
						recvbuf[tmp3*M + j] -= l*pline[j-k];
					}
				}
			}
		}
		free(pline);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {
		gettimeofday(&end, NULL);
		mtime = (end.tv_sec  - start.tv_sec) + (end.tv_usec - start.tv_usec)*0.000001 ;
		printf("Gauss 2: Elapsed Calculation time: %lf seconds\n", mtime);
	}
	
	for (j=0; j<M; j++) {	
		if (rank==(j%size)) {
			MPI_Send (&recvbuf[M*(j/size)], M, MPI_DOUBLE, 0, 50, MPI_COMM_WORLD);
		}	
		if (rank==0) {
			MPI_Recv (&A[M*j], M, MPI_DOUBLE, j%size, 50, MPI_COMM_WORLD, &stat);
		}			
	}
	
	//output = fopen("result.txt","w");

	if (rank == 0) {
		// print matrix
		/*for (i = 0; i < N; i++){
			for (j = 0; j < N; j++) {
				fprintf(output, "lf", A[N*i + j]);
			}
			fprintf(output, "\n");
		}*/
		// free memory
		free(recvbuf);
		free(A);
	}

	MPI_Finalize();

	gettimeofday(&tend, NULL);
	total = (tend.tv_sec  - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec)*0.000001;
	printf("Gauss 4: Elapsed total time: %lf seconds\n", total);

	return 0;
}
