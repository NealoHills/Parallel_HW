#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char** argv)
{
	int n = atoi(argv[1]);
	int t = atoi(argv[2]);
	char outFileName[11] = {0};
	snprintf(outFileName, 11, "%s%s", argv[1], ".txt");
	FILE* N;
	N = fopen(outFileName, "w");
	
	int* primes = (int*) malloc(sizeof(int)*(n/2));
	int i, j, num_to_cross;
	primes[0]=2;
	# pragma omp parallel for num_threads(t)
	for (i=1; i<(n/2); i++)
	{
		primes[i]=(i*2)+1;
	}
	
	int stop = (((n+1)/2)-1)/2;
	
	double tstart, ttaken;
	tstart = omp_get_wtime();
	
	
	# pragma omp parallel for num_threads(t) default(none) shared(i, stop, n, primes) private(j, num_to_cross) schedule(static, 1)
	for (i=1; i<=stop; i++)
	{
		if(primes[i] != 0)
		{
			num_to_cross = n/((i*2)+1);
			for (j=3; j<=num_to_cross; j+=2)
			{
				primes[(((((i*2)+1)*j))-1)/2] = 0;
			}
		}
	}
	
	ttaken = omp_get_wtime() - tstart;
	printf("Time taken for the main part: %f\n", ttaken);
	int num = 1;
	fprintf(N, "%d, %d, %d\n", num, 2, 0);
	int previous = 2;
	for (i=1; i<n/2; i++)
	{
		if(primes[i] != 0)
		{
			num++;
			fprintf(N, "%d, %d, %d\n", num, primes[i], primes[i]-previous);
			previous = primes[i];
		}
	}
	
	fclose(N);
	
	return 0;
}
