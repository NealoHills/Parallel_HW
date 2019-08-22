#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/*** Skeleton for Lab 1 ***/

/***** Globals ******/
float **a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float err; /* The absolute relative error */
int num = 0;  /* number of unknowns */
int comm_sz;
int my_rank;


/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */

/********************************/



/* Function definitions: functions are ordered alphabetically ****/
/*****************************************************************/

/* 
   Conditions for convergence (diagonal dominance):
   1. diagonal element >= sum of all other elements of the row
   2. At least one diagonal element > sum of all other elements of the row
 */
void check_matrix()
{
  int bigger = 0; /* Set to 1 if at least one diag element > sum  */
  int i, j;
  float sum = 0;
  float aii = 0;
  
  for(i = 0; i < num; i++)
  {
    sum = 0;
    aii = fabs(a[i][i]);
    
    for(j = 0; j < num; j++)
       if( j != i)
	 sum += fabs(a[i][j]);
       
    if( aii < sum)
    {
      printf("The matrix will not converge.\n");
      exit(1);
    }
    
    if(aii > sum)
      bigger++;
    
  }
  
  if( !bigger )
  {
     printf("The matrix will not converge\n");
     exit(1);
  }
}


/******************************************************/
/* Read input from file */
/* After this function returns:
 * a[][] will be filled with coefficients and you can access them using a[i][j] for element (i,j)
 * x[] will contain the initial values of x
 * b[] will contain the constants (i.e. the right-hand-side of the equations
 * num will have number of variables
 * err will have the absolute error that you need to reach
 */
void get_input(char filename[])
{
  FILE * fp;
  int i,j;  
 
  fp = fopen(filename, "r");
  if(!fp)
  {
    printf("Cannot open file %s\n", filename);
    exit(1);
  }

 fscanf(fp,"%d ",&num);
 fscanf(fp,"%f ",&err);

 /* Now, time to allocate the matrices and vectors */
 a = (float**)malloc(num * sizeof(float*));
 if( !a)
  {
	printf("Cannot allocate a!\n");
	exit(1);
  }

 for(i = 0; i < num; i++) 
  {
    a[i] = (float *)malloc(num * sizeof(float)); 
    if( !a[i])
  	{
		printf("Cannot allocate a[%d]!\n",i);
		exit(1);
  	}
  }
 
 x = (float *) malloc(num * sizeof(float));
 if( !x)
  {
	printf("Cannot allocate x!\n");
	exit(1);
  }


 b = (float *) malloc(num * sizeof(float));
 if( !b)
  {
	printf("Cannot allocate b!\n");
	exit(1);
  }

 /* Now .. Filling the blanks */ 

 /* The initial values of Xs */
 for(i = 0; i < num; i++)
	fscanf(fp,"%f ", &x[i]);
 
 for(i = 0; i < num; i++)
 {
   for(j = 0; j < num; j++)
     fscanf(fp,"%f ",&a[i][j]);
   
   /* reading the b element */
   fscanf(fp,"%f ",&b[i]);
 }
 
 fclose(fp); 

}


/************************************************************/

int shittyRREF(int nit) {
  /*
  printf("num in process %d: %d\n", my_rank, num);
  int i;
  
   for(i=0; i<num; i++)
 {
	 printf("%f\n",x[i]);
 }
 */

  // determine number of elements each process will work with (need to add its share of remainder)
  
  int *counts = (int*)malloc(num * sizeof(int));
  int *displs = (int*)malloc(num * sizeof(int));
  
  int remainder = num % comm_sz;
  int distr = num / comm_sz; 
  int start;
  if (my_rank < remainder) {
  	start = my_rank * distr + my_rank;
  } else {
  	start = my_rank * distr + remainder;
  }
  int finish = start + distr-1;
  if (my_rank < remainder) finish += 1;
  
  int loc_count = finish - start + 1;

  MPI_Allgather(&loc_count, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&start, 1, MPI_INT, displs, 1, MPI_INT, MPI_COMM_WORLD);
  
  //printf("%d\n", loc_count);
  
  float* new_x = (float*)malloc(loc_count * sizeof(float)); // keep track of newly calculated x values; use this in call to Allgather to create new x

  // a and b don't change so we can use their global values and not have to distribute them
  // each process has an array x

  // put new x's in a new array
  // do error calculation in each process
  // if the percentage is over, set a flag that says to loop again
  // use reduce to find the max of those flags (1 means loop again)
  // Allgather x-news and put them into x-olds
  // Allgather automatically synchronizes everything!
  int i;
  int j;
  float sum_x;
  float maxerr = 0.0;
  int repeat = 1;
  while (repeat) { // change this to accept a flag
    maxerr = 0.0;
  	if (my_rank == 0) nit++;
    
    // calculate new x values based on old ones
    for (i = 0; i < loc_count; i++) {
      // subtract all x's (except the x at j) * their corresponding a[j] and divide by a[i]
      sum_x = 0;
      for (j = 0; j < num; j++) {
      	if (j != i+start) {
      	  sum_x -= x[j] * a[i+start][j]; // follow calculations from step 1 in instructions
      	}
      }
      new_x[i] = (b[i+start] + sum_x) / a[i+start][i+start];
	  //printf("%d: %f\n", my_rank, new_x[i]);
    }

    // calculate percent error for each x in local_x, keeping track of the maximum
    // use this maximum to set the flag for whether or not to continue the while loop
	float cur_err;
    for (i = 0; i < loc_count; i++) {
        cur_err = fabs((new_x[i] - x[i + start]) / new_x[i]);
      if (cur_err > maxerr) { 
      	maxerr = cur_err;
      }
    }

	//printf("error: %f\n", maxerr);
	
    int flag = 0;
    if (maxerr >= err) flag = 1; // if the highest percent error is greater than err, we need to keep looping

    // set error to the result of reducing the flags, based on whether or not maxerr <= err
    MPI_Allreduce(&flag, &repeat, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); // use flag from each process and * distribute result to each process *
    // if any of the processes returns 1 for the flag, repeat will be set to 1 and the loop will continue,
    // otherwise, repeat will be set to 0 and the loop will not continue

    //printf("%d\n", loc_count);
    MPI_Allgatherv(new_x, loc_count, MPI_FLOAT, x, counts, displs, MPI_FLOAT, MPI_COMM_WORLD); // concatenate all our new x values (from new_x[]) and put them into x[] (so they will be treated as the initial values of our next iteration)
  }
  
  free(new_x); // values from new_x have been put into x so it is safe to free it

  // x should now hold all values within the proper percent error
  // output these values to stdout (in main)
  return nit;

}


/************************************************************/

int main(int argc, char *argv[])
{
    
    int i;
    int nit = 0; /* number of iterations */
    
    
    
     
    if( argc != 2)
    {
      printf("Usage: gsref filename\n");
      exit(1);
    }
    
       
    /* Read the input file and fill the global data structure above */ 
    //if(my_rank == 0){
    get_input(argv[1]);
    
    /* Check for convergence condition */
    /* This function will exit the program if the coffeicient will never converge to 
     * the needed absolute error. 
     * This is not expected to happen for this programming assignment.
     */
    	check_matrix();
    //}
    /*
    MPI_Bcast(a, sizeof(a), MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, sizeof(b), MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(x, sizeof(x), MPI_FLOAT, 0, MPI_COMM_WORLD);
    */
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    double local_start, local_finish, local_elapsed, elapsed;
    
    MPI_Barrier(MPI_COMM_WORLD);
    local_start = MPI_Wtime();
    //printf("num in main:%d\n", num);
    if(num < comm_sz && my_rank==0)
    {
    	printf("There are only %d unknowns, you don't need %d processes to solve this.\n", num, comm_sz); 
    }
    else
    {
        nit = shittyRREF(nit);
        
        local_finish = MPI_Wtime();
        local_elapsed = local_finish-local_start;
        MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if(my_rank==0){
            /* Writing to the stdout */
            /* Keep that same format */
            for( i = 0; i < num; i++)
                printf("%f\n",x[i]);
            
            printf("Total number of iterations: %d\n", nit);
            printf("Elapsed time = %e seconds\n", elapsed);
        }
    }
    
    MPI_Finalize();
    
    exit(0);
    
}