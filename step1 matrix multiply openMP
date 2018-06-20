//ECE 563 Project 1 Step 1 Matrix Multiply OpenMP version 2018/4/1
//Songrui Li
//0025338817
//li884

#include <stdio.h>
#include <omp.h>

#define cols1 200
#define rows1 100

#define cols2 300
#define rows2 200

void main(int argc, char **argv)
{

	int i, j, k;
	if (cols1 != rows2)
	{
		printf("matrix dimensions must agree!");
		return;
	}

	int a[rows1][cols1] = { 0 };
	for (i = 0; i < rows1; i++) 
	{
		for (j = 0; j < cols1; j++) 
		{
			a[i][j] = rand();
		}
	}

	int b[rows2][cols2] = { 0 };
	for (i = 0; i < rows2; i++) 
	{
		for (j = 0; j < cols2; j++) 
		{
			b[i][j] = rand();
		}
	}

	int result[rows1][cols2] = { 0 };


	double time;
	time = -omp_get_wtime();
	omp_set_num_threads(4);


#pragma omp parallel shared(a,b,result) private(i,j,k)
	{
#pragma omp for schedule (static)
		for (i = 0; i < rows1; i++) 
		{
			for (j = 0; j < cols2; j++) 
			{
				result[i][j] = 0;
				for (k = 0; k < cols1; k++) 
				{
					result[i][j] += a[i][k] * b[k][j];
				}
			}

		}
	}
	time = omp_get_wtime() + time;

	printf("Elapsed time is %lf seconds.\n", time);

}
