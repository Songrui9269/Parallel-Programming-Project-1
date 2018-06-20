//ECE 563 Project 1 Step 2 Matrix Multiply MPI version 2018/4/25
//Songrui Li
//0025338817
//li884

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define cols1 100
#define rows1 100
#define cols2 100
#define rows2 100

void main(int argc, char **argv)
{

	int taskid, numtasks;
	int i, j, k;
	int *p;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	
	int a[rows1][cols1];
	p = (int*)a;
	for (i = 0; i<rows1*cols1; i++)
		*(p + i) = i;

	int b[rows2][cols2];
	p = (int*)b;
	for (i = 0; i<rows2*cols2; i++)
		*(p + i) = i;

	int result[rows2][cols2];
	int(*temp_a)[cols1] = (int(*)[cols1])malloc(sizeof(int)* rows1 / numtasks * cols1);
	int(*temp_result)[cols2] = (int(*)[cols2])malloc(sizeof(int)* rows2 / numtasks * cols2);

	double time;
	time = -MPI_Wtime();
	MPI_Scatter(a, (rows1*cols1) / numtasks, MPI_INT, temp_a, (rows1*cols1) / numtasks, MPI_INT, 0, MPI_COMM_WORLD);

	for (i = 0; i<rows1 / numtasks; i++)
	{
		for (j = 0; j < cols2; j++)
		{
			temp_result[i][j] = 0;
			for (k = 0; k < cols1; k++)
			{
				temp_result[i][j] += temp_a[i][k] * b[k][j];
			}
		}
	}
	
	MPI_Gather(temp_result, (rows1*cols2) / numtasks, MPI_INT, result, (rows1*cols2) / numtasks, MPI_INT, 0, MPI_COMM_WORLD);
	
	time += MPI_Wtime();
	free(temp_a);
	free(temp_result);
	if (!taskid) 
		printf("Time:%lf, numtasks:%d\n", time, numtasks);
	MPI_Finalize();
}
