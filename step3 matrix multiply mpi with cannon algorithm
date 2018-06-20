//ECE 563 Project 1 Step 3 Matrix Multiply MPI version with Cannon Algorithm 2018/4/27
//Songrui Li
//0025338817
//li884

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void scatter_matrix(int* fstream, int n1, int n2, int*Q, int root, int tag)
{
	int rows = (n1 + root - 1) / root;
	int cols = (n2 + root - 1) / root;
	int* tmp_matrix = (int*)malloc(rows*cols*sizeof(int));

	int i, j;
	memset(Q, 0, rows*cols*sizeof(int));
	for (i = 0; i<root; i++)
	{
		for (j = 0; j<root; j++)
		{
			int p = 0, q = 0;
			int imin = i*rows*n2;
			int jmin = j*cols;
			memset(tmp_matrix, 0, sizeof(tmp_matrix));
			for (p = 0; p<rows; p++, imin += n2)
			{
				for (q = 0; q<cols; q++)
				{
					tmp_matrix[p*cols + q] = fstream[imin + jmin + q];
				}
			}
			if (i == 0 && j == 0)
			{
				memcpy(Q, tmp_matrix, rows*cols*sizeof(int));
			}
			else
			{
				MPI_Send(tmp_matrix, rows*cols, MPI_INT, i*root + j, tag, MPI_COMM_WORLD);
			}
		}
	}
}

int get_index(int row, int col, int sp)
{
	int tmp = ((row + sp) % sp)*sp + (col + sp) % sp;
	return tmp;
}

void matrix_multi(int* A, int *B, int *C, int n1, int n2, int n3, int myid)
{
	int i = 0, j = 0, k = 0;
	int* tmp_C = (int*)malloc(n1*n3*sizeof(int));
	memset(tmp_C, 0, sizeof(int)*n1*n3);

	for (i = 0; i<n1; i++)
	{
		for (j = 0; j<n3; j++)
		{
			for (k = 0; k<n2; k++)
			{
				tmp_C[i*n3 + j] += A[i*n2 + k] * B[k*n3 + j];
			}
			C[i*n3 + j] += tmp_C[i*n3 + j];
		}
	}
}

void shuffle(int*A, int*buf_A, int buf_A_size, int *B, int*buf_B, int buf_B_size, int root, int myid)
{
	int i, j;
	MPI_Status status;
	int cur_col = 0;
	int cur_row = 0;
	cur_row = myid / root;
	cur_col = myid - cur_row*root;
	for (i = 0; i<cur_row; i++)
	{
		MPI_Sendrecv(A, buf_A_size, MPI_INT, get_index(cur_row, cur_col - 1, root), 102,
			buf_A, buf_A_size, MPI_INT, get_index(cur_row, cur_col + 1, root), 102, MPI_COMM_WORLD, &status);
		memcpy(A, buf_A, buf_A_size*sizeof(int));
		memset(buf_A, 0, buf_A_size*sizeof(int));
	}
	for (j = 0; j<cur_col; j++)
	{
		MPI_Sendrecv(B, buf_B_size, MPI_INT, get_index(cur_row - 1, cur_col, root), 103,
			buf_B, buf_B_size, MPI_INT, get_index(cur_row + 1, cur_col, root), 103, MPI_COMM_WORLD, &status);
		memcpy(B, buf_B, buf_B_size*sizeof(int));
		memset(buf_B, 0, buf_B_size*sizeof(int));
	}
	/*printf("I have shuffled!\n");*/
}

////////////////////////Cannon Algorithm //////////////////////////////////////////////////////////
void cannon(int*A, int*buf_A, int buf_A_size, int *B, int*buf_B, int buf_B_size, int *C, int buf_C_size, int row_a, int col_a, int col_b, int root, int myid)
{	
	MPI_Status status;
	double elapsed_time, multiply_time = 0, passdata_time = 0;
	int i, j;
	memset(C, 0, sizeof(int)*buf_C_size);
	int cur_col = 0;
	int cur_row = 0;
	cur_row = myid / root;
	cur_col = myid - cur_row*root;


	for (i = 0; i<root; i++)
	{   /*root=sqrt(nprocs)*/
		elapsed_time = MPI_Wtime();
		matrix_multi(A, B, C, row_a, col_a, col_b, myid);
		elapsed_time = MPI_Wtime() - elapsed_time;
		multiply_time += elapsed_time;
		/*elapsed_time=MPI_Wtime();     */
		MPI_Sendrecv(A, buf_A_size, MPI_INT, get_index(cur_row, cur_col - 1, root), 102,
			buf_A, buf_A_size, MPI_INT, get_index(cur_row, cur_col + 1, root), 102, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(B, buf_B_size, MPI_INT, get_index(cur_row - 1, cur_col, root), 103,
			buf_B, buf_B_size, MPI_INT, get_index(cur_row + 1, cur_col, root), 103, MPI_COMM_WORLD, &status);
		/*elapsed_time=MPI_Wtime()-elapsed_time; passdata_time+=elapsed_time;*/
		memcpy(B, buf_B, buf_B_size*sizeof(int));
		memcpy(A, buf_A, buf_A_size*sizeof(int));

	}
	MPI_Send(C, row_a*col_b, MPI_INT, 0, 104, MPI_COMM_WORLD);
	/*printf("proc:%d, passdata time:%lf    multiply time:%lf\n", myid, passdata_time, multiply_time);*/
}

void gather_matrix(int *fstream, int n1, int n3, int*C, int root)
{
	MPI_Status status;
	int rows = (n1 + root - 1) / root;
	int cols = (n3 + root - 1) / root;
	int* tmp_matrix = (int*)malloc(rows*cols*sizeof(int));
	int i, j;

	for (i = 0; i<root; i++)
	{
		for (j = 0; j<root; j++)
		{
			int p, q;
			int imin = i*rows*n3;
			int jmin = j*cols;
			memset(tmp_matrix, 0, sizeof(tmp_matrix));
			MPI_Recv(tmp_matrix, rows*cols, MPI_INT, i*root + j, 104, MPI_COMM_WORLD, &status);
			/*printf("I am passed proc:%d \n",i*root+j);*/
			for (p = 0; p<rows; p++, imin += n3)
			{
				for (q = 0; q<cols; q++)
				{
					fstream[imin + jmin + q] = tmp_matrix[p*cols + q];
					/*printf("%d ",((int*)fstream)[imin+jmin+q]);*/
				}
			}
		}
	}

	for (i = 0; i<n1; i++)
	{
		for (j = 0; j<n3; j++)
		{
			printf( "%d ", fstream[i*n3 + j]);
		}
		printf("\n");
	}
}

int main(int argc, char**argv)
{
	int myid, numprocs;
	int i, j;
	MPI_Status status;
	int root = 0;
	int dim[3];
	double elapsed_time = 0;
	int max_rows_a, max_cols_a, max_rows_b, max_cols_b;
	int buf_A_size, buf_B_size, buf_C_size;
	FILE* fhc;
	/*suppose A:n1*n2 ,B:n2*n3;n1,n2,n3 are read from input file*/
	int n1, n2, n3;
	/*buffer for matrix A,B,C will be shifted ,so they each have two buffer*/
	int *A, *B, *C, *buf_A, *buf_B;

	/*on proc0,buffers to cache matrix files of A,B and C*/
	int *fstream_a = NULL, *fstream_b = NULL, *fstream_c = NULL;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	root = (int)sqrt(numprocs*1.0);
	if (numprocs != root*root)
	{
		printf("process number must be a squre!\n");
		exit(-1);
	}


	if (myid == 0)
	{
		FILE *file_a, *file_b, *file_c;
		int n1, n2, n3;
		n1 = 100;
		n2 = 120;
		n3 = 90;

	

		dim[0] = n1, dim[1] = n2, dim[2] = n3;
		fstream_a = (int*)malloc(n1*n2*sizeof(int));
		fstream_b = (int*)malloc(n2*n3*sizeof(int));

		for (i = 0; i < n1; i++)//get matrix A
		for (j = 0; j < n2; j++)
			((int*)fstream_a)[i*n2 + j] = 1;

		for (i = 0; i<n2; i++)//get matrix B
		for (j = 0; j < n3; j++)
			((int*)fstream_b)[i*n3 + j] = 2;

	}
	MPI_Bcast(dim, 3, MPI_INT, 0, MPI_COMM_WORLD);
	n1 = dim[0], n2 = dim[1], n3 = dim[2];

	max_rows_a = (n1 + root - 1) / root;
	max_cols_a = (n2 + root - 1) / root;
	max_rows_b = max_cols_a;
	max_cols_b = (n3 + root - 1) / root;
	buf_A_size = max_rows_a*max_cols_a;
	buf_B_size = max_rows_b*max_cols_b;
	buf_C_size = max_rows_a*max_cols_b;

	A = (int*)malloc(sizeof(int)*buf_A_size);
	buf_A = (int*)malloc(sizeof(int)*buf_A_size);
	B = (int*)malloc(sizeof(int)*buf_B_size);
	buf_B = (int*)malloc(sizeof(int)*buf_B_size);
	C = (int*)malloc(sizeof(int)*buf_C_size);
	if (A == NULL || buf_A == NULL || B == NULL || buf_B == NULL || C == NULL)
	{
		printf("Memory allocation failed!\n");
		exit(-1);
	}

	/*proc 0 scatter A,B to other procs in a 2D block distribution fashion*/
	if (myid == 0)
	{
		scatter_matrix((int*)fstream_a, n1, n2, A, root, 100);
		scatter_matrix((int*)fstream_b, n2, n3, B, root, 101);
	}
	else
	{
		MPI_Recv(A, max_rows_a*max_cols_a, MPI_INT, 0, 100, MPI_COMM_WORLD, &status);
		MPI_Recv(B, max_rows_b*max_cols_b, MPI_INT, 0, 101, MPI_COMM_WORLD, &status);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//////////////////////*compute C=A*B by Cannon algorithm*///////////////////////////

	shuffle(A, buf_A, buf_A_size, B, buf_B, buf_B_size, root, myid);
	elapsed_time = MPI_Wtime();
	cannon(A, buf_A, buf_A_size, B, buf_B, buf_B_size,
		C, buf_C_size, max_rows_a, max_cols_a, max_cols_b, root, myid);
	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = MPI_Wtime() - elapsed_time;

	MPI_Barrier(MPI_COMM_WORLD);

	int fsize_c = sizeof(int)*n1*n3;
	if (myid == 0)
	{
		fstream_c = (int*)malloc(fsize_c);
		gather_matrix(fstream_c, n1, n3, C, root);
	}

	MPI_Barrier(MPI_COMM_WORLD);    /*make sure proc 0 read all it needs*/

	if (myid == 0)
	{
		int i, j;
		/*printf("Cannon algorithm :multiply a %d* %d with a %d*%d, use %lf(s)\n",
			n1, n2, n2, n3, elapsed_time);
		printf("I have finished!\n");*/

		//fclose(fhc);

		free(fstream_a);
		free(fstream_b);
		free(fstream_c);
	}

	free(A); free(buf_A);
	free(B); free(buf_B);
	free(C);
	MPI_Finalize();
	return 0;
}
