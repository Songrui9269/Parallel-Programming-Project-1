//ECE 563 Project1 Step 1 Binary Tree Search OpenMP version 2018/4/25
//Songrui Li
//0025338817
//li884

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define ID 0
#define AGE 1 
#define LEFT 2
#define RIGHT 3
#define rows 100

typedef struct
{
	int id;
	int age;
}BTree;

void locate(int index, int* node_array, BTree * node)
{
	node_array[index * 4 + ID] = node->id;
	node_array[index * 4 + AGE] = node->age;
	node_array[index * 4 + LEFT] = -1;
	node_array[index * 4 + RIGHT] = -1;
}

void InsertNode(int root, BTree * node, int index, int * node_array)
{
	if (node->age > node_array[root * 4 + AGE])
	{
		if (node_array[root * 4 + RIGHT] == -1)
		{
			node_array[root * 4 + RIGHT] = index;
			locate(index, node_array, node);
		}
		else
		{
			InsertNode(node_array[root * 4 + RIGHT], node, index, node_array);
		}
	}
	else
	{
		if (node_array[root * 4 + LEFT] == -1)
		{
			node_array[root * 4 + LEFT] = index;
			locate(index, node_array, node);
		}
		else
		{
			InsertNode(node_array[root * 4 + LEFT], node, index, node_array);
		}
	}
}


void  CreateTree(int size, int * node_array)
{
	int i;
	BTree * tree = (BTree *)malloc(sizeof(BTree));
	tree->id = rand() % 10000;
	tree->age = rand() % 40 + 10;
	locate(0, node_array, tree);
	for (i = 1; i < size; i++)
	{
		BTree * node = (BTree *)malloc(sizeof(BTree));
		node->id = rand() % size;
		node->age = rand() % 60 + 10;
		InsertNode(0, node, i, node_array);
	}

}

void search(int age, int * node_array, int index)
{
	if (node_array[index * 4 + AGE] == age)
	{
		printf("ID = %d, Age= %d\n", node_array[index * 4 + ID], node_array[index * 4 + AGE]);
	}
	if (node_array[index * 4 + LEFT] != -1)
	{
		//#pragma omp task
		search(age, node_array, node_array[index * 4 + LEFT]);
	}
	if (node_array[index * 4 + RIGHT] != -1)
	{
		//#pragma omp task
		search(age, node_array, node_array[index * 4 + RIGHT]);
	}

}

void check(int age, int * node_array, int index)
{
	if (node_array[index * 4 + AGE] == age)
	{
		printf("ID = %d, Age= %d\n", node_array[index * 4 + ID], node_array[index * 4 + AGE]);
	}
}

void get_more_work(int * work_list, int * num, int numP, int age, int * node_array)
{
	if ((*num) == numP)
	{
		return;
	}
	else
	{
		(*num)--;
		int root = work_list[*num];
		check(age, node_array, root);
		int left = node_array[root * 4 + LEFT];
		int right = node_array[root * 4 + LEFT];
		if (left != -1)
		{
			work_list[*num] = left;
			(*num)++;
		}
		if (right != -1)
		{
			work_list[*num] = right;
			(*num)++;
		}
		get_more_work(work_list, num, numP, age, node_array);
	}
}


void main(int argc, char **argv)
{
	int tree[rows][4]; // 0 :ID, 1 :age, 2 :left , 3 : right
	CreateTree(rows, (int *)tree);

	int tragetAGE = 30;

	double time;
	time = -omp_get_wtime();
	
	int* work_list = (int*)malloc(sizeof(int));
	int work = 1;
	work_list[0] = 0;
	get_more_work(work_list, &work, 1, tragetAGE, (int *)tree);


	int begin_node = 0;
	search(50, (int *)tree, begin_node);
	time = omp_get_wtime() + time;

	printf("Elapsed time is %lf seconds.\n", time);

}

