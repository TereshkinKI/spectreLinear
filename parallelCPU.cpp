#include "pch.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <omp.h>

using namespace std;

int weigth(char* a, int N)
{
	int cnt = 0;
	for (int i = 0; i < N; i++)
		if (a[i] == '1')
			cnt++;
	return cnt;
}


int main(int argc, char** argv)
{
	//MPI
	int size, rank;
	MPI_Status st;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int param[2];
	int N, K;

	char* A = nullptr;

	if (rank == 0)
	{
		ifstream fin;
		fin.open("C:\\Users\\Dell\\source\\repos\\intern_research_parallel_cpu\\in.txt");
		if (fin.is_open())
		{
			fin >> N >> K;
			A = new char[K * N];
			for (int i = 0; i < K; i++)
			{
				for (int j = 0; j < N; j++)
					fin >> A[i * N + j];
			}

			fin.close();
		}
		else
			return -1337;

		param[0] = N;
		param[1] = K;
	}
	MPI_Bcast(param, 2, MPI_INT, 0, MPI_COMM_WORLD);

	N = param[0];
	K = param[1];

	if (rank)
		A = new char[K * N];

	MPI_Bcast(A, K * N, MPI_CHAR, 0, MPI_COMM_WORLD);

	int pow2k = pow(2, K);
	int* gist = new int[N + 1];
	for (int i = 0; i < N + 1; i++)
		gist[i] = 0;

	int temp_i, mult;
	
	int w;
	char* tmp;
	char* tmp_vec;
	int part = pow2k / size;

	int* gist_res = new int[N + 1];
	for (int i = 0; i < N + 1; i++)
		gist_res[i] = 0;

	double t1 = omp_get_wtime();

#pragma omp parallel for shared(gist) private(mult, temp_i, w, tmp, tmp_vec) num_threads(1) 
	for (int i = rank * part; i < (rank + 1) * part; i++)
	{
		tmp = new char[N];
		tmp_vec = new char[N];
		mult = i % 2;
		temp_i = i / 2;
		if (mult)
			for (int z = 0; z < N; z++)
			{
				tmp[z] = A[z];
			}
		else
		{
			for (int z = 0; z < N; z++)
			{
				tmp[z] = '0';
			}
		}

		for (int j = 1; j < K; j++)
		{
			mult = temp_i % 2;
			temp_i /= 2;
			if (mult)
			{
				for (int z = 0; z < N; z++)
					tmp_vec[z] = A[N*j + z];
				for (int z = 0; z < N; z++)
				{
					if (tmp[z] == tmp_vec[z])
						tmp[z] = '0';
					else
						tmp[z] = '1';
				}
			}
				
		}
		w = weigth(tmp, N);
#pragma omp atomic
		gist[w] += 1;
	}

	MPI_Reduce(gist, gist_res, N + 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	double t2 = omp_get_wtime();
	if (rank == 0)
	{
		cout << t2 - t1 << endl;
		ofstream fout;
		fout.open("C:\\Users\\Dell\\source\\repos\\intern_research_parallel_cpu\\out.txt");

		w = 0;

		for (int i = 0; i < N + 1; i++)
		{
			if (gist_res[i] != 0)
			{
				fout << i << '\t' << gist_res[i] << endl;
				w += gist_res[i];
			}
		}
		//cout << w;
		fout.close();
	}
	MPI_Finalize();
}

