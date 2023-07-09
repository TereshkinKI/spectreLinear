#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <omp.h>

using namespace std;


__device__ char* XOR(char* a, char* b, int N)
{
	char* c = new char[N];
	for (int i = 0; i < N; i++)
	{
		if (a[i] != b[i])
			c[i] = '1';
		else
			c[i] = '0';
	}

	return c;
}


__device__ int weigth(char* a, int N)
{
	int cnt = 0;
	for (int i = 0; i < N; i++)
		if (a[i] == '1')
			cnt++;
	return cnt;
}



__global__ void kernel(char* A, int N, int K, int pow2k, int* gist, char* tmp, char* tmp_vec)
{
	int id = blockDim.x * blockIdx.x + threadIdx.x;
	int mult, temp_i, w;
	
	if (id < pow2k)
	{
		mult = id % 2;
		temp_i = id / 2;
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
				tmp = XOR(tmp, tmp_vec, N);
			}
		}
		gist[id] = weigth(tmp, N);
	}
}


int main()
{
	ifstream fin;
	fin.open("in.txt");

	int N, K;
	fin >> N >> K;

	char* A = new char[N * K];

	for (int i = 0; i < N * K; i++)
		fin >> A[i];
	fin.close();

	int pow2k = pow(2, K);

	int threads = 32;
	int blocks = pow2k / threads + 1;

	int* gist = new int[pow2k];

	for (int i = 0; i < pow2k; i++)
		gist[i] = 0;

	double t1 = omp_get_wtime();

	char* A_d = new char [N * K];
	cudaMalloc((void**)&A_d, N * K * sizeof(char));

	int* gist_d = new int  [pow2k];
	cudaMalloc((void**)&gist_d, (pow2k) * sizeof(int));

	char* tmp = new char[N];
	char* tmp_vec = new char[N];

	char* tmp_d;
	cudaMalloc((void**)&tmp_d, N * sizeof(char));
	char* tmp_vec_d;
	cudaMalloc((void**)&tmp_vec_d, N * sizeof(char));

	cudaMemcpy(A_d, A, N * K * sizeof(char), cudaMemcpyHostToDevice);
	cudaMemcpy(gist_d, gist, (pow2k) * sizeof(int), cudaMemcpyHostToDevice);

	kernel <<< blocks, threads >>> (A_d, N, K, pow2k, gist_d, tmp_d, tmp_vec_d);

	cudaError_t cuerr;
	cuerr = cudaGetLastError();
	if (cuerr != cudaSuccess) {
		cout << "ERROR1!" << cudaGetErrorString(cuerr) << endl;
	}

	int* gist_res = new int[pow2k];

	cudaMemcpy(gist, gist_d, (pow2k) * sizeof(int), cudaMemcpyDeviceToHost);

	cuerr = cudaGetLastError();
	if (cuerr != cudaSuccess) {
		cout << "ERROR2!" << cudaGetErrorString(cuerr) << endl;
	}

	double t2 = omp_get_wtime();

	cout << t2 - t1 << endl;
	ofstream fout;
	fout.open("out.txt");

	int check = 0;

	int* res = new int[N + 1];

	for (int i = 0; i < N + 1; i++)
		res[i] = 0;

	for (int i = 0; i < pow2k; i++)
	{
		res[gist[i]] += 1;
	}

	for (int i = 0; i < N + 1; i++)
	{	if (res[i] != 0)
			fout << i << "\t" << res[i] << endl;
		check += res[i];
	}
	//cout << check;
	fout.close();
}

