#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include <set>

using namespace std;

string XOR(string a, string b, int N)
{
	string c(N, '0');
	for (int i = 0; i < N; i++)
	{
		if (a[i] != b[i])
			c[i] = '1';
		else
			c[i] = '0';
	}

	return c;
}


int weigth(string a, int N)
{
	int cnt = 0;
	for (int i = 0; i < N; i++)
		if (a[i] == '1')
			cnt++;
	return cnt;
}


int main()
{
	ifstream fin;
	fin.open("in.txt");
	
	int N, K;
	fin >> N >> K;

	string * A = new string[K];
	for (int i = 0; i < K; i++)
		fin >> A[i];

	set <string> res;
	string zero(N, '0');

	int pow2k = pow(2, K);

	int temp_i, mult;
	string tmp;
	for (int i = 0; i < pow2k; i++)
	{
		mult = i % 2;
		temp_i = i / 2;
		tmp = mult ? A[0] : zero;

		for (int j = 1; j < K; j++)
		{
			mult = temp_i % 2;
			temp_i /= 2;
			if (mult)
				tmp = XOR(tmp, A[j], N);
		}
		res.insert(tmp);
	}

	int* gist = new int[N + 1];
	for (int i = 0; i < N + 1; i++)
		gist[i] = 0;

	set <string> ::iterator iter = res.begin();

	while (iter != res.end())
	{
		gist[weigth(*iter, N)] += 1;
		iter++;
	}
	ofstream fout;
	fout.open("out.txt");


	for (int i = 0; i < N + 1; i++)
	{
		if (gist[i] != 0)
			fout << i << '\t' << gist[i] << endl;
	}
}

