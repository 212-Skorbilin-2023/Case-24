#include <iostream>
#include <cmath>
#include <pthread.h>
#include <stdio.h>

struct tMatr{
	double** A;
	double** B;
	double** C;
	int n, m, q, i, count;
};

struct dMatr{
	double** A;
	double** B;
	int n, q, k, j, count;
};

struct Matr{
	double** A;
	double result;
	int n, k, i, s, count;
};

struct ttMatr{
	double** A;
	double** x;
	double result;
	int n, i, ii, count;
};

void* refl(void* arg);
void* rev(void* arg);
void* norm(void* arg);
void* prod(void* arg);

int matrix_product(double** A, double** B, double** C, int n, int m, int p, int q, int count);
int reflection_matrix_product(double** A, double** B, int n, int p, int q, int k, int flag, int count);
double norm_of_the_vector(double** a, int n, int k, int start, int count);
void reverse_Gauss_method(double** A, double** b, double** x, int n, int count);
void function(double** A, double** b, double** x, double** temp, int n, int count);
