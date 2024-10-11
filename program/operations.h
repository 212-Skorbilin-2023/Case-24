#include <iostream>
#include <cmath>

int matrix_product(double** A, double** B, double** C, int n, int m, int p, int q);
int reflection_matrix_product(double** A, double** B, int n, int p, int q, int k, int flag);
double norm_of_the_vector(double** a, int n, int k);
void reverse_Gauss_method(double** A, double** b, double** x, int n);
void function(double** A, double** b, double** x, double** temp, int n);
