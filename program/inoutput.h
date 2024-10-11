#include <iostream>
#include <fstream>
#include <cstdio>

int input_from_file(double** A, int n, char* filename);
int initialization_of_the_matrix_by_the_formula(double** A, int n, int k);
void initialization_of_the_right_part(double** A, double** b, int n);
void printing_the_matrix(double** A, int n, int m);
