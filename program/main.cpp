#include <iostream>
#include <ctime>
#include "inoutput.h"
#include "operations.h"


int main(int argc, char** argv){
	int n, m, k; // the size of the matrix, the size of the printed part of the matrix, the flag
	char* filename; // the filename of input file
	double** A, ** b, ** x, ** temp; // the matrix (n * n), the right part of the equation (n * 1) and the solution
	
	if (argc < 4){
		std::cout << "There are too few arguments";
		return 0;
	}
        //try {
		n = std::stoi(argv[1]);
		m = std::stoi(argv[2]);
		k = std::stoi(argv[3]);
        /*}
	catch (std::invalid_argument){
		std::cout << "There is wrong type of some arguments";
		return 0;
        }*/
	
	A = new double*[n];
	b = new double*[n];
	x = new double*[n];
	temp = new double*[n];
	for (int i = 0; i < n; i++){
		A[i] = new double[n];
		b[i] = new double[1];
		x[i] = new double[1];
		temp[i] = new double[1];
	}
	
	clock_t t = clock();
	if (k == 0){
		if (argc < 5){
			std::cout << "There are too few arguments (filename is absent)";
			return 0;
		}
		else{
			filename = argv[4];
			switch (input_from_file(A, n, filename)){
				case 1: std::cout << "There is wrong data in the file";
						return 0;
				case 2: std::cout << "The file wasn't found";
						return 0;
			}
		}
	} else {
		if (initialization_of_the_matrix_by_the_formula(A, n, k)){
			std::cout << "There is wrong value of the argument k";
			return 0; 
		}
		
	}
	initialization_of_the_right_part(A, b, n);
	
	printing_the_matrix(A, m, m);
	printing_the_matrix(b, m, 1);
	//clock_t t = clock();
	function(A, b, x, temp, n);
	t = clock() - t;
	printing_the_matrix(A, m, m);
	printing_the_matrix(b, m, 1);
	std::cout << std::endl;
	printing_the_matrix(x, m, 1);
	
	matrix_product(A, x, temp, n, n, n, 1);
	for (int i = 0; i < n; i++){
		temp[i][0] -= b[i][0];
	}
	
	printf("\nNorm of nevyazki %10.3e \n", norm_of_the_vector(temp, n, 0) / norm_of_the_vector(b, n, 0));
	std::cout << "Time of algorithm working: " << double(t) / 1000;
	
	matrix_product(A, x, temp, n, n, n, 1);
	for (int i = 0; i < n; i++){
		if (i % 2){
			temp[i][0] = x[i][0];
		} else {
			temp[i][0] = x[i][0] - 1;
		}
	}
	printf("\nNorm of pogreshnosti: %10.3e\n", norm_of_the_vector(temp, n, 0));
	
	
	return 0;
}
