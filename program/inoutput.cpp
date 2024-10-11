#include "inoutput.h"

int input_from_file(double** A, int n, char* filename){
	FILE* input = fopen(filename, "r");
	if (!input){
		return 2;
	}
	double temp;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			if (fscanf(input, "%lf", &temp) != 1){
				return 1;
			}
			A[i][j] = temp;
		}
	}
	fclose(input);
	return 0;
}

int initialization_of_the_matrix_by_the_formula(double** A, int n, int k){
	if (k == 1){
		for (int i = 0; i < n; i++){
			for (int j = 0; j < n; j++){
				A[i][j] = n - std::max(i + 1, j + 1) + 1;
			}
		}
	} else if (k == 2){
		for (int i = 0; i < n; i++){
			for (int j = 0; j < n; j++){
				A[i][j] = std::max(i + 1, j + 1);
			}
		}
	} else if (k == 3){
		for (int i = 0; i < n; i++){
			for (int j = 0; j < n; j++){
				A[i][j] = std::abs(i - j);
			}
		}
	} else if (k == 4){
		for (int i = 0; i < n; i++){
			for (int j = 0; j < n; j++){
				A[i][j] = 1 / double(i + j + 1);
			}
		}
	} else {
		return 1;
	}
	return 0;
}

void initialization_of_the_right_part(double** A, double** b, int n){
	for (int i = 0; i < n; i++){
		b[i][0] = 0;
		for (int k = 0; k < n; k += 2){
			b[i][0] += A[i][k];
		}
	}
}

void printing_the_matrix(double** A, int n, int m){
    for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                    printf("%10.3e ", A[i][j]);
            }
            printf("\n");
    }
    printf("\n");
}
