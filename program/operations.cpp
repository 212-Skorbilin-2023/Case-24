#include "operations.h"
#include "inoutput.h"

int matrix_product(double** A, double** B, double** C, int n, int m, int p, int q){
	if (m != p){
		return 1;
	}
	for (int i = 0; i < n; i++){
		for (int j = 0; j < q; j++){
			C[i][j] = 0;
			for (int k = 0; k < m; k++){
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return 0;
}

// A is the vector forming reflection matrix
// n and m are sizes of reflection matrix so they are equal
int reflection_matrix_product(double** A, double** B, int n, int p, int q, int k, int flag){
	if (n != p){
		return 1;
	}
	double temp;
	if (flag){
		for (int j = k; j < q; j++){
			temp = 0;
			for (int i = k; i < n; i++){
				temp += A[i][0] * B[i][j];
			}
			temp *= 2;
			for (int i = k; i < n; i++){
				B[i][j] -= temp * A[i][0];
			}
		}
	} else {
		for (int j = 0; j < q; j++){
			temp = 0;
			for (int i = k; i < n; i++){
				temp += A[i][0] * B[i][j];
			}
			temp *= 2;
			for (int i = k; i < n; i++){
				B[i][j] -= temp * A[i][0];
			}
		}
	}
	
	return 0;
}

double norm_of_the_vector(double** a, int n, int k){ // k is a number of the column
	double sum = 0;
	for (int i = 0; i < n; i++){
		sum += a[i][k] * a[i][k];
	}
	return sqrt(sum);
}

void reverse_Gauss_method(double** A, double** b, double** x, int n){
	for (int i = n - 1; i >= 0; i--){
		x[i][0] = b[i][0];
		for (int j = i + 1; j < n; j++){
			x[i][0] -= A[i][j] * x[j][0];
		}
		x[i][0] /= A[i][i];
	}
}

void function(double** A, double** b, double** x, double** temp, int n){
	double s, norm_a, norm_x;
	for (int k = 0; k < n - 1; k++){
		s = 0;
		for (int j = k + 1; j < n; j++){
			s += A[j][k] * A[j][k];
		}
		norm_a = sqrt(A[k][k] * A[k][k] + s);
		temp[k][0] = A[k][k] - norm_a;
		for (int j = k + 1; j < n; j++){
			temp[j][0] = A[j][k];
		}
		norm_x = sqrt(temp[k][0] * temp[k][0] + s);
		for (int j = k; j < n; j++){
			temp[j][0] /= norm_x;
		}
		reflection_matrix_product(temp, A, n, n, n, k, 1);
		reflection_matrix_product(temp, b, n, n, 1, k, 0);
	}
	reverse_Gauss_method(A, b, x, n);
}
