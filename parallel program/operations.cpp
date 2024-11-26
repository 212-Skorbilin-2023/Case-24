#include "operations.h"
#include "inoutput.h"

void* refl(void* arg){
	struct dMatr* dM = (struct dMatr*)arg;
	double temp;
        for (int J = dM->j; J < dM->q; J += dM->count){
		temp = 0;
		for (int i = dM->k; i < dM->n; i++){
			temp += dM->A[i][0] * dM->B[i][J];
		}
		temp *= 2;
		for (int i = dM->k; i < dM->n; i++){
			dM->B[i][J] -= temp * dM->A[i][0];
		}
	}
	
	pthread_exit(NULL);
}

void* rev(void* arg){
        struct ttMatr* ttM = (struct ttMatr*)arg;
        double temp = 0;
        //std::cout << "3:) " << ttM->i << "\n";
        for (int J = ttM->i; J < ttM->n; J += ttM->count){
            //std::cout << "3:) " << ttM->i << " " << J << "\n";
            temp += ttM->A[ttM->ii][J] * ttM->x[J][0];
        }
        ttM->result = temp;
        pthread_exit(NULL);
}

void* norm(void* arg){
	struct Matr* M = (struct Matr*)arg;
	double sum = 0;
        for (int J = M->s + M->i; J < M->n; J += M->count){
		sum += M->A[J][M->k] * M->A[J][M->k];
	}
	M->result = sum;
	pthread_exit(NULL);
}

void *prod(void* arg){
	struct tMatr* tM = (struct tMatr*)arg;
        for (int i = tM->i; i < tM->n; i += tM->count){
		for (int j = 0; j < tM->q; j++){
			tM->C[i][j] = 0;
			for (int k = 0; k < tM->m; k++){
				tM->C[i][j] += tM->A[i][k] * tM->B[k][j];
			}
		}
	}
	pthread_exit(NULL);
}

int matrix_product(double** A, double** B, double** C, int n, int m, int p, int q, int count){
	if (m != p){
		return 1;
	}
	
        pthread_t* id = new pthread_t[count];
        for (int i = 0; i < count; i++){
		id[i] = 0;
	}
        struct tMatr* triple = new struct tMatr[count];
        for (int i = 0; (i < n) && (i < count); i++){
		triple[i].A = A;
		triple[i].B = B;
		triple[i].C = C;
		triple[i].n = n;
		triple[i].m = m;
		triple[i].q = q;
		triple[i].i = i;
                triple[i].count = count;
	}
        for (int i = 0; (i < n) && (i < count); i++){
		pthread_create(id + i, NULL, *prod, (void*)(triple + i));
	}
        for (int i = 0; (i < n) && (i < count); i++){
		pthread_join(id[i], NULL);
	}
	delete[] triple;
        delete[] id;
	return 0;
}

// A is the vector forming reflection matrix
// n and m are sizes of reflection matrix so they are equal
int reflection_matrix_product(double** A, double** B, int n, int p, int q, int k, int flag, int count){
	if (n != p){
		return 1;
	}
	
	double temp;
	
	if (flag){
                pthread_t* id = new pthread_t[count];
                for (int i = 0; i < count; i++){
			id[i] = 0;
		}
                struct dMatr* couple = new struct dMatr[count];
                for (int i = 0; (i < count) && (k + i < q); i++){
			couple[i].A = A;
			couple[i].B = B;
			couple[i].n = n;
			couple[i].q = q;
			couple[i].k = k;
                        couple[i].count = count;
			couple[i].j = k + i;
		}
                for (int j = 0; (j < q - k) && (j < count); j++){
			pthread_create(id + j, NULL, *refl, (void*)(couple + j));
		}
                for (int i = 0; (i < count) && (k + i < q); i++){
			pthread_join(id[i], NULL);
		}
		delete[] couple;
                delete[] id;
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

double norm_of_the_vector(double** a, int n, int k, int start, int count){ // k is a number of the column
        pthread_t* id = new pthread_t[count];
        for (int i = 0; i < count; i++){
		id[i] = 0;
	}
	double sum = 0;
        struct Matr* once = new struct Matr[count];
        for (int i = 0; (i < count) && (i < n - start); i++){
		once[i].A = a;
		once[i].n = n;
		once[i].k = k;
		once[i].s = start;
		once[i].i = i;
                once[i].count = count;
                once[i].result = 0;
	}
        for (int i = 0; (i < count) && (i < n - start); i++){
		pthread_create(id + i, NULL, *norm, (void*)(once + i));
	}
        for (int i = 0; (i < count) && (i < n); i++){
            if (id[i]) pthread_join(id[i], NULL);
	}
        for (int i = 0; (i < count) && (i < n); i++){
                if (id[i]) sum += once[i].result;
	}
	delete[] once;
        delete[] id;
	return sum;
}

void reverse_Gauss_method(double** A, double** b, double** x, int n, int count){
    struct ttMatr* ttriple = new struct ttMatr[count];
    pthread_t* id = new pthread_t[count];
    for (int ii = n - 1; ii >= 0; ii--){
        for (int i = 0; i < count; i++){
                id[i] = 0;
        }
        //std::cout << "1:)\n";
        for (int i = 0; (ii + i + 1 < n) && (i < count); i++){
                ttriple[i].A = A;
                ttriple[i].x = x;
                ttriple[i].n = n;
                ttriple[i].ii = ii;
                ttriple[i].i = ii + i + 1;
                ttriple[i].result = 0;
                ttriple[i].count = count;
        }
        //std::cout << "2:)\n";
        for (int i = 0; (ii + i + 1 < n) && (i < count); i++){
                pthread_create(id + i, NULL, *rev, (void*)(ttriple + i));
        }
        //std::cout << "3:)\n";
        for (int i = 0; (ii + i + 1 < n) && (i < count); i++){
            pthread_join(id[i], NULL);
        }
        //std::cout << "4:)\n";

        x[ii][0] = b[ii][0];
        for (int i = 0; (ii + i + 1 < n) && (i < count); i++){
                x[ii][0] -= ttriple[i].result;
        }
        x[ii][0] /= A[ii][ii];
        //std::cout << "5:)\n";
    }
    delete[] ttriple;
    delete[] id;
}

void function(double** A, double** b, double** x, double** temp, int n, int count){
	double s, norm_a, norm_x;
	for (int k = 0; k < n - 1; k++){
            s = norm_of_the_vector(A, n, k, k + 1, count);
            norm_a = sqrt(A[k][k] * A[k][k] + s);
            temp[k][0] = A[k][k] - norm_a;
            for (int j = k + 1; j < n; j++){
                    temp[j][0] = A[j][k];
            }
            norm_x = sqrt(temp[k][0] * temp[k][0] + s);
            for (int j = k; j < n; j++){
                    temp[j][0] /= norm_x;
            }
            reflection_matrix_product(temp, A, n, n, n, k, 1, count);
            reflection_matrix_product(temp, b, n, n, 1, k, 0, count);
	}
        reverse_Gauss_method(A, b, x, n, count);
}
