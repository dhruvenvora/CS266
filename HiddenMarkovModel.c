/*
 ============================================================================
 Name        : CS266_Assignment1.c
 Author      : Dhruven
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void print_array(int a, int b, double **x, short isTranspose) {
	double sum;
	if(isTranspose == 0){
		for (int i = 0; i < a; i++) {
			sum = 0;
			for (int j = 0; j < b; j++) {
				printf("%lf ", x[i][j]);
				sum += x[i][j];
			}
			printf(" --> TOTAL = %lf\n", sum);
		}
	}
	else {
		for (int i = 0; i < b; i++) {
			for (int j = 0; j < a; j++) {
				printf("%lf    ", x[j][i]);
			}
			printf("\n");
		}
	}
}

int main(int argc, char *argv[]) {
	time_t t;
	int N = 2;
	int T = 2000;
	int M = 27; // take M = 27 if one want to include space.
	int MAX_ITER = 400;

	if(argc == 6){
		char* file_path = argv[1];
		sscanf(argv[2], "%d", &N);
		sscanf(argv[3], "%d", &M);
		sscanf(argv[4], "%d", &T);
		sscanf(argv[5], "%d", &MAX_ITER);


		//printf("Enter values of N, M, T : ");
		//scanf("%d %d %d", &N, &M, &T);

		double **A = malloc(sizeof *A * N);
		for(int i=0;i<N;i++)
			A[i] = malloc(sizeof *A[i] * N);

		double **B = malloc(sizeof *B * N);
		for(int i=0;i<N;i++)
			B[i] = malloc(sizeof *B[i] * M);

		double *pi = malloc(sizeof *pi * N);

		int *O = malloc(sizeof *O * T);

		double **alpha = malloc(sizeof *alpha * T);

		for(int i=0;i<T;i++)
			alpha[i] = malloc(sizeof *alpha[i] * N);

		double **beta = malloc(sizeof *beta * T);
		for(int i=0;i<T;i++)
			beta[i] = malloc(sizeof *beta[i] * N);

		double **gamma = malloc(sizeof *gamma * T);
		for(int i=0;i<T;i++)
			gamma[i] = malloc(sizeof *gamma[i] * N);

		double ***gamma2 = malloc(sizeof **gamma2 * T);

		for(int t=0;t<T;t++){
			gamma2[t] = malloc(sizeof *gamma2[t] * N);
			for(int i=0;i<N;i++)
				gamma2[t][i] = malloc(sizeof *gamma2[t][i] * N);
		}

		/*
		 * Reading file for observations
		 */
		FILE *file;
		file = fopen(file_path, "r");
		if (!file)
			return 0;
		char d;
		int i=0;
		while(i<T){
			d = fgetc(file);
			if(d >= 'A' && d <= 'Z')
				O[i] = (int)d + 32 - 97; // First convert the capital alphabet to small one and then normalize it to index of array.
			else if(d >= 'a' && d <= 'z')
				O[i] = (int)d - 97;
			else if(d == ' ' || d == '\n')
				O[i] = M-1;
			else continue; //Don't allow special characters like punctuation marks, numbers and other symbols.
			//printf("%c", d);
			i++;//increament i only when character is alphabet or space.
		}
		fclose(file);

		srand((unsigned) time(&t));// seed value to generate random values

		/*
		 * Assigning initial probabilities. All equal.
		 * TODO : make probabilities unequal.
		 */
		double p = 1.0 / N;
		double sum = 0;
		for (int i = 0; i < N; i++) {
			sum = 0;
			for (int j=0;j<N-1;j++) {
				double r = (double)(rand()%10)/(double)RAND_MAX;
				A[i][j] = p - r;
				sum += r;
			}
			A[i][N-1] = p + sum;
		}

		sum = 0;
		for (int i = 0; i < N; i++) {
			double r = (double)(rand()%100)/(double)RAND_MAX;
			pi[i] = p - r;
			sum+=r;
		}
		pi[N-1] = p+sum;

		p = 1.0 / M;
		for (int i = 0; i < N; i++) {
			sum = 0;
			for (int j = 0; j < M-1; j++) {
				double r = (double)(rand()%10)/(double)RAND_MAX;
				B[i][j] = p-r;
				sum += r;
			}
			B[i][M-1] = p + sum;
		}

		printf("1. Initialized......\n");
		printf("2. N, M, T values are %d %d %d......\n", N, M, T);

		printf("\nA MATRIX---------------------------------------------\n");
		print_array(N, N, A, 0);

		printf("\nB MATRIX---------------------------------------------\n");
		print_array(N, M, B, 0);

		printf("\nPI MATRIX---------------------------------------------\n");
		for(int i=0;i<N;i++)
			printf("%lf  ", pi[i]);
		printf("\n");

		int ITER = 1;
		double logProb = -999999998;
		double oldLogProb = -999999999;

		double c[T]; //Scaling factor

		while(ITER < MAX_ITER && logProb > oldLogProb){
			oldLogProb = logProb;
			/*
			 * ALPHA PASS---------------------------------------------------------------------------------------
			 */
			c[0] = 0;

			for (int i = 0; i < N; i++) {
				alpha[0][i] = pi[i] * B[i][O[0]];
				c[0] = c[0] + alpha[0][i];
			}
			c[0] = 1.0 / c[0];
			for (int i = 0; i < N; i++) {
				alpha[0][i] = c[0] * alpha[0][i]; // Normalization of ALPHA.
			}

			for (int t = 1; t < T; t++) {
				c[t] = 0;
				for (int i = 0; i < N; i++) {
					alpha[t][i] = 0;
					for (int j = 0; j < N; j++) {
						alpha[t][i] += alpha[t - 1][j] * A[j][i];
					}
					alpha[t][i] *= B[i][O[t]];
					c[t] += alpha[t][i];
				}

				c[t] = 1.0 / c[t];
				for (int i = 0; i < N; i++) {
					alpha[t][i] *= c[t];
				}
			}

			//printf("\n3. ALPHA Pass complete.---------------------------------------------\n");
			//print_array(T, N, alpha, 0);

			/*
			 * BETA PASS-----------------------------------------------------------------------------------------
			 */
			for (int i = 0; i < N; i++) {
				beta[T - 1][i] = c[T - 1];
			}

			for (int t = T - 2; t >= 0; t--) {
				for (int i = 0; i < N; i++) {
					beta[t][i] = 0;
					for (int j = 0; j < N; j++) {
						beta[t][i] += A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
					}
					beta[t][i] *= c[t]; // Normalization of beta.
				}
			}

			//printf("\n4. BETA pass complete---------------------------------------------\n");
			//print_array(T, N, beta);

			/*
			 * GAMMA PASS----------------------------------------------------------------------------------------
			 */
			double denom;
			for (int t = 0; t < T - 1; t++) {
				denom = 0;
				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						denom += alpha[t][i] * A[i][j] * B[j][O[t + 1]]
								* beta[t + 1][j];
					}
				}

				for (int i = 0; i < N; i++) {
					gamma[t][i] = 0;
					for (int j = 0; j < N; j++) {
						gamma2[t][i][j] = (alpha[t][i] * A[i][j] * B[j][O[t + 1]]
								* beta[t + 1][j]) / denom;
						gamma[t][i] += gamma2[t][i][j];
					}
				}
			}

			denom = 0;
			for (int i = 0; i < N; i++) {
				denom += alpha[T - 1][i];
			}
			for (int i = 0; i < N; i++) {
				gamma[T - 1][i] = alpha[T - 1][i] / denom;
			}
			//printf("\nGAMMA---------------------------------------------\n");
			//print_array(T, N, gamma);

			/*
			 * RE-ESTIMATION------------------------------------------------------------------------------------
			 */
			for (int i = 0; i < N; i++) {
				pi[i] = gamma[0][i]; // Re-estimation of PI.
			}

			double numer = 0;
			denom = 0;
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					numer = 0;
					denom = 0;
					for (int t = 0; t < T - 1; t++) {
						numer += gamma2[t][i][j];
						denom += gamma[t][i];
					}
					A[i][j] = numer / denom;
				}
			}

			for (int i = 0; i < N; i++) {
				for (int j = 0; j < M; j++) {
					numer = 0;
					denom = 0;
					for (int t = 0; t < T; t++) {
						if (O[t] == j)
							numer += gamma[t][i];
						denom += gamma[t][i];
					}
					B[i][j] = numer / denom;
				}
			}
			//printf("\n%d Re-estimation complete.---------------------------------------------\n", ITER);
			/*
			 * COMPUTING LOG PROBABILITY-------------------------------------------------------------------
			 */

			//printf("\nA MATRIX---------------------------------------------\n");
			//print_array(N, N, A);

			logProb = 0;
			for (int i = 0; i < T; i++) {
				logProb += log(c[i]);
			}
			logProb = -logProb;
			printf("PASS %d = %lf\n", ITER, logProb);
			ITER++;
		}

		printf("\n6. HMM Model updated....\n");
		printf("\nUPDATED A MATRIX---------------------------------------------\n");
		print_array(N, N, A, 0);
		printf("\nUPDATED B MATRIX---------------------------------------------\n");
		print_array(N, M, B, 1);
	}
	return 0;
}
