#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fenv.h>
#include <time.h>
#include <sys/time.h>
#include "utils.h"

#define _POSIX_C_SOURCE 199309L


double** inicializarmatriz(int ordem){
    double** matriz;
    int j = 0;

    matriz = (double **)malloc(ordem * sizeof(double*)); 
    if(!matriz){
        printf("Erro ao alocar memoria\n");
        return NULL;
    }

    for(j = 0; j<ordem; j++){
        matriz[j] = (double *)malloc(sizeof(double) * ordem); 
        if (!(matriz[j])){
            printf("Erro ao alocar memoria\n");
            return NULL;
        }
    }

    return matriz;
}

void zeramatriz(double **matriz, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matriz[i][j] = 0.0;
        }
    }
}

void zerarvetor(double *v, int N){
    for(int i = 0; i < N; i++ )
        v[i] = 0.0;
}

double* inicializavetor(int ordem){
    double *vetor;
    vetor = (double*) malloc(ordem * sizeof(double)); 
    if(!vetor){
        printf("Erro ao alocar memoria\n");
        return NULL;
    }
    return vetor;
}

double* lelinha(int ordem, double* linha){
    for (int i = 0; i< ordem; i++)
        scanf("%lf", &linha[i]);
    
    return linha;
}

void funcaoleitura(double** matriz, int ordem){
    for (int i = 0; i < ordem; i++)
        lelinha(ordem, matriz[i]);
}

void imprimevetor(int ordem, double* vetor){
    for (int i = 0; i < ordem; i++)
        printf("%.15e ", vetor[i]);
    printf("\n");
}

void funcaoimpressao(double** matriz, int ordem){
    for (int i = 0; i<ordem; i++)
        imprimevetor(ordem, matriz[i]);
}

int decomposicaoLU(double** A, double** L, double** U, int N) {
    fesetround(FE_DOWNWARD);

    for (int i = 0; i < N; i++) {

        for (int k = i; k < N; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++) 
                sum += (L[i][j] * U[j][k]);
            U[i][k] = A[i][k] - sum;
        }

        if (U[i][i] == 0){
            printf("Essa matriz não tem inversa\n");
            return 0; 
        }

        for (int k = i; k < N; k++) {
            if (i == k) 
                L[i][i] = 1; 
            else {
                double sum = 0;
                for (int j = 0; j < i; j++) 
                    sum += (L[k][j] * U[j][i]);

                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }

    }
    return 1; 
}

void primeirasubstituicao(double** L, double* b, double* y, int N) {
    fesetround(FE_DOWNWARD);

    for (int i = 0; i < N; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
    }
}

void segundasubstituicao(double** U, double* y, double* x, int N) {
    fesetround(FE_DOWNWARD);

    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}

int calcularinversa(double **A, double **L, double **U, double **invA, int N) {
    fesetround(FE_DOWNWARD); 

    if (!decomposicaoLU(A, L, U, N)) {
        printf("A matriz nao pode ser decomposta.\n");
        return 0;
    }

    double *e = inicializavetor(N);
    double *y = inicializavetor(N);
    double *x = inicializavetor(N);

    for (int i = 0; i < N; i++) {
        zerarvetor(e, N);
        e[i] = 1; 

        primeirasubstituicao(L, e, y, N);
        segundasubstituicao(U, y, x, N);

        for (int j = 0; j < N; j++) {
            invA[j][i] = x[j];
        }
    }
    free(e);
    free(y);
    free(x);
    return 1;
}

void freematriz(double** matriz, int ordem) {
    for (int i = 0; i < ordem; i++) {
        free(matriz[i]); 
    }
    free(matriz); 
}

void multiplicamatrizes(double **A, double **B, double **C, int n) {
    int i, j, k;
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void matrizAAmenosI(double** A, int N) {
    for (int i = 0; i < N; i++) 
            A[i][i] = A[i][i] - 1.0; 
}


double calculonorma(double **A, int N) {
    double somona = 0.0;
    
    for (int i = 0; i < N; i++) {
        double soma = 0.0; 
        for (int j = 0; j < N; j++) 
            soma += A[j][i] * A[j][i]; 
        somona += sqrt(soma); 
    }
    
    return somona / N; 
}

/*rtime_t timestamp (void)
{
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
  return ( (rtime_t) tp.tv_sec*1.0e3 + (rtime_t) tp.tv_nsec*1.0e-6 );
}*/


int main (){
    int ordem;
    double **matrizA, **matrizU, **matrizL, 
    **matriz_invA, **matrizAA;
    double norma;
    scanf("%d", &ordem);
    rtime_t tempo;

    matrizA = inicializarmatriz(ordem);
    matrizU = inicializarmatriz(ordem);
    matrizL = inicializarmatriz(ordem);
    matriz_invA = inicializarmatriz(ordem);

    zeramatriz(matrizU, ordem);
    zeramatriz(matrizL, ordem);

    funcaoleitura(matrizA, ordem);

    tempo = timestamp();
    printf("Vamos calcular a matriz invA\n");
    if (!calcularinversa(matrizA, matrizL, matrizU, matriz_invA, ordem))
        return 0;
    tempo = timestamp() - tempo;
    printf("Tempo: %.15e\n", tempo);
    
    printf("Matriz Inv A:\n");
    funcaoimpressao(matriz_invA, ordem);
    printf("\n");

    matrizAA = inicializarmatriz(ordem);
    zeramatriz(matrizAA, ordem);

    multiplicamatrizes(matrizA, matriz_invA, matrizAA, ordem);
    printf("Matriz AA:\n");
    funcaoimpressao(matrizAA, ordem);
    printf("\n");

    matrizAAmenosI(matrizAA, ordem);
    printf("Matriz AA - I :\n");
    funcaoimpressao(matrizAA, ordem);
    printf("\n");

    norma = calculonorma(matrizAA, ordem);
    printf("Calculo norma = %.8e\n", norma);
    
    freematriz(matrizA , ordem);
    freematriz(matrizL , ordem);
    freematriz(matrizU , ordem);
    freematriz(matriz_invA , ordem);
    freematriz(matrizAA, ordem);

    return 0;
}
