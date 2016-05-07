#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>
#include "mpi.h"
#include "timing.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct vector {
  int *indices;
  double *values;
  int length;
} Vector;

// function headers
double *dot_product(Vector *col, Vector *row);
void get_counts(int *indices, int size, int *send_counts); // TODO, send_counts is buffer
void get_displacements(int *send_counts, int size, int *displacements); // TODO, displacements is buffer

// initializations
Vector *generate_vector(int n); // TODO: fix randomness
void destroy_vector(Vector *vec);
Vector **generate_matrix(int n);
void destroy_matrix(Vector **matrix, int n);

// debugging
Vector **serial(Vector *columns); // TODO
int are_matrices_same(Vector **serial, Vector **parallel);
void print_vector(Vector *vec);
void print_ints(int *array, int length);
void print_doubles(double *array, int length);

int main (int argc, char **argv) {
  srand(12345);

  MPI_Init(&argc, &argv);

  // number of processes
  int num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  // rank of process
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Finalize();
  return 0;
}

void print_vector(Vector *vec) {
  for (int i = 0; i < vec->length; i++) {
    printf("vec[%d] = %f\n", vec->indices[i], vec->values[i]);
  }
  printf("\n");
}

void print_ints(int *array, int length) {
  for (int i = 0; i < length; i++) {
    printf("%d ", array[i]);
  }
  printf("\n");
}

void print_doubles(double *array, int length) {
  for (int i = 0; i < length; i++) {
    printf("%f ", array[i]);
  }
  printf("\n");
}

Vector *generate_vector(int n) {
  Vector *vec = malloc(sizeof(Vector));
  
  if (n == 0) {
    vec->length = 0;
  } else {
    vec->length = rand() % n;
    
    if (vec->length != 0) { 
      vec->indices = malloc(sizeof(int) * vec->length);
      vec->values = malloc(sizeof(double) * vec->length);

      for (int i = 0; i < vec->length; i++) {
        vec->indices[i] = i; // TODO: random list of increasing indices?
        vec->values[i] = rand() / RAND_MAX;
      }
    }
  }

  return vec;
}

Vector **generate_matrix(int n) {
  Vector **matrix = malloc(n * sizeof(Vector *));

  // randomly generate column
  for (int i = 0; i < n; i++) {
    int length = rand() % n;
    matrix[i] = generate_vector(length);
  }

  return matrix;
}

void destroy_vector(Vector *vec) {
  if (vec->length != 0) {
    free(vec->indices);
    free(vec->values);
  }

  free(vec);
}

void destroy_matrix(Vector **matrix, int n) {
  for (int i = 0; i < n; i++) {
    destroy_vector(matrix[i]);
  }

  free(matrix);
}

int are_matrices_same(Vector **a, Vector **b, int n) {
  for (int i = 0; i < n; i++) {
    if (a[i]->length != b[i]->length) {
      printf("row %d length doesn't match: %d, %d\n", i, a[i]->length, b[i]->length);
      return 0;
    }

    for (int j = 0; j < a[i]->length; j++) {
      // TODO: change error constant?
      if (a[i]->indices[j] != b[i]->indices[j] || a[i]->values[j] != b[i]->values[j] > 0.001) {
        printf("different at element %d, %d\n", i, j);
        return 0;
      }
    }
  }

  return 1;
}

double dot_product(Vector *a, Vector *b) {
  double result = 0.;
  int index_a = 0;
  int index_b = 0;

  while (index_a < a->length  && index_b < b->length) {
    if (a->indices[index_a] == b->indices[index_b]) {
      result += a->values[index_a] * b->values[index_b];
      index_a++;
      index_b++;
    } else if (a->indices[index_a] > b->indices[index_b]) {
      index_b++;
    } else {
      index_a++;
    }
  }

  return result;
}