#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>
#include "mpi.h"
#include "timing.h"

typedef struct vector {
  int *indices;
  double *values;
  int length;
} Vector;

// function headers
void print_ints(int *array, int length);
void print_ints(int *array, int length);
Vector *generate_matrix(int n);
Vector *serial(Vector *columns);
double *dot_product(Vector *col, Vector *row);
void get_counts(int *indices, int size, int *send_counts); // send_counts is buffer
void get_displacements(int *send_counts, int size, int *displacements); // displacements is buffer
void compare(Vector *serial, Vector *parallel);

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