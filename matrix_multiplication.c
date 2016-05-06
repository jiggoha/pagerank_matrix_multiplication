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
  int size;
} Vector;

// function headers
void print_ints(int *array, int length);
void print_ints(int *array, int length);
Vector *generate(int n);
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
