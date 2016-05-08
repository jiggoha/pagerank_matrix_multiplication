#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>
#include "mpi.h"
#include "timing.h"
#include "random_list.h"

#define DEBUG (0)

#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct vector {
  int *indices;
  double *values;
  int length;
} Vector;

typedef struct matrix {
  Vector **vectors;
  int n;
} Matrix;

// function headers
double *dot_product(Vector *col, Vector *row);
void get_counts(int *indices, int length, int *send_counts, int n, int num_procs); 
void get_displacements(int *send_counts, int *displacements, int num_procs); 

// initializations
Vector *generate_vector(int n, int length, int debug);
void destroy_vector(Vector *vec);
Matrix* newMatrix(int n);
Matrix *generate_matrix(int n, int debug);
void destroy_matrix(Matrix *matrix);

// debugging
Matrix *serial(Matrix *matrix, int p); // TODO
Matrix *transpose_representation(Matrix *matrix);

int are_matrices_same(Matrix *a, Matrix *b);
void print_vector(Vector *vec);
void print_matrix(Matrix *matrix);
void print_ints(int *array, int length);
void print_doubles(double *array, int length);

// global
int n;

int main (int argc, char **argv) {
  srand(12345);

  //get n from args
  n = atoi(argv[1])

  Matrix* matrix = generate_matrix(n, (DEBUG < 0));

  MPI_Init(&argc, &argv);

  // number of processes
  int num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  // rank of process
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Finalize();

  Matrix* serialResult = serial(matrix, p);

  printf("Are the matrices the same?\n");
  if(are_matrices_same(serialResult, parallelResult)) {
    printf("Yes!\n");
  } else {
    printf("No :(\n");
  }
  return 0;
}

void print_vector(Vector *vec) {
  for (int i = 0; i < vec->length; i++) {
    printf("vec[%d] = %f\n", vec->indices[i], vec->values[i]);
  }
  printf("\n");
}

void print_matrix(Matrix *matrix) {
  for (int i = 0; i < matrix->n; i++) {
    printf("column %d:\n", i);
    // printf("length: %d\n", matrix[i]->length);
    for (int j = 0; j < matrix->vectors[i]->length; j++) {
      printf("matrix[%d][%d] =  %f\n", i, matrix->vectors[i]->indices[j], matrix->vectors[i]->values[j]);
    }

    printf("\n");
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

Vector *generate_vector(int n, int length, int debug) {
  Vector *vec = malloc(sizeof(Vector));
  vec->length = length;

  if (length != 0) {
    vec->length = length;
    
    if (vec->length != 0) { 
      // make enough room for all n
      vec->indices = random_increasing_ints(n, length);
      vec->values = malloc(sizeof(double) * n);

      for (int i = 0; i < vec->length; i++) {
        if (debug) {
          vec->values[i] = i;
        } else {
          vec->values[i] = (double) rand() / RAND_MAX;
        }
      }
    }
  }

  return vec;
}

Matrix *generate_matrix(int n, int debug) {
  Matrix *matrix = malloc(sizeof(Matrix));
  matrix->vectors = malloc(sizeof(Vector *) * n);
  matrix->n = n;

  // randomly generate column
  for (int i = 0; i < n; i++) {
    int length = rand() % (n + 1);

    matrix->vectors[i] = generate_vector(n, length, debug);
  }

  return matrix;
}


Matrix* newMatrix(int n) {
  Matrix *new = malloc(sizeof(Matrix));
  new->n = n;
  new->vectors = malloc(sizeof(Vector *) * n);

  // memory stuff
  for (int i = 0; i < n; i++) {
    // actually doesn't need to be this long, but it's an upper bound
    new->vectors[i] = malloc(sizeof(Vector));
    new->vectors[i]->length = 0;
    new->vectors[i]->indices = malloc(sizeof(int) * n);
    new->vectors[i]->values = malloc(sizeof(double) * n);
  }

  return(new);
}

void destroy_vector(Vector *vec) {
  if (vec->length != 0) {
    free(vec->indices);
    free(vec->values);
  }

  free(vec);
}

void destroy_matrix(Matrix *matrix) {
  for (int i = 0; i < matrix->n; i++) {
    destroy_vector(matrix->vectors[i]);
  }

  free(matrix->vectors);
  free(matrix);
}

int are_matrices_same(Matrix *a, Matrix *b) {
  if (a->n != b-> n) {
    printf("not the same dimension\n");
    return 0;
  }

  for (int i = 0; i < a->n; i++) {
    if (a->vectors[i]->length != b->vectors[i]->length) {
      printf("row %d length doesn't match: %d, %d\n", i, a->vectors[i]->length, b->vectors[i]->length);
      return 0;
    }

    for (int j = 0; j < a->vectors[i]->length; j++) {
      // TODO: change error constant?
      if (a->vectors[i]->indices[j] != b->vectors[i]->indices[j]  
            || abs(a->vectors[i]->values[j] - b->vectors[i]->values[j]) > 0.001) {
        printf("different at element %d, %d by %f\n", i, j, a->vectors[i]->values[j] - b->vectors[i]->values[j]);
        return 0;
      }
    }
  }

  return 1;
}


Matrix *transpose_representation(Matrix *matrix) {

  Matrix* transposed = newMatrix(matrix->n);
  // fill it up
  for (int i = 0; i < matrix->n; i++) {
    for (int j = 0; j < matrix->vectors[i]->length; j++) {
      int new_i = matrix->vectors[i]->indices[j];
      double value = matrix->vectors[i]->values[j];

      int next_j = transposed->vectors[new_i]->length;
      transposed->vectors[new_i]->length++;

      transposed->vectors[new_i]->indices[next_j] = i;
      transposed->vectors[new_i]->values[next_j] = value;
    }
  }

  return transposed;
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

void get_counts(int *indices, int size, int *send_counts, int n, int num_procs) {

  int k = 0;
  for(int i = 0; i < num_procs; i++) {
    send_counts[i] = 0;
  }

  for(int i = 0; i < num_procs; i++) {
    while(indices[k] < (i+1) * n/num_procs && k < size) {
      send_counts[i]++;
      k++;
    }
  }
}

void get_displacements(int *send_counts, int *displacements, int num_procs) {

  displacements[0] = 0;
  for(int i = 1; i < num_procs; i++) {
    displacements[i] = displacements[i - 1] + send_counts[i - 1];
  }
}

Matrix *serial(Matrix *matrix_by_cols, int p) {

  Matrix* matrix_by_rows = transpose_representation(matrix);
  Matrix* temp_by_cols = newMatrix(matrix->n);
  Matrix* temp_by_rows = newMatrix(matrix->n);
  Matrix* temp;

  double res;
  int col_count;    //temporary variable to store growing length of column
  int row_counts[n];    //temporary array to store growing length of rows
  

  for(int i = 0; i < p; i++) {  //powers

    //init
    col_count = 0;
    for(int i = 0; i < n; i++) {
      row_counts[i] = 0;
    }

    //compute
    for(int j = 0; j < n; j++) {  //cols
      for(int k = 0; k < n; k++) {  //rows
        if(res = dot_product(matrix->vectors[j], matrix_by_rows->vectors[k]) != 0) {
          temp_by_cols->vectors[j]->indices[col_count] = k;
          temp_by_cols->vectors[j]->values[col_count] = res;
          count++;
          temp_by_rows->vectors[k]->indices[row_counts[k]] = j;
          temp_by_rows->vectors[k]->values[row_counts[k]] = res;
          row_counts[k]++;
        }
      }
      //save column counts
      temp_by_cols->vectors[j]->length = count;
      col_count = 0;
    } //end computation

    //save row counts
    for(int l = 0; l < n; l++) {
      temp_by_rows->vectors[l]->length = row_counts[l];
    }

    //swap pointers
    temp = matrix_by_cols;
    matrix_by_cols = temp_by_cols;  //save answer as new base matrix
    temp_by_cols = temp;  //to be overwritten

    temp = matrix_by_rows;
    matrix_by_rows = temp_by_rows;  //save answer as new base matrix
    temp_by_rows = temp;    //to be overwritten
    
  } //end powers

  destroy_matrix(temp_by_rows);
  destroy_matrix(temp_by_cols);
  destroy_matrix(matrix_by_rows);
  return(matrix_by_cols);
}
