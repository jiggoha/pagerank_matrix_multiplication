#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>
#include "mpi.h"
#include "timing.h"
#include "matrix_multiplication.h"
#include "prints.h"
#include "random_list.h"

#define DEBUG (0)
#define INITIAL_SEND_COL_TAG (1)
#define SEND_ROW_TAG (2)

#define MAX(a,b) (((a)>(b))?(a):(b))


int main (int argc, char **argv) {
  srand(12345);

  if (argc != 3) {
    fprintf(stderr, "Usage: matrix_multiplication [matrix_length] [power] [vecs_per_proc]\n");
    exit(1);
  }
  int n = atoi(argv[1]);
  int p = atoi(argv[2]) >> 1;
  int vecs_per_proc = atoi(argv[3]);

  MPI_Init(&argc, &argv);

  // number of processes
  int num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  // rank of process
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Status status;

  Matrix *matrix;
  if (rank == 0) {
    // generate and distribute 
    matrix = generate_matrix(n, (DEBUG  >= 0));
  }

  // call twice:
  Matrix *col_block = newMatrix(vecs_per_proc, n);
  Matrix *row_block = newMatrix(vecs_per_proc, n);
  Matrix *result_block = newMatrix(vecs_per_proc, n);   //results stored by cols

  // receive/distribute cols
  for (int i = 0; i < n; i++) {
    int to_rank = i / num_procs;
    MPI_Send(&matrix->vectors[i]->length, 1, MPI_INT, to_rank, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD);
    MPI_Send(matrix->vectors[i]->indices, matrix->vectors[i]->length, MPI_INT, to_rank, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD);
    MPI_Send(matrix->vectors[i]->values, matrix->vectors[i]->length, MPI_DOUBLE, to_rank, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD);
  }

  for (int i = 0; i < vecs_per_proc; i++) {
    int count;
    MPI_Recv(&count, 1, MPI_INT, 0, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(col_block->vectors[i]->indices, count, MPI_INT, 0, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(col_block->vectors[i]->values, count, MPI_DOUBLE, 0, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD, &status);
  }

  // buffers
  int *send_counts = malloc(sizeof(int) * num_procs);
  int *send_displacements = malloc(sizeof(int) * num_procs);

  int **receive_counts = malloc(sizeof(int *) * vecs_per_proc);
  for (int i = 0; i < vecs_per_proc; i++) {
    receive_counts[i] = malloc(sizeof(int) * num_procs);
  }
  int *receive_displacements = malloc(sizeof(int) * num_procs);

  // iterate through power number of iterations
  for (int it = 0; it < p; it++) {
    // distribute rows

    //TODO: LENGTH COMPUTATION
    for (int i = 0; i < vecs_per_proc; i++) {
      // send_counts
      get_counts(col_block->vectors[i]->indices, col_block->vectors[i]->length, send_counts, n, num_procs);
      get_displacements(send_counts, send_displacements, num_procs);

      // send/receive counts for number of elements
      MPI_Alltoall(send_counts, 1, MPI_INT, receive_counts[i], 1, MPI_INT, MPI_COMM_WORLD);
      get_displacements(receive_counts[i], receive_displacements, num_procs);

      // send and receive indices
      MPI_Alltoallv(col_block->vectors[i]->indices, send_counts, send_displacements, MPI_INT, row_block->vectors[i]->indices, receive_counts[i], receive_displacements, MPI_INT, MPI_COMM_WORLD);

      // send and receive values
      MPI_Alltoallv(col_block->vectors[i]->values, send_counts, send_displacements, MPI_DOUBLE, row_block->vectors[i]->values, receive_counts[i], receive_displacements, MPI_DOUBLE, MPI_COMM_WORLD);
    }

    // sort rows
    for (int i = 0; i < vecs_per_proc; i++) {
      jay_sort(row_block->vectors[i], receive_counts[i]);
    }

    int col_count = 0; //number of elements so far in result column

    // calculate on rows and round robin
    for (int i = 0; i < num_procs; i++) {
      // calculate dot product
      for (int x = 0; x < vecs_per_proc; x++) {   // cols
        for (int y = 0; y < vecs_per_proc; y++) {    // rows
          double result = dot_product(col_block->vectors[x], row_block->vectors[y]);
          if (fabs(result) > 0.001) {
            result_block->vectors[x]->indices[col_count] = y;
            result_block->vectors[x]->values[col_count] = result;
            col_count++;
          }
        } // end rows
        result_block->vectors[x]->length = col_count;
        col_count = 0;
      } // end cols

      // swap rows  
      // TODO: possibly use sendrecv?
      for (int j = 0; j < vecs_per_proc; j++) {
        int next_rank = next_rank(rank, num_procs);
        int prev_rank = prev_rank(rank, num_procs);

        // send count
        int send_count = row_block->vectors[j]->length;
        MPI_Send(&send_count, 1, MPI_INT, next_rank, SEND_ROW_TAG, MPI_COMM_WORLD);

        // send actual row
        MPI_Send(row_block->vectors[j]->indices, send_count, MPI_INT, next_rank, SEND_ROW_TAG, MPI_COMM_WORLD);
        MPI_Send(row_block->vectors[j]->values, send_count, MPI_DOUBLE, next_rank, SEND_ROW_TAG, MPI_COMM_WORLD);

        // receive count
        int recv_count;
        MPI_Recv(&recv_count, 1, MPI_INT, prev_rank, SEND_ROW_TAG, MPI_COMM_WORLD, &status);

        // receive actual row
        MPI_Recv(row_block->vectors[j]->indices, recv_count, MPI_INT, prev_rank, SEND_ROW_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(row_block->vectors[j]->values, recv_count, MPI_DOUBLE, prev_rank, SEND_ROW_TAG, MPI_COMM_WORLD, &status);
      }
    }
  }

  //TODO: GATHER RESULTS

  MPI_Finalize();

  print_matrix(matrix);
  Matrix *serialResult = serial(matrix, p);
  print_matrix(serialResult);

  printf("Are the matrices the same?\n");
  if (are_matrices_same(serialResult, parallelResult)) {
    printf("Yes!\n");
  } else {
    printf("No :(\n");
  }
  return 0;
}

int next_rank(int rank, int num_procs) {
  if (rank == num_procs - 1) {
    return 0;
  } else {
    return rank + 1;
  }
}

int prev_rank(int rank, int num_procs) {
  if (rank == 0) {
    return num_procs - 1;
  } else {
    return rank - 1;
  }
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


Matrix* newMatrix(int n, int m) {
  Matrix *new = malloc(sizeof(Matrix));
  new->n = n;
  new->vectors = malloc(sizeof(Vector *) * n);

  // memory stuff
  for (int i = 0; i < n; i++) {
    // actually doesn't need to be this long, but it's an upper bound
    new->vectors[i] = malloc(sizeof(Vector));
    new->vectors[i]->length = 0;
    new->vectors[i]->indices = malloc(sizeof(int) * m);
    new->vectors[i]->values = malloc(sizeof(double) * m);
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
            || fabs(a->vectors[i]->values[j] - b->vectors[i]->values[j]) > 0.001) {
        printf("different at element %d, %d by %f\n", i, j, a->vectors[i]->values[j] - b->vectors[i]->values[j]);
        return 0;
      }
    }
  }

  return 1;
}


Matrix *transpose_representation(Matrix *matrix) {
  Matrix* transposed = newMatrix(matrix->n, matrix->n);

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

double dot_product(Vector *col, Vector *row) {
  double result = 0.;
  int index_col = 0;
  int index_row = 0;

  while (index_col < col->length  && index_row < row->length) {
    if (col->indices[index_col] == row->indices[index_row]) {
      result += col->values[index_col] * row->values[index_row];
      index_col++;
      index_row++;
    } else if (col->indices[index_col] > row->indices[index_row]) {
      index_row++;
    } else {
      index_col++;
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

  int n = matrix_by_cols->n;
  Matrix* matrix_by_rows = transpose_representation(matrix_by_cols);
  Matrix* temp_by_cols = newMatrix(n, n);
  Matrix* temp_by_rows = newMatrix(n, n);
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
        if(fabs(res = dot_product(matrix_by_cols->vectors[j], matrix_by_rows->vectors[k])) > 0.001) {
          temp_by_cols->vectors[j]->indices[col_count] = k;
          temp_by_cols->vectors[j]->values[col_count] = res;
          col_count++;
          temp_by_rows->vectors[k]->indices[row_counts[k]] = j;
          temp_by_rows->vectors[k]->values[row_counts[k]] = res;
          row_counts[k]++;
        }
      }
      //save column counts
      temp_by_cols->vectors[j]->length = col_count;
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
