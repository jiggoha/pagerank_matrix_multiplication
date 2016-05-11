#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>
#include "mpi.h"
#include "timing.h"
#include "matrix_multiplication.h"
#include "prints.h"T

#define DEBUG (-1)
#define INITIAL_SEND_COL_TAG (1)
#define SEND_ROW_TAG (2)

#define MAX(a,b) (((a)>(b))?(a):(b))


int main (int argc, char **argv) {
  srand(12345);
  double serial_start_time, serial_end_time, cputime;
  double parallel_start_time, parallel_end_time;

  if (argc != 3) {
    fprintf(stderr, "Usage: matrix_multiplication [matrix_length] [power]\n");
    exit(1);
  }
  int n = atoi(argv[1]);
  int p = atoi(argv[2]);

  MPI_Init(&argc, &argv);

  // number of processes
  int num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  int vecs_per_proc = n / num_procs;
  
  // rank of process
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Status status;
  MPI_Request sendrequest[vecs_per_proc], recvrequest[vecs_per_proc];

  Matrix *matrix;
  Matrix *serial_result;
  Matrix *parallel_result;
  if (rank == 0) {
    // generate and distribute 
    matrix = generate_matrix(n, (DEBUG  >= 0));
    // print_matrix(matrix);
    serial_result = newMatrix(n, n);
    parallel_result = newMatrix(n, n);
  }

  // call twice:
  Matrix *col_block = newMatrix(vecs_per_proc, n);
  Matrix *row_block = newMatrix(vecs_per_proc, n);
  Matrix *nextrow_block = newMatrix(vecs_per_proc, n);
  Matrix *result_block = newMatrix(vecs_per_proc, n);   //results stored by cols

  // distribute columns
  if(rank == 0) {
    timing(&parallel_start_time, &cputime);
    for (int i = 0; i < n; i++) {
      int to_rank = i / vecs_per_proc;
      MPI_Send(&matrix->vectors[i]->length, 1, MPI_INT, to_rank, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD);
      MPI_Send(matrix->vectors[i]->indices, matrix->vectors[i]->length, MPI_INT, to_rank, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD);
      MPI_Send(matrix->vectors[i]->values, matrix->vectors[i]->length, MPI_DOUBLE, to_rank, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD);
    }
  }

  for (int i = 0; i < vecs_per_proc; i++) {
    int count;
    MPI_Recv(&count, 1, MPI_INT, 0, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD, &status);
    col_block->vectors[i]->length = count;
    MPI_Recv(col_block->vectors[i]->indices, count, MPI_INT, 0, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(col_block->vectors[i]->values, count, MPI_DOUBLE, 0, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD, &status);
  }

  // buffers
  int *send_counts = malloc(sizeof(int) * num_procs);
  int *send_displacements = malloc(sizeof(int) * num_procs);

  int *receive_counts = malloc(sizeof(int) * num_procs);
  int *receive_displacements = malloc(sizeof(int) * num_procs);
  int *receive_idx_buf = malloc(sizeof(int) * n);
  double *receive_val_buf = malloc(sizeof(double) * n);

  int new_row_idx;
  int length;

  MPI_Barrier(MPI_COMM_WORLD);

  // iterate through power number of iterations
  for(int it = 0; it < p; it++) {
    // distribute rows
    for(int i = 0; i < vecs_per_proc; i++) {
      // send_counts
      get_counts(col_block->vectors[i]->indices, col_block->vectors[i]->length, send_counts, n, num_procs);
      get_displacements(send_counts, send_displacements, num_procs);

      // send/receive counts for number of elements
      MPI_Alltoall(send_counts, 1, MPI_INT, receive_counts, 1, MPI_INT, MPI_COMM_WORLD);
      get_displacements(receive_counts, receive_displacements, num_procs);

      // send and receive indices
      MPI_Alltoallv(col_block->vectors[i]->indices, send_counts, send_displacements, MPI_INT, receive_idx_buf, receive_counts, receive_displacements, MPI_INT, MPI_COMM_WORLD);

      // send and receive values
      MPI_Alltoallv(col_block->vectors[i]->values, send_counts, send_displacements, MPI_DOUBLE, receive_val_buf, receive_counts, receive_displacements, MPI_DOUBLE, MPI_COMM_WORLD);

      // copy indices and values correctly into rows:
      for (int j = 0; j < num_procs; j++) {    //process number
        for (int k = 0; k < receive_counts[j]; k++) {  //number of elements received from process j

          // index to location in the receive_idx_buf and receive_val_buf buffers
          int buf_i = receive_displacements[j] + k;

          // index of row in which to store, subtracting so that row indices start at 0 for each proc
          new_row_idx = receive_idx_buf[buf_i] - vecs_per_proc * rank;   
          length = row_block->vectors[new_row_idx]->length++;             //length of row thus far

          row_block->vectors[new_row_idx]->indices[length] = j * vecs_per_proc + i;   // store index
          row_block->vectors[new_row_idx]->values[length] = receive_val_buf[buf_i]; //store value
        }
      }
    } // end distribute rows

    // sort rows
    for (int i = 0; i < vecs_per_proc; i++) {
      jay_sort(row_block->vectors[i]);
    }

    // calculate on rows and round robin
    for (int i = 0; i < num_procs; i++) {

      // swap rows  
      for (int j = 0; j < vecs_per_proc; j++) {
        int next = prev_rank(rank, num_procs);
        int prev = next_rank(rank, num_procs);

        // send count
        int send_count = row_block->vectors[j]->length;
        MPI_Send(&send_count, 1, MPI_INT, next, SEND_ROW_TAG, MPI_COMM_WORLD);

        // send actual row
        MPI_Isend(row_block->vectors[j]->indices, send_count, MPI_INT, next, SEND_ROW_TAG, MPI_COMM_WORLD, &sendrequest[j]);
        MPI_Isend(row_block->vectors[j]->values, send_count, MPI_DOUBLE, next, SEND_ROW_TAG, MPI_COMM_WORLD, &sendrequest[j]);

        // receive count
        int recv_count;
        MPI_Recv(&recv_count, 1, MPI_INT, prev, SEND_ROW_TAG, MPI_COMM_WORLD, &status);
        nextrow_block->vectors[j]->length = recv_count;

        // receive actual row
        MPI_Irecv(nextrow_block->vectors[j]->indices, recv_count, MPI_INT, prev, SEND_ROW_TAG, MPI_COMM_WORLD, &recvrequest[j]);
        MPI_Irecv(nextrow_block->vectors[j]->values, recv_count, MPI_DOUBLE, prev, SEND_ROW_TAG, MPI_COMM_WORLD, &recvrequest[j]);
      }

      // calculate dot product
      for (int x = 0; x < vecs_per_proc; x++) {   // cols
        for (int y = 0; y < vecs_per_proc; y++) {    // rows
          double result = dot_product(col_block->vectors[x], row_block->vectors[y]);
          if (fabs(result) > 0.001) {
            int length = result_block->vectors[x]->length;
            result_block->vectors[x]->indices[length] = (y + rank * vecs_per_proc + vecs_per_proc * i) % n;
            result_block->vectors[x]->values[length] = result;
            result_block->vectors[x]->length++;
          }
        } // end rows
      } // end cols

      // wait on send and receive
      for (int j = 0; j < vecs_per_proc; j++) {
        MPI_Wait(&sendrequest[j], &status);
        MPI_Wait(&recvrequest[j], &status);
      }

      // swap out next row block and current row block
      Matrix* temp = nextrow_block;
      nextrow_block = row_block;
      row_block = temp;

    } // end round robin

    for(int w = 0; w < vecs_per_proc; w++) {
      jay_sort(result_block->vectors[w]); 
    }

    // swap pointers to prepare for next iteration
    Matrix *temp = col_block;
    col_block = result_block;
    result_block = temp;

    // wipe clean the result block
    for (int w = 0; w < vecs_per_proc; w++) {
      result_block->vectors[w]->length = 0;
      row_block->vectors[w]->length = 0;
    }
  } // end iterative computation (powers)

  // gather results
  for (int i = 0; i < vecs_per_proc; i++) {
    MPI_Send(&col_block->vectors[i]->length, 1, MPI_INT, 0, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD);
    MPI_Send(col_block->vectors[i]->indices, col_block->vectors[i]->length, MPI_INT, 0, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD);
    MPI_Send(col_block->vectors[i]->values, col_block->vectors[i]->length, MPI_DOUBLE, 0, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD);
  }


  if (rank == 0) {
    for (int i = 0; i < n; i++) {
      int count;
      int from_rank = i / vecs_per_proc;
      MPI_Recv(&count, 1, MPI_INT, from_rank, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD, &status);
      
      parallel_result->vectors[i]->length = count;
      MPI_Recv(parallel_result->vectors[i]->indices, count, MPI_INT, from_rank, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(parallel_result->vectors[i]->values, count, MPI_DOUBLE, from_rank, INITIAL_SEND_COL_TAG, MPI_COMM_WORLD, &status);
    }
    timing(&parallel_end_time, &cputime);

    timing(&serial_start_time, &cputime);
    Matrix *serial_result = serial(matrix, p);
    timing(&serial_end_time, &cputime);

    printf("Are the matrices the same?\n");
    if (are_matrices_same(serial_result, parallel_result)) {
      printf("Yes!\n");
    } else {
      printf("No :(\n");
    }

    printf("Serial time: %f\n", serial_end_time - serial_start_time);
    printf("Parallel time: %f\n", parallel_end_time - parallel_start_time);
  }

  MPI_Finalize();
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
          vec->values[i] = i + 1;
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

void jay_sort(Vector* vec) {

  int n = vec->length;
  int d, temp_idx;
  double temp_val;

  //actually insertion sort right now
  for (int i = 1 ; i <= n - 1; i++) {
    d = i;
 
    while (d > 0 && vec->indices[d] < vec->indices[d-1]) {
      temp_idx = vec->indices[d];
      temp_val = vec->values[d];

      vec->indices[d] = vec->indices[d-1];
      vec->values[d] = vec->values[d-1];
      
      vec->indices[d-1] = temp_idx;
      vec->values[d-1] = temp_val;
 
      d--;
    }
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
