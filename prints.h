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