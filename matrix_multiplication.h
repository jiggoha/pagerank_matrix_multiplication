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
double dot_product(Vector *col, Vector *row);
void get_counts(int *indices, int length, int *send_counts, int n, int num_procs); // output: send_counts
void get_displacements(int *send_counts, int *displacements, int num_procs); // output: displacements
int next_rank(int rank, int num_procs);
int prev_rank(int rank, int num_procs);
void jay_sort(Vector* row);

// initializations
int *random_increasing_ints(int max, int k);
Vector *generate_vector(int n, int length, int debug);
void destroy_vector(Vector *vec);
Matrix *newMatrix(int n, int m);
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
