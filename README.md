#PageRank Matrix Multiplication
Parallelizing raising sparse matrices to high powers with OpenMPI

CPSC 424 Parallel Programming Final Project by Amelia Holcomb and Jay Hou

### Smart things that we did
- dot product
- log2 of power instead of power # of iterations
- To change columns to rows between processors, we wanted to use scatterv. But then we realized that this isn't so nice because if we assign more than one column per process, we don't know where to put the entries in the receiving process. Instead, we'll receive them out of order. To ameliorate this, we are taking it out of order, but doing a fast n*log(procs) merge-sort on the vectors.
