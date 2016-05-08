#include <stdlib.h> 

static int rand_in_range(int min, int max) {
  return rand() % (max - min) + min;
}

// swaps first k numbers
static void shuffle(int *arr, int length, int k) {
  for (int i = 0; i < k; i++) {
    int to_swap_with = rand_in_range(i, length);

    // swap
    int temp = arr[i];
    arr[i] = arr[to_swap_with];
    arr[to_swap_with] = temp;
  }
}

static int *amelia_sort(int *arr, int length, int k) {
  int *ones_and_zeroes = calloc(length, sizeof(int));

  for (int i = 0; i < k; i++) {
    ones_and_zeroes[arr[i]] = 1;
  }

  int *sorted = malloc(sizeof(int) * length);
  int count = 0;
  for (int i = 0; i < length; i++) {
    if (ones_and_zeroes[i] == 1) {
      sorted[count] = i;
      count++;
    }
  }

  free(ones_and_zeroes);
  return sorted;
}

int *random_increasing_ints(int max, int k) {
  int *full_list = malloc(sizeof(int) * max);

  for (int i = 0; i < max; i++) {
    full_list[i] = i;
  }

  shuffle(full_list, max, k);
  int *sorted = amelia_sort(full_list, max, k);

  free(full_list);
  return sorted;
}