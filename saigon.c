#include "saigon_base.c"

int main(int argc, char *argv[]) {
  init_static_vars();
  start = clock();

  int start_k = 2; // default value
  if (argc > 1) {
    start_k = atoi(argv[1]);
    if (start_k < 2 || start_k >= MAX_N) {
      printf("Starting k must be between 2 and %d\n", MAX_N - 1);
      return 1;
    }
  }

  for (int k = start_k; k < MAX_N; k++) {
    gen_all(k);
  }

  return 0;
}
