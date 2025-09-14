#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <stdint.h>
#include <assert.h>

#define MAX_N 23

#define base_row (row + 1)
#define fwd_row_count (size - base_row)

// Static storage for temporary variables (up to MAX_NxMAX_N matrices)
static int64_t Mi[MAX_N][MAX_N];
static int64_t rhsi[MAX_N];
static mpq_t temp;
static mpq_t V[MAX_N];
static mpq_t upper_r, lower_r;
static mpq_t R;

static int64_t G[MAX_N][MAX_N] = {0};
static int64_t lower_G[MAX_N][MAX_N];
static int64_t upper_G[MAX_N][MAX_N];
static int lower_size;
static int upper_size;
static clock_t start;


void init_static_vars() {
  mpq_inits(temp, upper_r, lower_r, R, NULL);
  for (int i = 0; i < MAX_N; i++) {
    mpq_init(V[i]);
  }
}

// ---------- Bareiss elimination (fraction-free on integers) ----------

/*
 * Integer Bareiss elimination followed by mpq_t back-substitution.
 * Uses global Mi (n x n) and rhsi (n), writes results to V (mpq_t).
 * Mi and rhsi are modified in-place during elimination.
 */
void bareiss_solve(int n) {
  // Elimination (fraction-free on int64_t)
  int64_t prev_pivot = 1;
  for (int k = 0; k < n - 1; k++) {
    int64_t pivot = Mi[k][k];
    assert(pivot != 0);
    for (int i = k + 1; i < n; i++) {
      for (int j = k + 1; j < n; j++) {
        __int128 val = (__int128) Mi[i][j] * pivot - (__int128) Mi[i][k] * Mi[k][j];
        Mi[i][j] = (int64_t) (val / prev_pivot);
      }
      // RHS update
      __int128 val_rhs = (__int128) rhsi[i] * pivot - (__int128) Mi[i][k] * rhsi[k];
      rhsi[i] = (int64_t) (val_rhs / prev_pivot);
    }
    prev_pivot = pivot;
  }

  // Back substitution with mpq_t
  for (int i = n - 1; i >= 0; i--) {
    mpq_set_si(V[i], rhsi[i], 1);
    for (int j = i + 1; j < n; j++) {
      mpq_set_si(temp, Mi[i][j], 1);
      mpq_mul(temp, temp, V[j]); // temp = Mi[i][j] * V[j]
      mpq_sub(V[i], V[i], temp);
    }
    mpq_set_si(temp, Mi[i][i], 1);
    mpq_div(V[i], V[i], temp); // V[i] /= Mi[i][i]
  }
}

/*
 * Computes node voltages for an integer conductance matrix C (NxN).
 * - sink is node N-1 (grounded to 0V)
 * - source is node 0 with 1A injection
 * Writes the result into preallocated, initialized vector V of length N (V[i] must be mpq_init'ed by the caller).
 * All temporaries are allocated on the stack.
 */
void calculate_voltage_vector_from_conductance(const int N) {
  if (N < 2) return;

  int n = N - 1;

  // Reduced conductance (Laplacian) matrix Lprime directly from C (exclude sink row/col)
  for (int i = 0; i < n; i++) {
    // Diagonal: sum of all conductances connected to i
    long long diag_sum = 0;
    for (int k = 0; k < N; k++) {
      if (k == i) continue;
      diag_sum += G[i][k];
    }
    Mi[i][i] = diag_sum;

    // Off-diagonals: -conductance between originalI and originalJ
    for (int j = 0; j < n; j++) {
      if (i == j) continue;
      Mi[i][j] = -G[i][j];
    }
  }

  // Current injection vector b (size n): +1A at sourceIndex (source=0 -> sourceIndex=0)
  for (int i = 0; i < n; i++)
    rhsi[i] = 0;
  rhsi[0] = 1;

  // Solve L' v = b
  bareiss_solve(n);
}

void populate_conductance_matrix(int n, const char *edge_list) {
  // Clear matrix
  memset(G, 0, sizeof(G));

  const char *p = edge_list;
  while (*p) {
    // Parse [from,to]^conductance
    int from, to, conductance = 1;
    if (sscanf(p, "[%d,%d]^%d", &from, &to, &conductance) >= 2) {
      // Convert 1-based to 0-based indices
      from--;
      to--;
      // printf("%d -> %d = %d\n", from, to, conductance);

      // Set conductances
      if (from >= 0 && from < n && to >= 0 && to < n) {
        G[from][to] = conductance;
        G[to][from] = conductance; // Symmetric matrix
      }
    }

    // Move to next edge
    p++;
    while (*p && *p != '[') {
      p++;
    }
  }
}

/*
 * Fraction-free determinant (Bareiss) on int64_t matrix (in-place).
 * Uses __int128 intermediates to avoid overflow in the product/subtraction step.
 * Returns the exact determinant as int64_t (beware potential overflow for larger values).
 */
static int64_t bareiss_det_int64(int n, int64_t A[MAX_N][MAX_N]) {
  int64_t prev_pivot = 1;
  for (int k = 0; k < n - 1; k++) {
    int64_t pivot = A[k][k];
    assert(pivot != 0);
    for (int i = k + 1; i < n; i++) {
      for (int j = k + 1; j < n; j++) {
        __int128 val = (__int128)A[i][j] * pivot - (__int128)A[i][k] * A[k][j];
        A[i][j] = (int64_t)(val / prev_pivot);
      }
    }
    prev_pivot = pivot;
  }
  return A[n - 1][n - 1];
}

/*
 * Compute the equivalent resistance Geq between node 0 (source) and node N-1 (sink)
 * for the current conductance matrix G (NxN, symmetric, integer).
 * Writes Geq as an exact rational into out_geq (caller must mpq_init it).
 *
 * Geq = det(L') / det(L' without row/col 0), where L' is the reduced Laplacian
 * (L with sink row/col removed). This avoids solving the full system.
 */
void compute_equivalent_resistance_int64(int N) {
  const int n = N - 1;

  // Build reduced Laplacian L' into A (size n x n)
  int64_t A[MAX_N][MAX_N];
  for (int i = 0; i < n; i++) {
    long long diag_sum = 0;
    for (int k = 0; k < N; k++) {
      if (k == i) continue;
      diag_sum += G[i][k];
    }
    A[i][i] = (int64_t)diag_sum;
    for (int j = 0; j < n; j++) {
      if (i == j) continue;
      A[i][j] = - (int64_t)G[i][j];
    }
  }

  // Copy A to compute D = det(A)
  int64_t A_det[MAX_N][MAX_N];
  memcpy(A_det, A, sizeof(A_det));


  int64_t D = bareiss_det_int64(n, A_det);

  // Build A_minor = A with row/col 0 removed (size (n-1) x (n-1))
  int64_t D0 = 1;
  if (n>1) {
    int64_t A_minor[MAX_N][MAX_N];
    for (int i = 1; i < n; i++) {
      memcpy(&A_minor[i-1][0], &A[i][1], (n-1) * sizeof(int64_t));
    }
    D0 = bareiss_det_int64(n - 1, A_minor);
  }
  mpq_set_si(R, D0, D);
  mpq_canonicalize(R);
}


// ---------- Example usage ----------

void calc(int n, char *edge_list) {
  populate_conductance_matrix(n, edge_list);
  calculate_voltage_vector_from_conductance(n);
  gmp_printf("Equivalent resistance = %Qd for %s\n", V[0], edge_list);
}

void set_g(int i, int j, int64_t g) {
  G[i][j] = g;
  G[j][i] = g;
}

void gen_m(int size, int64_t k, int row, int64_t max_sink_links);

void gen_m2(int size, int64_t k, int row, const int groups[MAX_N], int group_index, int group_count, int64_t max_sink_links) {
  if (k <= 0) return;
  // find indices for the group
  int group_indices[MAX_N];
  int group_size = 0;
  for (int i = 0; i < fwd_row_count; i++) {
    if (groups[i] == group_index) {
      group_indices[group_size++] = base_row + i;
    }
  }

  //reset the group
  for (int j = 0; j < group_size; j++) {
    assert(G[row][group_indices[j]] == 0 && "group must be empty");
  }

  int64_t max_to_use = k - (size - row - 2);
  assert(max_to_use > 0 && "max_to_use must be positive");

  int64_t current_sum = 0;
  if (group_index == group_count - 1) {
    // last group
    int used_any = 0;
    for (int i = row + 1; i < size; i++) {
      if (G[row][i] != 0) {
        used_any = 1;
        break;
      }
    }
    if (!used_any) {
      // if we have not used any links in this row yet, then
      // the last group should have at least one non-zero element
      set_g(row, group_indices[0], 1);
      current_sum = 1;
    }
  }

  // iterate over all unique partitions (order unimportant) using between min_to_use and max_to_use among the group)
  int i = 0;
  int partition_count = 0;
  while (1) {
    // ----
    // if (group_size > 1) {
    //   printf("partition (%d) for group size %d max_to_use %d: ",
    //          ++partition_count,
    //          group_size,
    //          max_to_use);
    //   for (int j = 0; j < group_size; j++) {
    //     printf(" %d", G[row][group_indices[j]]);
    //   }
    //   printf("\n");
    // }
    // ----
    // continue recursively
    int64_t max_sink_link_next = group_index == 0 ? max_sink_links - current_sum : max_sink_links;
    if (group_index == group_count - 1 || current_sum == max_to_use) {
      // if this was the last group or if we used max already - go to the next row
      if (row == 0) {
        // to exclude inversions, we allow sink to have no more links than the source,
        // so we set msln to the sum of the first row excluding the last (direct links source to sink)
        max_sink_link_next = 0;
        for (int j = 1; j < size - 1; j++) {
          max_sink_link_next += G[0][j];
        }
      }
      gen_m(size, k - current_sum, row + 1, max_sink_link_next);
    } else {
      gen_m2(size, k - current_sum, row, groups, group_index + 1, group_count, max_sink_link_next);
    }

    if (current_sum == max_to_use
      || (group_index == 0 && current_sum == max_sink_links - 1)
      || (i > 0 && i == group_size - 1
        && G[row][group_indices[i]] == G[row][group_indices[i - 1]])) {
      // used up, backtrack
      if (i == 0) {
        // no more partitions
        break;
      }
      do {
        int64_t to_return = G[row][group_indices[i]];
        set_g(row, group_indices[i], 0);
        current_sum -= to_return;
        i--;
      } while (i > 0 && G[row][group_indices[i]] == G[row][group_indices[i - 1]]);
      set_g(row, group_indices[i], G[row][group_indices[i]] + 1);
      current_sum++;
    } else if (i < group_size - 1 && G[row][group_indices[i]] > 0) {
      // not the last element - move forward
      set_g(row, group_indices[++i], 1);
      current_sum++;
    } else {
      // increment current (no need to worry about backtracking?!)
      set_g(row, group_indices[i], G[row][group_indices[i]] + 1);
      current_sum++;
    }
  }
  // clean up
  for (int j = 0; j < group_size; j++) {
    set_g(row, group_indices[j], 0);
  }
}


void gen_m(int size, int64_t k, int row, int64_t max_sink_links) {
  int row_connected = (row == 0);
  int i = 0;
  while (i < row && !row_connected) {
    row_connected = (G[row][i] != 0);
    i++;
  }
  if (!row_connected) return;
  if (row == size - 2) {
    // last row before sink - connect remaining resistors to sink

    if (k > max_sink_links) return;

    set_g(row, row + 1, k);

    // ----
    // printf("Conductance matrix G for size=%d:\n", size);
    // for (int i = 0; i < size; i++) {
    //   for (int j = 0; j < size; j++) {
    //     printf("%d ", G[i][j]);
    //   }
    //   printf("\n");
    // }
    // printf("---\n");
    // ----

    compute_equivalent_resistance_int64(size);
    // calculate_voltage_vector_from_conductance(size);

    int cmp = mpq_cmp_si(R, 1, 1);
    if (cmp < 0) {
      int cmp1 = mpq_cmp(R, lower_r);
      if (cmp1 > 0) {
        mpq_set(lower_r, R);
        memcpy(lower_G, G, sizeof(G));
        lower_size = size;
      }
    } else if (cmp > 0) {
      int cmp1 = mpq_cmp(R, upper_r);
      if (cmp1 < 0) {
        mpq_set(upper_r, R);
        memcpy(upper_G, G, sizeof(G));
        upper_size = size;
      }
    }

    // clear
    set_g(row, row + 1, 0);
  } else {
    // group following rows if their connections are identical
    // rows from row+1
    int row_offset_to_group[MAX_N];
    memset(row_offset_to_group, 0, sizeof(row_offset_to_group));
    // sink is group 0
    int nxt_grp = 1;
    for (int row_offset = 0; row_offset < fwd_row_count - 1; row_offset++) {
      int max_group_seen = 0;
      for (int i = 0; row_offset_to_group[row_offset] == 0 && i < row_offset; i++) {
        if (row_offset_to_group[i] <= max_group_seen) continue;
        int diff = 0;
        for (int j = 0; diff == 0 && j < base_row; j++) {
          if (G[base_row + row_offset][j] != G[base_row + i][j]) {
            diff = 1;
          }
        }
        if (diff == 0) {
          // found match
          row_offset_to_group[row_offset] = row_offset_to_group[i];
        } else if (row_offset_to_group[i] == nxt_grp - 1) {
          // seen all groups already
          break;
        } else {
          max_group_seen = row_offset_to_group[i];
        }
      }
      // if no match found create a  new group
      if (row_offset_to_group[row_offset] == 0) {
        row_offset_to_group[row_offset] = nxt_grp++;
      }
    }

    //
    gen_m2(size, k, row, row_offset_to_group, 0, nxt_grp, max_sink_links);
  }
}

void clear_g(int size) {
  memset(G, 0, sizeof(G));
}

void print_adjacency_vector(int size, const int64_t matrix[MAX_N][MAX_N]) {
  int first = 1;
  for (int i = 0; i < size; i++) {
    for (int j = i + 1; j < size; j++) {
      if (matrix[i][j] > 0) {
        if (!first) printf(",");
        first = 0;
        printf("[%d,%d]", i + 1, j + 1);
        if (matrix[i][j] > 1) printf("^%lld", matrix[i][j]);
      }
    }
  }
}

void gen_all(int k) {
  mpq_set_si(lower_r, 0, 1);
  mpq_set_si(upper_r, k + 1, 1);
  for (int size = 2; size <= k + 1; size++) {
    memset(G, 0, sizeof(G));
    gen_m(size, k, 0, k);
  }

  gmp_printf("RESULT k=%d lower %Qd circuit: ", k, lower_r);
  print_adjacency_vector(lower_size, lower_G);
  printf(" voltages=");
  memcpy(G, lower_G, sizeof(G));
  calculate_voltage_vector_from_conductance(lower_size);
  for (int i = 0; i < lower_size - 1; i++) {
    if (i > 0) printf(","); else printf("[");
    gmp_printf("%Qd", V[i]);
  }
  printf("] t=%lus\n", (clock() - start) / CLOCKS_PER_SEC);

  gmp_printf("RESULT k=%d upper %Qd circuit: ", k, upper_r);
  print_adjacency_vector(upper_size, upper_G);
  printf(" voltages=");
  memcpy(G, upper_G, sizeof(G));
  calculate_voltage_vector_from_conductance(upper_size);
  for (int i = 0; i < upper_size - 1; i++) {
    if (i > 0) printf(","); else printf("[");
    gmp_printf("%Qd", V[i]);
  }
  printf("] t=%lus\n", (clock() - start) / CLOCKS_PER_SEC);

  fflush(stdout);
}


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
