#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <stdint.h>
#include <assert.h>

#define MAX_N 30

#define base_row (row + 1)
#define fwd_row_count (size - base_row)

// Static storage for temporary variables (up to MAX_NxMAX_N matrices)
static int64_t Mi[MAX_N][MAX_N];
static int64_t rhsi[MAX_N];
static mpq_t temp;
static mpq_t V[MAX_N];
static mpq_t lower_V[MAX_N];
static mpq_t upper_V[MAX_N];
static mpq_t upper_r, lower_r;

static int G[MAX_N][MAX_N] = {0};
static int lower_G[MAX_N][MAX_N];
static int upper_G[MAX_N][MAX_N];
static int lower_size;
static int upper_size;
static clock_t start;


void init_static_vars() {
  mpq_inits(temp, upper_r, lower_r, NULL);
  for (int i = 0; i < MAX_N; i++) {
    mpq_init(V[i]);
    mpq_init(lower_V[i]);
    mpq_init(upper_V[i]);
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

typedef struct {
  int from;
  int to;
} Edge;

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


// ---------- Example usage ----------

void calc(int n, char *edge_list) {
  populate_conductance_matrix(n, edge_list);
  calculate_voltage_vector_from_conductance(n);
  gmp_printf("Equivalent resistance = %Qd for %s\n", V[0], edge_list);
}

void set_g(int i, int j, int g) {
  G[i][j] = g;
  G[j][i] = g;
}

void gen_m(int size, int k, int row, int max_sink_links);

void gen_m2(int size, int k, int row, const int groups[MAX_N], int group_index, int group_count, int max_sink_links) {
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

  int max_to_use = k - (size - row - 2);
  assert(max_to_use > 0 && "max_to_use must be positive");

  int current_sum = 0;
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
    int max_sink_link_next = group_index == 0 ? max_sink_links - current_sum : max_sink_links;
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
        int to_return = G[row][group_indices[i]];
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


void gen_m(int size, int k, int row, int max_sink_links) {
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

    calculate_voltage_vector_from_conductance(size);

    int cmp = mpq_cmp_si(V[0], 1, 1);
    if (cmp < 0) {
      int cmp1 = mpq_cmp(V[0], lower_r);
      if (cmp1 > 0) {
        mpq_set(lower_r, V[0]);
        memcpy(lower_G, G, sizeof(G));
        lower_size = size;
        for (int i = 0; i < size; i++) {
          mpq_set(lower_V[i], V[i]);
        }
      }
    } else if (cmp > 0) {
      int cmp1 = mpq_cmp(V[0], upper_r);
      if (cmp1 < 0) {
        mpq_set(upper_r, V[0]);
        memcpy(upper_G, G, sizeof(G));
        upper_size = size;
        for (int i = 0; i < size; i++) {
          mpq_set(upper_V[i], V[i]);
        }
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

void gen_all(int k) {
  mpq_set_si(lower_r, 0, 1);
  mpq_set_si(upper_r, k + 1, 1);
  for (int size = 2; size <= k + 1; size++) {
    memset(G, 0, sizeof(G));
    gen_m(size, k, 0, k);
  }
  gmp_printf("for k = %d lr = %Qd ur = %Qd t=%ds\n", k, lower_r, upper_r, (long) (clock() - start) / CLOCKS_PER_SEC);
  printf("Lower G matrix (size %d):\n", lower_size);
  for (int i = 0; i < lower_size; i++) {
    for (int j = 0; j < lower_size; j++) {
      printf("%d ", lower_G[i][j]);
    }
    printf("\n");
  }
  printf("Lower voltage vector:\n");
  for (int i = 0; i < lower_size - 1; i++) {
    gmp_printf("%Qd ", lower_V[i]);
  }
  printf("\nUpper G matrix (size %d):\n", upper_size);
  for (int i = 0; i < upper_size; i++) {
    for (int j = 0; j < upper_size; j++) {
      printf("%d ", upper_G[i][j]);
    }
    printf("\n");
  }

  printf("Upper voltage vector:\n");
  for (int i = 0; i < upper_size - 1; i++) {
    gmp_printf("%Qd ", upper_V[i]);
  }
  printf("\n");

  fflush(stdout);
}


int main() {
  init_static_vars();
  start = clock();

  // 0 1 1 0 0 0
  // 1 0 0 1 1 0
  // 1 0 0 1 0 1

  // set_g(0, 1, 1);
  // set_g(0, 2, 1);
  // set_g(1, 3, 1);
  // set_g(1, 4, 1);
  // set_g(2, 3, 1);
  // set_g(2, 5, 1);


  // gen_m(4, 4, 0, 4);

  // gen_all(4);

  for (int k = 1; k < MAX_N; k++) {
    gen_all(k);
  }


  //
  // calc(2, "[1,2]^2");
  // calc(3, "[1,2],[2,3]");
  // calc(3, "[1,2],[1,3],[2,3]");
  // calc(3, "[1,2],[2,3]^2");
  // calc(4, "[1,2],[1,4],[2,3],[3,4]");
  // calc(3, "[1,2],[2,3]^3");
  // calc(4, "[1,2],[1,3],[2,4],[3,4]^2");
  // calc(4, "[1,2],[1,3],[2,3],[3,4]^2");
  // calc(5, "[1,2],[1,3],[2,5],[3,4],[3,5],[4,5]");
  // calc(4, "[1,2],[1,3],[2,3]^2,[3,4]^2");
  // calc(5, "[1,2],[1,3],[2,3],[2,4],[3,5]^2,[4,5]");
  // calc(5, "[1,2],[1,3],[2,4],[3,4],[3,5],[4,5]^2");
  // calc(6, "[1,2],[1,3],[2,4],[2,5],[3,4],[3,6],[4,6],[5,6]");
  // calc(5, "[1,2],[1,3],[2,3],[2,4],[3,4],[3,5],[4,5]^2");
  // calc(6, "[1,2],[1,3],[2,4]^2,[3,5],[3,6],[4,5],[4,6],[5,6]");
  // calc(6, "[1,2],[1,3],[2,3],[2,4],[3,4],[3,5],[3,6],[4,6],[5,6]");
  // calc(7, "[1,2],[1,3],[2,4],[2,5],[3,6],[3,7],[4,5],[4,6],[5,7],[6,7]");
  // calc(6, "[1,2],[1,3],[2,3],[2,4],[3,4],[3,5],[3,6],[4,5],[4,6],[5,6]");
  // calc(7, "[1,2],[1,3],[2,3],[2,4],[3,5],[3,7],[4,6],[4,7],[5,6],[5,7],[6,7]");
  // calc(7, "[1,2],[1,3],[2,4],[2,5],[3,5],[3,7],[4,5],[4,6],[5,6],[5,7],[6,7]");
  // calc(7, "[1,2],[1,3],[2,3],[2,4],[2,5],[3,4],[3,5],[3,7],[4,6],[5,6],[5,7],[6,7]");
  // calc(7, "[1,2],[1,3],[2,4],[2,5],[3,4],[3,5],[3,7],[4,5],[4,6],[5,6],[5,7],[6,7]");
  // calc(8, "[1,2],[1,3],[2,4],[2,5],[2,6],[3,4],[3,5],[3,8],[4,6],[5,7],[6,7],[6,8],[7,8]");
  // calc(8, "[1,2],[1,3],[2,3],[2,4],[2,6],[3,5],[3,8],[4,5],[4,7],[5,6],[5,7],[6,8],[7,8]");
  // calc(8, "[1,2],[1,3],[2,3],[2,4],[2,5],[3,4],[3,6],[3,8],[4,6],[4,7],[5,6],[5,7],[6,8],[7,8]");
  // calc(9, "[1,2],[1,3],[2,4],[2,5],[2,7],[3,4],[3,5],[3,9],[4,6],[5,8],[6,7],[6,8],[7,9],[8,9]");
  // calc(9,
  //      "[1,2],[1,3],[2,4],[2,7],[3,4],[3,5],[3,9],[4,6],[4,8],[5,7],[5,8],[6,7],[6,8],[7,9],[8,9]");
  // calc(9,
  //      "[1,2],[1,3],[2,4],[2,5],[2,7],[3,4],[3,5],[3,9],[4,5],[4,6],[5,8],[6,7],[6,8],[7,9],[8,9]");
  // calc(9,
  //      "[1,2],[1,3],[2,4],[2,5],[2,7],[3,4],[3,5],[3,9],[4,5],[4,6],[5,8],[6,7],[6,8],[7,8],[7,9],[8,9]");
  // calc(10,
  //      "[1,2],[1,3],[2,4],[2,5],[2,7],[3,4],[3,6],[3,10],[4,8],[5,6],[5,9],[6,8],[7,8],[7,9],[8,10],[9,10]");
  // calc(10,
  //      "[1,2],[1,3],[2,4],[2,5],[2,7],[3,4],[3,6],[3,10],[4,8],[5,6],[5,7],[5,9],[6,8],[7,8],[7,9],[8,10],[9,10]");
  // calc(10,
  //      "[1,2],[1,3],[2,4],[2,5],[2,7],[3,4],[3,6],[3,10],[4,6],[4,8],[5,6],[5,8],[6,9],[7,8],[7,9],[8,10],[9,10]");
  // calc(10,
  //      "[1,2],[1,3],[2,4],[2,5],[2,6],[3,4],[3,5],[3,10],[4,8],[5,7],[5,9],[6,7],[6,8],[6,9],[7,8],[7,9],[8,10],[9,10]");
  // calc(11,
  //      "[1,2],[1,3],[2,4],[2,6],[2,7],[3,5],[3,8],[3,11],[4,5],[4,9],[5,7],[6,8],[6,9],[7,9],[7,10],[8,10],[9,11],[10,11]");
  // calc(11,
  //      "[1,2],[1,3],[2,5],[2,6],[2,7],[3,4],[3,5],[3,11],[4,6],[4,7],[5,8],[5,10],[6,8],[6,9],[7,9],[7,10],[8,9],[9,11],[10,11]");
  // calc(11,
  //      "[1,2],[1,3],[2,4],[2,5],[2,7],[3,5],[3,6],[3,11],[4,6],[4,9],[5,8],[5,10],[6,7],[6,8],[7,9],[7,10],[8,9],[9,11],[10,11]");

  return 0;
}
