#include <ctype.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "additional_fct.h"
#include "fips202.h"
#include "packing.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "sign.h"
#include "symmetric.h"

#include "lp_lib.h"

int main(int argc, char const *argv[]) {
  uint32_t NB_INEQ;
  int sk_index;
  float C_low, C_up;

  if (argc == 4) {
    // Case 1: Two arguments   (int NB_inequalities, float C)
    sk_index = atoi(argv[1]);
    NB_INEQ = atoi(argv[2]);
    C_up = strtof(argv[3], NULL);
    C_low = -C_up;
  } else if (argc == 5) {
    // Case 2: Three arguments (int NB_inequalities, float C_low, float C_up)
    sk_index = atoi(argv[1]);
    NB_INEQ = atoi(argv[2]);
    C_low = -strtof(argv[3], NULL);
    C_up = strtof(argv[4], NULL);
  } else {
    // Invalid number of arguments
    printf("Usage:\n");
    printf("  %s <int> <float>\n", argv[0]);
    printf("  %s <int> <float1> <float2>\n", argv[0]);
    return DATA_ERROR;
  }

  char t0_guess_filename[60];
  FILE *t0_guess_file;
  double t0_guess[K][N];
  double t0_guess_updated[K][N];

  char sign_compressed_file_name[128];
  FILE *sign_compressed_file;
  uint8_t sm_compressed[MLEN + COMP_CRYPTO_BYTES];

  uint32_t i, j;

  float C = C_up;
  int8_t ineq_line[N];
  uint32_t cpt_ineq[K], cpt_ineq1[K], cpt_ineq0[K];
  uint32_t poly_index, coeff_index, ineq_index, cpt_signs;

  double bound;
  double guess;

  polyveck r0, h;
  poly c;
  uint8_t c_seed[CTILDEBYTES];

  clock_t start, end;
  double cpu_time_used;

  char name_var[14];
  char name_lp_file[58];
  lprec *LPs[K];
  int *colno = NULL, ret = 0;
  REAL *row = NULL;

  int count = 0;
  int ch;
  unsigned char ich;
  int started = 0;

  char directory[32];
  struct stat st = {0};

  sprintf(directory, "../Guess/");
  if (stat(directory, &st) == -1) {
    mkdir(directory, 0700);
  }

  sprintf(directory, "../Guess/%.16s/", CRYPTO_ALGNAME);
  if (stat(directory, &st) == -1) {
    mkdir(directory, 0700);
  }

  sprintf(directory, "../Guess/%.16s/key%d/", CRYPTO_ALGNAME, sk_index);
  if (stat(directory, &st) == -1) {
    mkdir(directory, 0700);
  }

  sprintf(directory, "../Lps/");
  if (stat(directory, &st) == -1) {
    mkdir(directory, 0700);
  }

  sprintf(directory, "../Lps/%.16s/", CRYPTO_ALGNAME);
  if (stat(directory, &st) == -1) {
    mkdir(directory, 0700);
  }

  sprintf(directory, "../Lps/%.16s/key%d/", CRYPTO_ALGNAME, sk_index);
  if (stat(directory, &st) == -1) {
    mkdir(directory, 0700);
  }

  sprintf(sign_compressed_file_name,
          "../Signs/%.16s/key%d/PQCsignKAT_%.16s_compressed.rsp",
          CRYPTO_ALGNAME, sk_index, CRYPTO_ALGNAME);
  if ((sign_compressed_file = fopen(sign_compressed_file_name, "r")) == NULL) {
    printf("Couldn't open <%s> for read\n", sign_compressed_file_name);
    return DATA_ERROR;
  }

  sprintf(t0_guess_filename,
          "../Guess/%.16s/key%d/t0_guess_file.bin",
          CRYPTO_ALGNAME, sk_index);
  t0_guess_file = fopen(t0_guess_filename, "rb");
  if (t0_guess_file == NULL) {
    perror("Error opening file");
    return FILE_OPEN_ERROR;
  }

  fread(t0_guess, sizeof(double), K * N, t0_guess_file);
  fclose(t0_guess_file);

  for (poly_index = 0; poly_index < K; poly_index++) {
    cpt_ineq[poly_index] = 0;
  }

  /* We create models with 0 rows, N columns because we build the model row by
   * row */
  for (poly_index = 0; poly_index < K; poly_index++) {
    LPs[poly_index] = make_lp(0, N);
    if (LPs[poly_index] == NULL) {
      ret = 1;
    }
  }

  if (ret == 0) {
    for (poly_index = 0; poly_index < K; poly_index++) {
      for (coeff_index = 0; coeff_index < N; coeff_index++) {
        sprintf(name_var, "s%d", coeff_index);
        if (!set_col_name(LPs[poly_index], coeff_index + 1, name_var)) {
          ret = 1;
        }
      }
    }

    colno = (int *)malloc(N * sizeof(*colno));
    row = (REAL *)malloc(N * sizeof(*row));
    if ((colno == NULL) || (row == NULL)) {
      ret = 1;
    }
  }

  printf("C: %.2f\n", C_up);
  if (ret == 0) {
    for (poly_index = 0; poly_index < K; poly_index++) {
      for (coeff_index = 0; coeff_index < N; coeff_index++) {
        if (!set_lowbo(LPs[poly_index], coeff_index + 1,
                       get_max(t0_guess[poly_index][coeff_index] + C_low))) {
          ret = 1;
        }
        if (!set_upbo(LPs[poly_index], coeff_index + 1,
                      get_min(t0_guess[poly_index][coeff_index] + C_up))) {
          ret = 1;
        }
      }
    }
  }

  cpt_signs = 0;
  start = clock();
  while (find_min(cpt_ineq) < NB_INEQ) {
    started = 0;
    memset(sm_compressed, 0x00, MLEN + COMP_CRYPTO_BYTES);
    while ((ch = fgetc(sign_compressed_file)) != EOF) {
      if (!isxdigit(ch)) {
        if (!started) {
          if (ch == '\n')
            break;
          else
            continue;
        } else
          break;
      }
      started = 1;
      if ((ch >= '0') && (ch <= '9'))
        ich = ch - '0';
      else if ((ch >= 'A') && (ch <= 'F'))
        ich = ch - 'A' + 10;
      else if ((ch >= 'a') && (ch <= 'f'))
        ich = ch - 'a' + 10;
      else // shouldn't ever get here
        ich = 0;

      for (i = 0; i < MLEN + COMP_CRYPTO_BYTES - 1; i++)
        sm_compressed[i] =
            (sm_compressed[i] << 4) | (sm_compressed[i + 1] >> 4);
      sm_compressed[MLEN + COMP_CRYPTO_BYTES - 1] =
          (sm_compressed[MLEN + COMP_CRYPTO_BYTES - 1] << 4) | ich;
    }

    unpack_sig_compressed(c_seed, &r0, &h, sm_compressed);
    poly_challenge(&c, c_seed);

    for (poly_index = 0; poly_index < K; poly_index++) {
      for (coeff_index = 0; coeff_index < N; coeff_index++) {
        if (h.vec[poly_index].coeffs[coeff_index] == 1) {
          if (r0.vec[poly_index].coeffs[coeff_index] > 0) {
            creat_ineq(&c, coeff_index, ineq_line);
            guess = scalar_product(ineq_line, t0_guess[poly_index]);
            bound = -GAMMA2 - BETA - 1 + r0.vec[poly_index].coeffs[coeff_index];
            if (bound + guess < C * TAU) {
              if (ret == 0) {
                j = 0;
                for (ineq_index = 0; ineq_index < N; ineq_index++) {
                  if (ineq_line[ineq_index] != 0) {
                    colno[j] = ineq_index + 1;
                    row[j++] = ineq_line[ineq_index];
                  }
                }

                /* add the row to lpsolve */
                if (!add_constraintex(LPs[poly_index], j, row, colno, LE,
                                      bound)) {
                  ret = 1;
                }
              }
              cpt_ineq1[poly_index]++;
              cpt_ineq[poly_index]++;
            }
          }

          if (r0.vec[poly_index].coeffs[coeff_index] < 0) {
            creat_ineq(&c, coeff_index, ineq_line);
            guess = scalar_product(ineq_line, t0_guess[poly_index]);
            bound = -GAMMA2 - BETA - 1 - r0.vec[poly_index].coeffs[coeff_index];
            // if(bound < C*TAU){
            if (bound - guess < C * TAU) {
              if (ret == 0) {
                j = 0;
                for (ineq_index = 0; ineq_index < N; ineq_index++) {
                  if (ineq_line[ineq_index] != 0) {
                    colno[j] = ineq_index + 1;
                    row[j++] = -ineq_line[ineq_index];
                  }
                }

                /* add the row to lpsolve */
                if (!add_constraintex(LPs[poly_index], j, row, colno, LE,
                                      bound)) {
                  ret = 1;
                }
              }
              cpt_ineq1[poly_index]++;
              cpt_ineq[poly_index]++;
            }
          }
        } else {
          bound = GAMMA2 - BETA - 1 - r0.vec[poly_index].coeffs[coeff_index];
          creat_ineq(&c, coeff_index, ineq_line);
          guess = scalar_product(ineq_line, t0_guess[poly_index]);

          if (bound - guess < C * TAU) {
            if (ret == 0) {
              j = 0;
              for (ineq_index = 0; ineq_index < N; ineq_index++) {
                if (ineq_line[ineq_index] != 0) {
                  colno[j] = ineq_index + 1;
                  row[j++] = -ineq_line[ineq_index];
                }
              }

              if (!add_constraintex(LPs[poly_index], j, row, colno, LE,
                                    bound)) {
                ret = 1;
              }
            }
            cpt_ineq0[poly_index]++;
            cpt_ineq[poly_index]++;
          }

          bound = GAMMA2 - BETA - 1 + r0.vec[poly_index].coeffs[coeff_index];
          if (bound + guess < C * TAU) {
            if (ret == 0) {
              j = 0;
              for (ineq_index = 0; ineq_index < N; ineq_index++) {
                if (ineq_line[ineq_index] != 0) {
                  colno[j] = ineq_index + 1;
                  row[j++] = ineq_line[ineq_index];
                }
              }

              if (!add_constraintex(LPs[poly_index], j, row, colno, LE,
                                    bound)) {
                ret = 1;
              }
            }
            cpt_ineq[poly_index]++;
            cpt_ineq0[poly_index]++;
          }
        }
      }
    }
    cpt_signs++;
    // if(cpt_signs%500 == 0){
    //   printf("%d/%d\r", find_min(cpt_ineq), cpt_signs);
    //   fflush(stdout);
    // }
  }

  end = clock();

  cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
  printf("Building %d LPs: %f sec/%d signs/", K, cpu_time_used, cpt_signs);
  for (poly_index = 0; poly_index < K; poly_index++) {
    printf("%d, ", cpt_ineq0[poly_index]);
  }
  printf(" 0 ineq/ ");
  for (poly_index = 0; poly_index < K; poly_index++) {
    printf("%d, ", cpt_ineq1[poly_index]);
  }
  printf(" 1 ineq\n");
  if (ret == 0) {
    for (poly_index = 0; poly_index < K; poly_index++) {
      set_add_rowmode(LPs[poly_index], FALSE);
      sprintf(name_lp_file, "../Lps/%.16s/key%d/poly%d.lp",
              CRYPTO_ALGNAME, sk_index, poly_index);
      write_lp(LPs[poly_index], name_lp_file);
    }

    for (poly_index = 0; poly_index < K; poly_index++) {
      set_verbose(LPs[poly_index], IMPORTANT);
      start = clock();
      /* Now let lpsolve calculate a solution */
      ret = solve(LPs[poly_index]);
      end = clock();
      cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
      printf("Solving LP#%d (code %d): %f sec\n", poly_index, ret,
             cpu_time_used);

      if (ret == OPTIMAL) {
        ret = 0;
      } else {
        ret = 5;
      }
    }
  }

  for (poly_index = 0; poly_index < K; poly_index++) {
    get_variables(LPs[poly_index], t0_guess_updated[poly_index]);
  }

  if (ret == 0) {
    sprintf(t0_guess_filename,
            "../Guess/%.16s/key%d/t0_guess_file.bin",
          CRYPTO_ALGNAME, sk_index);
    t0_guess_file = fopen(t0_guess_filename, "wb");
    if (t0_guess_file == NULL) {
      perror("Error opening file");
      return FILE_OPEN_ERROR;
    }

    fwrite(t0_guess_updated, sizeof(double), K * N, t0_guess_file);
    fclose(t0_guess_file);
  }

  if (row != NULL) { /* free allocated memory */
    free(row);
  }
  if (colno != NULL) {
    free(colno);
  }

  for (poly_index = 0; poly_index < K; poly_index++) {
    if (LPs[poly_index] != NULL) {
      delete_lp(LPs[poly_index]); /* clean up such that all used memory by
                                     lpsolve is freed */
    }
  }

  fclose(sign_compressed_file);
  return 0;
}
