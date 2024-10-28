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

#include "sign_rdm_msg_and_save.h"

void compute_lb_AZ_minus_ct12d(polyveck *r0_out, poly *c_out, polyveck *h_out,
                               const uint8_t *sig, const uint8_t *m,
                               size_t mlen, const uint8_t *ctx,
                               size_t ctxlen, const uint8_t *pk)
{
  unsigned int i, j;
  uint8_t rho[SEEDBYTES];
  uint8_t mu[CRHBYTES];
  uint8_t c[CTILDEBYTES];
  poly cp;
  polyvecl mat[K], z;
  polyveck Az, t1, r, h, r0;
  keccak_state state;

  unpack_pk(rho, &t1, pk);
  if (unpack_sig(c, &z, &h, sig))
  {
    printf("problem\n");
  }

  for (i = 0; i < K; i++)
  {
    for (j = 0; j < N; j++)
    {
      h_out->vec[i].coeffs[j] = h.vec[i].coeffs[j];
    }
  }

  /* Compute CRH(H(rho, t1), msg) */
  shake256(mu, TRBYTES, pk, CRYPTO_PUBLICKEYBYTES);
  shake256_init(&state);
  shake256_absorb(&state, mu, TRBYTES);
  mu[0] = 0;
  mu[1] = ctxlen;
  shake256_absorb(&state, mu, 2);
  shake256_absorb(&state, ctx, ctxlen);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  /* Matrix-vector multiplication; compute Az - c2^dt1 */
  poly_challenge(&cp, c);
  for (i = 0; i < N; i++)
  {
    c_out->coeffs[i] = cp.coeffs[i];
  }

  polyvec_matrix_expand(mat, rho);

  polyvecl_ntt(&z);
  polyvec_matrix_pointwise_montgomery(&Az, mat, &z);

  poly_ntt(&cp);
  polyveck_shiftl(&t1);
  polyveck_ntt(&t1);
  polyveck_pointwise_poly_montgomery(&t1, &cp, &t1);

  polyveck_sub(&r, &Az, &t1);
  polyveck_reduce(&r);
  polyveck_invntt_tomont(&r);

  /* Reconstruct w1 */
  polyveck_caddq(&r);
  polyveck_decompose(&r, &r0, &r);
  for (i = 0; i < K; i++)
  {
    for (j = 0; j < N; j++)
    {
      r0_out->vec[i].coeffs[j] = r0.vec[i].coeffs[j];
    }
  }
}

int main(int argc, char const *argv[])
{
  int NB_SIGNS, sk_index;

  if (argc == 3)
  {
    // Case 1: Two arguments   (int NB_inequalities, float C)
    sk_index = atoi(argv[1]);
    NB_SIGNS = atoi(argv[2]);
  }
  else
  {
    // Invalid number of arguments
    printf("Usage:\n");
    printf("  %s <int> <int>\n", argv[0]);
    return DATA_ERROR;
  }

  char directory[32];
  char sign_compressed_file_name[128];
  FILE *sign_compressed_file;
  uint8_t sm_compressed[MLEN + COMP_CRYPTO_BYTES];

  char fn_rsp[90];
  FILE *fp_rsp;
  int i;

  size_t smlen;
  uint8_t m[MLEN + CRYPTO_BYTES];
  uint8_t sm[MLEN + CRYPTO_BYTES];
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];

  polyveck r0, h;
  poly c;

  int count = 0;

  struct stat st = {0};

  sprintf(fn_rsp, "../../dilithium/ref/PQCsignKAT_%.16s.rsp",
          CRYPTO_ALGNAME);
  if ((fp_rsp = fopen(fn_rsp, "r")) == NULL)
  {
    printf("Couldn't open <%s> for read\n", fn_rsp);
    return DATA_ERROR;
  }

  sprintf(directory, "../../signs/");
  if (stat(directory, &st) == -1)
  {
    mkdir(directory, 0700);
  }

  sprintf(directory, "../../signs/%.16s/", CRYPTO_ALGNAME);
  if (stat(directory, &st) == -1)
  {
    mkdir(directory, 0700);
  }

  sprintf(directory, "../../signs/%.16s/key%d/", CRYPTO_ALGNAME, sk_index);
  if (stat(directory, &st) == -1)
  {
    mkdir(directory, 0700);
  }

  sprintf(sign_compressed_file_name,
          "../../signs/%.16s/key%d/PQCsignKAT_%.16s_compressed.rsp",
          CRYPTO_ALGNAME, sk_index, CRYPTO_ALGNAME);
  if ((sign_compressed_file = fopen(sign_compressed_file_name, "w")) == NULL)
  {
    printf("Couldn't open <%s> for write\n", sign_compressed_file_name);
    return DATA_ERROR;
  }

  do
  {
    if (FindMarker(fp_rsp, "count = "))
    {
      if (fscanf(fp_rsp, "%d", &count) != 1)
      {
        printf("Parse error\n");
        exit(-1);
      }
    }
    else
    {
      break;
    }

    if (!ReadHex(fp_rsp, pk, (int)CRYPTO_PUBLICKEYBYTES, "pk = "))
    {
      printf("ERROR: unable to read 'pk' from <%s>\n", fn_rsp);
      return DATA_ERROR;
    }
    if (!ReadHex(fp_rsp, sk, (int)CRYPTO_SECRETKEYBYTES, "sk = "))
    {
      printf("ERROR: unable to read 'sk' from <%s>\n", fn_rsp);
      return DATA_ERROR;
    }

  } while (count < sk_index);

  for (i = 0; i < NB_SIGNS; i++)
  {
    randombytes(m, MLEN);                          /* Random MLEN = 32 bytes msg */
    crypto_sign(sm, &smlen, m, MLEN, NULL, 0, sk); /* Get standard signature */
    compute_lb_AZ_minus_ct12d(&r0, &c, &h, sm, m, MLEN,
                              NULL, 0, pk); /* Retreive LowBits(A.z - c.t_{1}.2^{d}) */
    pack_sig_compressed(sm_compressed, sm, &r0);
    fprintBstr(sign_compressed_file, "", sm_compressed,
               MLEN + COMP_CRYPTO_BYTES);
    if (i % 1000 == 0)
    {
      printf("%d/%d\r", i, NB_SIGNS);
      fflush(stdout);
    }
  }
  fclose(sign_compressed_file);
  return SUCCESS;
}
