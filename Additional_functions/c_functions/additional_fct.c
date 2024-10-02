#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "additional_fct.h"
#include "fips202.h"
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "rounding.h"
#include "symmetric.h"

#define MAX_MARKER_LEN 50

//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
int FindMarker(FILE *infile, const char *marker) {
  char line[MAX_MARKER_LEN];
  int i, len;
  int curr_line;

  len = (int)strlen(marker);
  if (len > MAX_MARKER_LEN - 1)
    len = MAX_MARKER_LEN - 1;

  for (i = 0; i < len; i++) {
    curr_line = fgetc(infile);
    line[i] = curr_line;
    if (curr_line == EOF)
      return 0;
  }
  line[len] = '\0';

  while (1) {
    if (!strncmp(line, marker, len))
      return 1;

    for (i = 0; i < len - 1; i++)
      line[i] = line[i + 1];
    curr_line = fgetc(infile);
    line[len - 1] = curr_line;
    if (curr_line == EOF)
      return 0;
    line[len] = '\0';
  }

  // shouldn't get here
  return 0;
}

//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
int ReadHex(FILE *infile, unsigned char *a, int Length, char *str) {
  int i, ch, started;
  unsigned char ich;

  if (Length == 0) {
    a[0] = 0x00;
    return 1;
  }
  memset(a, 0x00, Length);
  started = 0;
  if (FindMarker(infile, str))
    while ((ch = fgetc(infile)) != EOF) {
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

      for (i = 0; i < Length - 1; i++)
        a[i] = (a[i] << 4) | (a[i + 1] >> 4);
      a[Length - 1] = (a[Length - 1] << 4) | ich;
    }
  else
    return 0;

  return 1;
}

void fprintBstr(FILE *fp, char *s, unsigned char *a, unsigned long long l) {
  unsigned long long i;

  fprintf(fp, "%s", s);

  for (i = 0; i < l; i++)
    fprintf(fp, "%02X", a[i]);

  if (l == 0)
    fprintf(fp, "00");

  fprintf(fp, "\n");
}

void fprintBstr2(FILE *fp, int b, char *s, unsigned char *a,
                 unsigned long long l) {
  unsigned long long i;

  fprintf(fp, "%s", s);

  for (i = 0; i < l; i++)
    fprintf(fp, "%02X", a[i]);

  if (l == 0)
    fprintf(fp, "00");

  fprintf(fp, ";%d", b);
  fprintf(fp, "\n");
}

void polyw0_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint32_t t[4];
#if GAMMA2 == (Q - 1) / 88
  for (i = 0; i < N / 4; ++i) {
    t[0] = GAMMA2 - a->coeffs[4 * i + 0];
    t[1] = GAMMA2 - a->coeffs[4 * i + 1];
    t[2] = GAMMA2 - a->coeffs[4 * i + 2];
    t[3] = GAMMA2 - a->coeffs[4 * i + 3];

    r[9 * i + 0] = t[0];
    r[9 * i + 1] = t[0] >> 8;
    r[9 * i + 2] = t[0] >> 16;
    r[9 * i + 2] |= t[1] << 2;
    r[9 * i + 3] = t[1] >> 6;
    r[9 * i + 4] = t[1] >> 14;
    r[9 * i + 4] |= t[2] << 4;
    r[9 * i + 5] = t[2] >> 4;
    r[9 * i + 6] = t[2] >> 12;
    r[9 * i + 6] |= t[3] << 6;
    r[9 * i + 7] = t[3] >> 2;
    r[9 * i + 8] = t[3] >> 10;
  }
#elif GAMMA2 == (Q - 1) / 32
  for (i = 0; i < N / 2; ++i) {
    t[0] = GAMMA2 - a->coeffs[2 * i + 0];
    t[1] = GAMMA2 - a->coeffs[2 * i + 1];

    r[5 * i + 0] = t[0];
    r[5 * i + 1] = t[0] >> 8;
    r[5 * i + 2] = t[0] >> 16;
    r[5 * i + 2] |= t[1] << 4;
    r[5 * i + 3] = t[1] >> 4;
    r[5 * i + 4] = t[1] >> 12;
  }
#endif
}

void polyw0_unpack(poly *r, const uint8_t *a) {
  unsigned int i;

#if GAMMA2 == (Q - 1) / 88
  for (i = 0; i < N / 4; ++i) {
    r->coeffs[4 * i + 0] = a[9 * i + 0];
    r->coeffs[4 * i + 0] |= (uint32_t)a[9 * i + 1] << 8;
    r->coeffs[4 * i + 0] |= (uint32_t)a[9 * i + 2] << 16;
    r->coeffs[4 * i + 0] &= 0x3FFFF;

    r->coeffs[4 * i + 1] = a[9 * i + 2] >> 2;
    r->coeffs[4 * i + 1] |= (uint32_t)a[9 * i + 3] << 6;
    r->coeffs[4 * i + 1] |= (uint32_t)a[9 * i + 4] << 14;
    r->coeffs[4 * i + 1] &= 0x3FFFF;

    r->coeffs[4 * i + 2] = a[9 * i + 4] >> 4;
    r->coeffs[4 * i + 2] |= (uint32_t)a[9 * i + 5] << 4;
    r->coeffs[4 * i + 2] |= (uint32_t)a[9 * i + 6] << 12;
    r->coeffs[4 * i + 2] &= 0x3FFFF;

    r->coeffs[4 * i + 3] = a[9 * i + 6] >> 6;
    r->coeffs[4 * i + 3] |= (uint32_t)a[9 * i + 7] << 2;
    r->coeffs[4 * i + 3] |= (uint32_t)a[9 * i + 8] << 10;
    r->coeffs[4 * i + 3] &= 0x3FFFF;

    r->coeffs[4 * i + 0] = GAMMA2 - r->coeffs[4 * i + 0];
    r->coeffs[4 * i + 1] = GAMMA2 - r->coeffs[4 * i + 1];
    r->coeffs[4 * i + 2] = GAMMA2 - r->coeffs[4 * i + 2];
    r->coeffs[4 * i + 3] = GAMMA2 - r->coeffs[4 * i + 3];
  }
#elif GAMMA2 == (Q - 1) / 32
  for (i = 0; i < N / 2; ++i) {
    r->coeffs[2 * i + 0] = a[5 * i + 0];
    r->coeffs[2 * i + 0] |= (uint32_t)a[5 * i + 1] << 8;
    r->coeffs[2 * i + 0] |= (uint32_t)a[5 * i + 2] << 16;
    r->coeffs[2 * i + 0] &= 0xFFFFF;

    r->coeffs[2 * i + 1] = a[5 * i + 2] >> 4;
    r->coeffs[2 * i + 1] |= (uint32_t)a[5 * i + 3] << 4;
    r->coeffs[2 * i + 1] |= (uint32_t)a[5 * i + 4] << 12;
    /* r->coeffs[2*i+1] &= 0xFFFFF; */ /* No effect, since we're anyway at 20
                                          bits */

    r->coeffs[2 * i + 0] = GAMMA2 - r->coeffs[2 * i + 0];
    r->coeffs[2 * i + 1] = GAMMA2 - r->coeffs[2 * i + 1];
  }
#endif
}

/*************************************************
 * Name:        pack_sig_compressed
 *
 * Description: Bit-pack signature sig = (c, z, h).
 *
 * Arguments:   - uint8_t sig[]: output byte array
 *              - const uint8_t *c: pointer to challenge hash length SEEDBYTES
 *              - const polyvecl *z: pointer to vector z
 *              - const polyveck *h: pointer to hint vector h
 **************************************************/
void pack_sig_compressed(uint8_t sig_compressed[MLEN + COMP_CRYPTO_BYTES],
                         const uint8_t sig[CRYPTO_BYTES], const polyveck *r0) {
  unsigned int i;

  for (i = 0; i < CTILDEBYTES; ++i)
    sig_compressed[i] = sig[i];
  sig += CTILDEBYTES;
  sig_compressed += CTILDEBYTES;

  for (i = 0; i < K; ++i)
    polyw0_pack(sig_compressed + i * POLYZ_PACKEDBYTES, &r0->vec[i]);
  sig += L * POLYZ_PACKEDBYTES;
  sig_compressed += K * POLYZ_PACKEDBYTES;

  /* Encode h */
  for (i = 0; i < OMEGA + K; ++i)
    sig_compressed[i] = sig[i];

  sig += OMEGA + K;
  sig_compressed += OMEGA + K;

  for (i = 0; i < MLEN; ++i)
    sig_compressed[i] = sig[i];
}

/*************************************************
 * Name:        unpack_sig_compressed
 *
 * Description: Unpack values used for the attack  sig = (c, r0, h).
 *
 * Arguments:   - uint8_t *c: pointer to output challenge hash
 *              - polyvecl *r0: pointer to output vector r0
 *              - polyveck *h: pointer to output hint vector h
 *              - const uint8_t sig[]: byte array containing
 *                bit-packed signature
 *
 * Returns 1 in case of malformed signature; otherwise 0.
 **************************************************/
int unpack_sig_compressed(
    uint8_t c[CTILDEBYTES], polyveck *r0, polyveck *h,
    const uint8_t sig_compressed[MLEN + COMP_CRYPTO_BYTES]) {
  unsigned int i, j, k;

  for (i = 0; i < CTILDEBYTES; ++i)
    c[i] = sig_compressed[i];
  sig_compressed += CTILDEBYTES;

  for (i = 0; i < K; ++i)
    polyw0_unpack(&r0->vec[i], sig_compressed + i * POLYZ_PACKEDBYTES);
  sig_compressed += K * POLYZ_PACKEDBYTES;

  /* Decode h */
  k = 0;
  for (i = 0; i < K; ++i) {
    for (j = 0; j < N; ++j)
      h->vec[i].coeffs[j] = 0;

    if (sig_compressed[OMEGA + i] < k || sig_compressed[OMEGA + i] > OMEGA)
      return 1;

    for (j = k; j < sig_compressed[OMEGA + i]; ++j) {
      /* Coefficients are ordered for strong unforgeability */
      if (j > k && sig_compressed[j] <= sig_compressed[j - 1])
        return 1;
      h->vec[i].coeffs[sig_compressed[j]] = 1;
    }

    k = sig_compressed[OMEGA + i];
  }

  /* Extra indices are zero for strong unforgeability */
  for (j = k; j < OMEGA; ++j)
    if (sig_compressed[j])
      return 1;

  return 0;
}

void creat_ineq(poly *c, uint8_t coefficient_index, int8_t *ineq_line) {
  uint16_t i;

  for (i = 0; i <= coefficient_index; i++) {
    ineq_line[i] = c->coeffs[coefficient_index - i];
  }

  for (i = coefficient_index + 1; i < N; i++) {
    ineq_line[i] = -c->coeffs[N + coefficient_index - i];
  }
}

uint32_t find_min(uint32_t arr[K]) {
  uint32_t min = arr[0];

  for (int i = 1; i < K; i++) {
    if (arr[i] < min) {
      min = arr[i];
    }
  }
  return min;
}

double scalar_product(int8_t ineq_line[N], double t0_guess_poly[N]) {
  int i;
  double result = 0.0;

  for (i = 0; i < N; i++) {
    result += (t0_guess_poly[i] * ineq_line[i]);
  }
  return -result;
}

double get_min(double val) {
  double up_bound = 1 << (D - 1);
  return (val < up_bound) ? val : up_bound;
}

double get_max(double val) {
  double low_bound = -(1 << (D - 1)) + 1;
  return (val > low_bound) ? val : low_bound;
}
