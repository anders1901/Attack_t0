#ifndef ADDITIONAL_FCT_H
#define ADDITIONAL_FCT_H

#include <ctype.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "poly.h"
#include "polyvec.h"

#define SUCCESS 0
#define FILE_OPEN_ERROR -1
#define DATA_ERROR -2
#define LP_ERROR -3
#define MLEN 32

#define MAX_MARKER_LEN 50

#define COMP_CRYPTO_BYTES                                                      \
  (CTILDEBYTES + K * POLYZ_PACKEDBYTES + POLYVECH_PACKEDBYTES)

int FindMarker(FILE *infile, const char *marker);

int ReadHex(FILE *infile, unsigned char *a, int Length, char *str);

void fprintBstr(FILE *fp, char *s, unsigned char *a, unsigned long long l);

void fprintBstr2(FILE *fp, int b, char *s, unsigned char *a,
                 unsigned long long l);

void polyw0_pack(uint8_t *r, const poly *a);

void polyw0_unpack(poly *r, const uint8_t *a);

void pack_sig_compressed(uint8_t sig_compressed[MLEN + COMP_CRYPTO_BYTES],
                         const uint8_t sig[CRYPTO_BYTES], const polyveck *r0);

int unpack_sig_compressed(
    uint8_t c[CTILDEBYTES], polyveck *r0, polyveck *h,
    const uint8_t sig_compressed[MLEN + COMP_CRYPTO_BYTES]);

void creat_ineq(poly *c, uint8_t coefficient_index, int8_t *ineq_line);

uint32_t find_min(uint32_t arr[K]);

double scalar_product(int8_t ineq_line[N], double t0_guess_poly[N]);

double get_min(double val);

double get_max(double val);

#endif // ADDITIONAL_FCT_H
