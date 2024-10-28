#ifndef SIGN_RDM_MSG_AND_SAVE_H
#define SIGN_RDM_MSG_AND_SAVE_H

// #include <stddef.h>
#include <stdint.h>
// #include <stdlib.h>
// #include <string.h>
// #include <time.h>

#include "additional_fct.h"
#include "fips202.h"
#include "packing.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "sign.h"
#include "symmetric.h"

#define COMP_CRYPTO_BYTES \
  (CTILDEBYTES + K * POLYZ_PACKEDBYTES + POLYVECH_PACKEDBYTES)

void compute_lb_AZ_minus_ct12d(polyveck *r0_out, poly *c_out, polyveck *h_out,
                               const uint8_t *sig, const uint8_t *m,
                               size_t mlen, const uint8_t *ctx,
                               size_t ctxlen, const uint8_t *pk);

#endif // SIGN_RDM_MSG_AND_SAVE_H