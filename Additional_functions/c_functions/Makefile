DILITHIUM_PATH = ../../dilithium-master/ref
LPLIB_PATH ?= ../..

CC ?= /usr/bin/cc
CFLAGS += -Wall -Wextra -Wmissing-prototypes -Wno-unused-variable\
  -Wshadow -Wvla -Wpointer-arith -O3 -fomit-frame-pointer \
  -I$(DILITHIUM_PATH)
SOURCES = sign.c packing.c polyvec.c poly.c ntt.c reduce.c rounding.c randombytes.c
HEADERS = config.h params.h api.h sign.h packing.h polyvec.h poly.h ntt.h \
  reduce.h rounding.h symmetric.h randombytes.h randombytes.h
KECCAK_SOURCES = $(SOURCES) fips202.c symmetric-shake.c
KECCAK_HEADERS = $(HEADERS) fips202.h

ALL_SOURCES = $(patsubst %,$(DILITHIUM_PATH)/%,$(KECCAK_SOURCES)) additional_fct.c
ALL_HEADERS = $(patsubst %,$(DILITHIUM_PATH)/%,$(KECCAK_HEADERS)) additional_fct.h

sign_rdm_msg_and_save2: sign_rdm_msg_and_save.c $(ALL_SOURCES) $(ALL_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=2 -o $@ $< $(ALL_SOURCES)

sign_rdm_msg_and_save3: sign_rdm_msg_and_save.c $(ALL_SOURCES) $(ALL_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=3 -o $@ $< $(ALL_SOURCES)

sign_rdm_msg_and_save5: sign_rdm_msg_and_save.c $(ALL_SOURCES) $(ALL_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=5 -o $@ $< $(ALL_SOURCES)

build_solve_t0_lp2: build_solve_t0_lp.c $(ALL_SOURCES) $(ALL_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=2 -o $@ $< $(ALL_SOURCES) -I$(LPLIB_PATH)/lp_solve_5.5/ -L$(LPLIB_PATH)/lp_solve_5.5/lpsolve55/bin/ux64 -llpsolve55

build_solve_t0_lp3: build_solve_t0_lp.c $(ALL_SOURCES) $(ALL_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=3 -o $@ $< $(ALL_SOURCES) -I/home/amina/Documents/kali/lp_solve_5.5/ -L/home/amina/Documents/kali/lp_solve_5.5/lpsolve55/bin/ux64 -llpsolve55

build_solve_t0_lp5: build_solve_t0_lp.c $(ALL_SOURCES) $(ALL_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=5 -o $@ $< $(ALL_SOURCES) -I/home/amina/Documents/kali/lp_solve_5.5/ -L/home/amina/Documents/kali/lp_solve_5.5/lpsolve55/bin/ux64 -llpsolve55

clean:
	rm -f sign_rdm_msg_and_save2 sign_rdm_msg_and_save3 sign_rdm_msg_and_save5 build_solve_t0_lp2 build_solve_t0_lp3 build_solve_t0_lp5

format:
	clang-format-15 -i sign_rdm_msg_and_save.[ch] additional_fct.[ch] build_solve_t0_lp.[ch]

.PHONY: clean format