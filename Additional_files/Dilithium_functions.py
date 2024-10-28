from Dilithium_parameters import *

try:
    from Crypto.Hash import SHAKE128
except ModuleNotFoundError:
    print("Module Crypto not found, trying with Cryptodome")
    try:
        from Cryptodome.Hash import SHAKE128            
    except ModuleNotFoundError:
        print("Module Crypto AND Cryptodom not found, please verify the installation")

from binascii import unhexlify


def montgomery_reduce(a):
    """
    *************************************************
    * Description: For finite field element a with 0 <= a <= Q*2^32,
    *              compute r \equiv a*2^{-32} (mod Q) such that 0 <= r < 2*Q.
    *
    * Arguments:   - int a: finite field element
    *
    * Returns:     - int r: such as in the Description.
    **************************************************
    """
    r = (a * MONT_INV) % Q
    return r - Q if r > Q_HALF else r


def reduce32(a):
    """
    *************************************************
    * Description: For finite field element a, compute r \equiv a (mod Q)
    *              such that 0 <= r < 2*Q.
    *
    * Arguments:   - int a: finite field element
    *
    * Returns:     - int r: such as in the Description.
    **************************************************
    """
    MAX_VAL = 2**31 - 2**22 - 1
    MIN_VAL = -2**31 + 2**22
    t = (a + (1 << 22)) >> 23  # Calculate t as floor((a + 2^22) / 2^23)
    t = a - t * Q  # Calculate t as a - t * Q
    t = ((t + Q) & MAX_VAL) - Q if t > MAX_VAL else t  # Check if t is greater than the maximum value, if so, subtract Q, else if t is less than the minimum value, add Q   
    return t


def poly_reduce(a):
    """
    *************************************************
    * Description: Inplace reduction of all coefficients of
    *              input polynomial to representative in [0,2*Q[.
    *
    * Arguments:   - array[N](int) a: input/output polynomial
    **************************************************
    """
    for i in range(N):
        a[i] = reduce32(a[i])


def polyt1_unpack(a):
    """
    *************************************************
    * Description: Unpack polynomial t1 with 9-bit coefficients.
    *              Output coefficients are standard representatives.
    *
    * Arguments:   - bytes[POLT1_SIZE_PACKED] a: byte array with bit-packed polynomial
    *
    * Returns      - array[N](int) r: output polynomial
    **************************************************
    """
    r = [0]*N
    for i in range(N//4):
        r[4*i+0] = ((a[5*i+0] >> 0) | (a[5*i+1] << 8)) & 0x3FF
        r[4*i+1] = ((a[5*i+1] >> 2) | (a[5*i+2] << 6)) & 0x3FF
        r[4*i+2] = ((a[5*i+2] >> 4) | (a[5*i+3] << 4)) & 0x3FF
        r[4*i+3] = ((a[5*i+3] >> 6) | (a[5*i+4] << 2)) & 0x3FF
    return r


def polyt0_unpack(a):
    """
    *************************************************
    * Description: Unpack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
    *              Output coefficients lie in ]Q-2^{D-1},Q+2^{D-1}].
    *
    * Arguments:   - bytes[POLT0_SIZE_PACKED] a: byte array with bit-packed polynomial
    *
    * Returns:     - array[N](int) r: output polynomial
    **************************************************
    """
    r = [0]*N
    for i in range(N//8):
        r[8*i+0]  = a[13*i+0]
        r[8*i+0] |= a[13*i+1] << 8
        r[8*i+0] &= 0x1FFF

        r[8*i+1]  = a[13*i+1] >> 5
        r[8*i+1] |= a[13*i+2] << 3
        r[8*i+1] |= a[13*i+3] << 11
        r[8*i+1] &= 0x1FFF

        r[8*i+2]  = a[13*i+3] >> 2
        r[8*i+2] |= a[13*i+4] << 6
        r[8*i+2] &= 0x1FFF

        r[8*i+3]  = a[13*i+4] >> 7
        r[8*i+3] |= a[13*i+5] << 1
        r[8*i+3] |= a[13*i+6] << 9
        r[8*i+3] &= 0x1FFF

        r[8*i+4]  = a[13*i+6] >> 4
        r[8*i+4] |= a[13*i+7] << 4
        r[8*i+4] |= a[13*i+8] << 12
        r[8*i+4] &= 0x1FFF

        r[8*i+5]  = a[13*i+8] >> 1
        r[8*i+5] |= a[13*i+9] << 7
        r[8*i+5] &= 0x1FFF

        r[8*i+6]  = a[13*i+9] >> 6
        r[8*i+6] |= a[13*i+10] << 2
        r[8*i+6] |= a[13*i+11] << 10
        r[8*i+6] &= 0x1FFF

        r[8*i+7]  = a[13*i+11] >> 3
        r[8*i+7] |= a[13*i+12] << 5
        r[8*i+7] &= 0x1FFF


        r[8*i+0] =  (1 << (D-1)) - r[8*i+0]
        r[8*i+1] =  (1 << (D-1)) - r[8*i+1]
        r[8*i+2] =  (1 << (D-1)) - r[8*i+2]
        r[8*i+3] =  (1 << (D-1)) - r[8*i+3]
        r[8*i+4] =  (1 << (D-1)) - r[8*i+4]
        r[8*i+5] =  (1 << (D-1)) - r[8*i+5]
        r[8*i+6] =  (1 << (D-1)) - r[8*i+6]
        r[8*i+7] =  (1 << (D-1)) - r[8*i+7]

    return r


def polyeta_unpack(a):
    """
    *************************************************
    * Description: Unpack polynomial with coefficients in [-ETA,ETA].
    *              Output coefficients lie in [Q-ETA,Q+ETA].
    *
    * Arguments:   - bytes[POLETA_SIZE_PACKED] a: byte array with bit-packed polynomial
    *
    * Returns:     - array[N](int) r: output polynomial
    **************************************************
    """
    r = [0]*N
    if ETA == 2:
        for i in range(N//8):
            r[8*i+0] = a[3*i+0] & 0x07
            r[8*i+1] = (a[3*i+0] >> 3) & 0x07
            r[8*i+2] = ((a[3*i+0] >> 6) | (a[3*i+1] << 2)) & 0x07
            r[8*i+3] = (a[3*i+1] >> 1) & 0x07
            r[8*i+4] = (a[3*i+1] >> 4) & 0x07
            r[8*i+5] = ((a[3*i+1] >> 7) | (a[3*i+2] << 1)) & 0x07
            r[8*i+6] = (a[3*i+2] >> 2) & 0x07
            r[8*i+7] = (a[3*i+2] >> 5) & 0x07

            r[8*i+0] = ETA - r[8*i+0]
            r[8*i+1] = ETA - r[8*i+1]
            r[8*i+2] = ETA - r[8*i+2]
            r[8*i+3] = ETA - r[8*i+3]
            r[8*i+4] = ETA - r[8*i+4]
            r[8*i+5] = ETA - r[8*i+5]
            r[8*i+6] = ETA - r[8*i+6]
            r[8*i+7] = ETA - r[8*i+7]

    else:
        for i  in range(N//2):
            r[2*i+0] = a[i] & 0x0F
            r[2*i+1] = a[i] >> 4
            r[2*i+0] = ETA - r[2*i+0]
            r[2*i+1] = ETA - r[2*i+1]

    return r


def unpack_pk(pk):
    """
    *************************************************
    * Description: Unpack public key pk = (rho, t1).
    *
    * Arguments:   - str(hex) pk: string of hex values containing bit-packed pk
    *
    * Returns:     - bytes[SEEDBYTES] rho: output byte array for rho
    *              - array[K][N](int) t1: output vector t1
    **************************************************
    """
    offset = 0
    # rho
    rho = unhexlify(pk[:SEEDBYTES])
    offset = SEEDBYTES

    # t1
    t1 = [ polyt1_unpack(unhexlify(pk[offset + index : offset + index + POLT1_SIZE_PACKED*2])) for index in range(0, (POLT1_SIZE_PACKED*2)*K, (POLT1_SIZE_PACKED*2))]

    return rho, t1


def unpack_sk(sk):
    """
    *************************************************
    * Description: Unpack secret key sk = (rho, key, tr, s1, s2, t0).
    *
    * Arguments:   - str(hex) sk: input byte array
    *
    * Returns:     - bytes[SEEDBYTES] rho: output byte array for rho
    *              - str(hex) key: string of hex values containing key
    *              - str(hex) tr: string of hex values containing tr
    *              - array[L][N](int) s1: vector s1
    *              - array[K][N](int) s2: vector s2
    *              - array[K][N](int) t0: vector t1
    **************************************************
    """
    offset = 0
    # rho
    rho = unhexlify(sk[:SEEDBYTES])
    offset = SEEDBYTES

    # key
    key = sk[offset : offset + SEEDBYTES]
    offset += SEEDBYTES

    # tr
    tr = sk[offset : offset + 2*TRBYTES]
    offset += 2*TRBYTES


    # s1
    s1 = [ polyeta_unpack(unhexlify(sk[offset + index : offset + index + POLETA_SIZE_PACKED*2])) for index in range(0, (POLETA_SIZE_PACKED*2)*L, (POLETA_SIZE_PACKED*2))]
    offset += (POLETA_SIZE_PACKED*2)*L

    # s2
    s2 = [ polyeta_unpack(unhexlify(sk[offset + index : offset + index + POLETA_SIZE_PACKED*2])) for index in range(0, (POLETA_SIZE_PACKED*2)*K, (POLETA_SIZE_PACKED*2))]
    offset += (POLETA_SIZE_PACKED*2)*K

    # t0
    t0 = [ polyt0_unpack(unhexlify(sk[offset + index : offset + index + POLT0_SIZE_PACKED*2])) for index in range(0, (POLT0_SIZE_PACKED*2)*K, (POLT0_SIZE_PACKED*2))]

    return rho, key, tr, s1, s2, t0


def rej_uniform(buf, buflen, len_ = N):
    """
    *************************************************
    * Description: Sample uniformly random coefficients in [0, Q-1] by
    *              performing rejection sampling using array of random bytes.
    *
    * Arguments:   - array[len_](int) A: output array (declared outside)
    *              - unsigned int len_: number of coefficients to be sampled (default: N)
    *              - str[buflen] buf: array of random bytes
    *              - int buflen: length of array of random bytes
    *
    * Returns:     - int ctr: number of sampled coefficients. Can be smaller than len_ if not enough
    *                random bytes were given.
    **************************************************
    """
    A_ = []
    ctr, pos = 0, 0
    while(ctr < len_ and pos + 3 <= buflen):
        t  = buf[pos]
        pos+= 1
        t |= (buf[pos] << 8)
        pos+=1
        t |= (buf[pos] << 16)
        pos+= 1
        t &= 0x7FFFFF

        if(t < Q):
            A_.append(t)
            ctr+=1
    return ctr, A_


def poly_uniform(seed, nonce):
    """
    *************************************************
    * Description: Sample polynomial with uniformly random coefficients
    *              in [-ETA,ETA] by performing rejection sampling using the
    *              output stream from SHAKE128(seed|nonce).
    *
    * Arguments:   - bytes[SEEDBYTES] seed: byte array with seed
    *              - int nonce: 2-byte nonce
    *
    * Returns:     - array[N](int) S: output polynomial
    **************************************************
    """
    ctr = 0
    nblocks = ((768 + STREAM128_BLOCKBYTES - 1)//STREAM128_BLOCKBYTES)

    buflen = nblocks*STREAM128_BLOCKBYTES
    
    m = bytearray(34)
    for i in range(32):
        m[i] = seed[i]
        
    m[33] = nonce >> 8
    m[32] = nonce ^ (m[33]<<8 )


    # initialise shake
    shake = SHAKE128.new(m)

    S = []
    out = b""
    
    # First Squeeze Block
    while (nblocks > 0):
        out+= shake.read(SHAKE128_RATE)
        nblocks -= 1

    ctr, S = rej_uniform(out, buflen)

    # if not enough random bytes to fill the polynomial coefficients
    out = bytearray(out)
    
    while(ctr < N):
        off = buflen % 3
        for i in range(off):
            out[i] = out[buflen - off + i]

        buflen = STREAM128_BLOCKBYTES + off
        nblocks = 1
        while (nblocks > 0):
            out[off:] = shake.read(SHAKE128_RATE)
            nblocks -= 1

        ctr_, a = rej_uniform(out, buflen, len_ = N - ctr)
        ctr += ctr_
        S += a
    return S


def polyvec_matrix_expand(rho):
    """
    **************************************************
    * Description: Implementation of ExpandA. Generates matrix A with uniformly
    *              random coefficients a_{i,j} by performing rejection
    *              sampling on the output stream of SHAKE128(rho|i|j).
    *              warning: C version uses poly_uniform function not implemented here (maybe code it ?)
    *
    * Arguments:   - str rho: byte array containing seed rho
    *
    * Returns:     - array[K][L][N](int) A: output matrix
    **************************************************
    """
    A = []
    for i in range(K):
        A_ = []
        for j in range(L):
            A_.append(poly_uniform(rho,  (i << 8) + j))
        A.append(A_)
    return A


def invntt_frominvmont(p):
    """
    *************************************************
    * Description: backward NTT, in-place. No modular reduction is performed after
    *              additions or subtractions. Hence output coefficients can be up
    *              to 16*Q larger than the coefficients of the input polynomial.
    *              Output vector is in bitreversed order.
    *              Elements of p must be at least 32 bits long or else everflow occurs.
    *
    * Arguments:   - array[N](int) p: input/output coefficient array
    **************************************************
    """
    k = 256
    l = 1
    f = 41978
    
    while l < 256:
        start = 0
        while start < 256:
            k -= 1
            zeta = -zetas[k]
            for j in range(start, start + l):
                t        = p[j]
                p[j]     = t + p[j + l]
                p[j + l] = t - p[j + l]
                p[j + l] = montgomery_reduce(zeta * p[j + l])
            start = j + l + 1
        l = l << 1
        
    for j in range(256):
        p[j] = montgomery_reduce(f * p[j])