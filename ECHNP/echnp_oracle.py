# from ecdsa import NIST256p
# from ecdsa.ecdsa import Public_key
# from ecdsa.ellipticcurve import Point
from sage.all import EllipticCurve, GF, ZZ, Zmod, PolynomialRing, prod, matrix, QQ
import random
from Crypto.Util.number import getPrime

def secp256r1_curve():
    # https://neuromancer.sk/std/nist/P-256
    p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
    K = GF(p)
    a = K(0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc)
    b = K(0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b)
    E = EllipticCurve(K, (a, b))
    G = E(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296, 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)
    E.set_order(0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551 * 0x1)
    return (E, G)
    
class DH():
    """
    A simple Diffie-Hellman and its oracle
    
    """
    def __init__(self, p, g):
        self.p = p
        self.nbit = ZZ(p).nbits()
        self.base_field = GF(p)
        self.g = self.base_field(g)
        self.a = random.randint(1, p-1)
        self.b = random.randint(1, p-1)
        self.A = self.g ** self.a
        self.B = self.g ** self.b
        self.shared_secret = self.g ** (self.a * self.b)

    def get_pubs(self):
        return (ZZ(self.A), ZZ(self.B))

    def get_shared_secret(self):
        return ZZ(self.shared_secret)
    
    def dh(self, C, alice=True):
        """Negotiate a shared secret between two parties

        Args:
            C (Integer): public key of the other party
            alice (bool, optional): If True, this party is Alice. Otherwise, Bob. Defaults to True. 

        Returns:
            Integer: shared secret
        """
        if alice:
            return ZZ(self.base_field(C) ** self.a)
        else:
            return ZZ(self.base_field(C) ** self.b)
    
    def oracle(self, C, kbit, msb=True, alice=True):
        """Oracle for the attacker to leak the dh output

        Args:
            C (Integer): public key of the other party
            kbit (Integer): number of bits to leak
            msb (bool, optional): leak from msb. Defaults to True. If False, leak from lsb.
            alice (bool, optional): If True, this party is Alice. Otherwise, Bob. Defaults to True. 
        """
        shared_secret = self.dh(C, alice)
        if msb:
            return ZZ(shared_secret >> (self.nbit - kbit))
        else:
            return ZZ(shared_secret % (2**kbit))
        
    def oracle_func(self, kbit, msb=True, alice=True):
        # return function for the oracle with fixed k, with input being only `C`
        return lambda C: self.oracle(C, kbit, msb, alice)
    
class ECDH:
    """
    A simple Elliptic Curve Diffie-Hellman and its oracle
    
    """
    
    def __init__(self, curve, G):
        self.curve = curve
        self.G = G
        self.nbit = ZZ(curve.base_field().order()).nbits()
        self.a = random.randint(1, curve.order()-1)
        self.b = random.randint(1, curve.order()-1)
        self.A = self.a * G
        self.B = self.b * G
        self.shared_secret = self.a * self.b * G
        
    def get_pubs(self):
        return (self.A, self.B)
    
    def get_shared_secret(self):
        return ZZ(self.shared_secret[0])
    
    def dh(self, C, alice=True):
        if alice:
            return ZZ((self.a * C)[0])
        else:
            return ZZ((self.b * C)[0])
        
    def oracle(self, C, kbit, msb=True, alice=True):
        shared_secret = self.dh(C, alice)
        if msb:
            return ZZ(shared_secret >> (self.nbit - kbit))
        else:
            return ZZ(shared_secret % (2**kbit))
        
    def oracle_func(self, kbit, msb=True, alice=True):
        # return function for the oracle with fixed k, with input being only `C`
        return lambda C: self.oracle(C, kbit, msb, alice)
    
    
def test_dh():
    print("[+] Testing DH-Oracle Implementation")
    p = getPrime(256)
    while p.bit_length() != 256:
        p = getPrime(256)
    g = 2
    dh = DH(p, g)
    A, B = dh.get_pubs()
    shared_secret = dh.get_shared_secret()
    assert shared_secret == dh.dh(B, True)
    assert shared_secret == dh.dh(A, False)
    print("[+] DH test passed")
    assert dh.oracle_func(128, True, True)(B) == shared_secret >> (dh.nbit - 128)
    assert dh.oracle_func(128, False, True)(B) == shared_secret % (2**128)
    assert dh.oracle_func(128, True, False)(A) == shared_secret >> (dh.nbit - 128)
    assert dh.oracle_func(128, False, False)(A) == shared_secret % (2**128)
    print("[+] DH-Oracle test passed")
    
def test_ecdh():
    print("[+] Testing ECDH-Oracle Implementation")
    E, G = secp256r1_curve()
    ecdh = ECDH(E, G)
    A, B = ecdh.get_pubs()
    shared_secret = ecdh.get_shared_secret()
    assert shared_secret == ecdh.dh(B, True)
    assert shared_secret == ecdh.dh(A, False)
    print("[+] ECDH test passed")
    assert ecdh.oracle_func(128, True, True)(B) == shared_secret >> (ecdh.nbit - 128)
    assert ecdh.oracle_func(128, False, True)(B) == shared_secret % (2**128)
    assert ecdh.oracle_func(128, True, False)(A) == shared_secret >> (ecdh.nbit - 128)
    assert ecdh.oracle_func(128, False, False)(A) == shared_secret % (2**128)
    print("[+] ECDH-Oracle test passed")
    
def EcdhOracle2ECHNP(oracleA, B, Q):
    # oracleA: function to leak the shared secret oracle_a(C) = a*C
    # A, B: public keys of Alice and Bob
    # Q: generator point on the curve
    # returns the ECHNP oralce
    return lambda r : oracleA(B + r * Q)

def GenECHNP(oracleA, CurveParas, n, kbit, msb=True):
    # oracleA: function to leak the shared secret oracle_a(C) = a*C
    # CurveParas:
    # Curve: elliptic curve
    #   p : the modulus of the curve
    #   A, B: public keys of Alice and Bob
    #   Q: generator point on the curve
    # n: number of queries (total 2*n + 1)
    # kbit: number of bits to leak
    # returns the ECHNP samples
    (Curve, p, A, B, Q) = CurveParas
    oracle_ECHNP = EcdhOracle2ECHNP(oracleA, B, Q)
    pbit = ZZ(p).nbits()
    hs = []
    phs = []
    mhs = []
    xqs = []
    h0 = oracle_ECHNP(0) << (int(msb) * (pbit - kbit))
    for i in range(1, n + 1):
        h1 = ZZ(oracle_ECHNP(i) << (int(msb) * (pbit - kbit)))
        h2 = ZZ(oracle_ECHNP(-i) << (int(msb) * (pbit - kbit)))
        xq = ZZ((i * A)[0])
        xqs.append(xq)
        phs.append(h1)
        mhs.append(h2)
        hs.append(h1 + h2)
    return h0, hs, xqs, phs, mhs

def GetECHNPPolys(h0, hs, xqs, CurveParas):
    # hs: leaked bits
    # xqs: x-coordinates of the queries
    # CurveParas:
    # Curve: elliptic curve
    #   p : the modulus of the curve
    #   A, B: public keys of Alice and Bob
    #   Q: generator point on the curve
    # returns the ECHNP polynomials
    assert len(hs) == len(xqs), "Invalid input"
    n = len(hs)
    (Curve, p, A, B, Q) = CurveParas
    p = ZZ(p)
    a, b = ZZ(Curve.a4()), ZZ(Curve.a6())
    Aix = lambda hi, h0, xqi, p, a, b: ZZ((hi*(h0 - xqi)**2  - 2 * h0 ** 2 * xqi - 2*(a + xqi**2)*h0 -\
                                        2*a*xqi - 4*b) % p)
    Bix = lambda hi, h0, xqi, p, a, b: ZZ(2 * (hi *(h0 - xqi) - 2* h0 * xqi - a - xqi**2) % p)
    Cix = lambda hi, h0, xqi, p, a, b: ZZ((hi - 2 * xqi) % p)
    Dix = lambda hi, h0, xqi, p, a, b: ZZ((h0 - xqi)**2 % p)
    Eix = lambda hi, h0, xqi, p, a, b: ZZ(2 * (h0 - xqi) % p)
    PR = PolynomialRing(GF(p), ["x", "y"])
    x, y = PR.gens()
    fxy = lambda hi, h0, xqi: Aix(hi, h0, xqi, p, a, b) * 1 +\
                              Bix(hi, h0, xqi, p, a, b) * x +\
                              Cix(hi, h0, xqi, p, a, b) * x**2 +\
                              Dix(hi, h0, xqi, p, a, b) * y + \
                              Eix(hi, h0, xqi, p, a, b) * x * y +\
                              x**2*y
    polys = []
    for i in range(n):
        hi = ZZ(hs[i])
        xqi = ZZ(xqs[i])
        f = -(2 * (xqi  *(h0 + x) ** 2 + (a + xqi**2)*(h0 + x) + a*xqi + 2*b) - (hi + y) * (h0 + x - xqi) ** 2)
        polys.append(PR(f))
        # polys.append(fxy(hi, h0, xqi))
        assert (fxy(hi, h0, xqi) == f)
    return polys

def test_GetECHNPPolys():
    E, G = secp256r1_curve()
    ecdh = ECDH(E, G)
    A, B = ecdh.get_pubs()
    p = ZZ(E.base_field().order())
    nbit = p.nbits()
    k = 140
    msb = True
    n = 10
    print("[+] Testing ECHNP Polynomials")
    print(f"[+] {msb = }")
    print(f"[+] {n = } {k = }")
    CurveParas = (E, E.base_field().order(), A, B, G)
    oracle = ecdh.oracle_func(k, msb, True)
    h0, hs, xqs, phs, mhs = GenECHNP(oracle, CurveParas, n, k, msb)
    polys = GetECHNPPolys(h0, hs, xqs, CurveParas)
    # full values
    oracle = ecdh.oracle_func(256, msb, True)
    H0, Hs, Xqs, pHs, mHs = GenECHNP(oracle, CurveParas, n, 256, msb)
    assert H0 >> (nbit - k) == h0 >> (nbit - k), "Invalid h0"
    assert all([h >> (nbit - k) == hh  >> (nbit - k) for h, hh in zip(phs, pHs)]), "Invalid phs"
    assert all([h >> (nbit - k) == hh  >> (nbit - k) for h, hh in zip(mhs, mHs)]), "Invalid mhs"
    assert all([x == xx for x, xx in zip(xqs, Xqs)]), "Invalid xqs"

    ebit = nbit - k
    e0 = H0 % 2**(ebit)
    es = [ZZ(h1 % 2**ebit) + ZZ(h2 % 2**ebit) for h1, h2 in zip(pHs, mHs)]
    # check the polynomial small roots
    for i in range(n):
        assert polys[i](e0, es[i]) % p == 0,  "Invalid polynomial"
    print("[+] All MSB ECHNP Polynomials Checked")
    
    msb = False
    print(f"[+] {msb = }")
    print(f"[+] {n = } {k = }")
    oracle = ecdh.oracle_func(k, msb, True)
    h0, hs, xqs, phs, mhs = GenECHNP(oracle, CurveParas, n, k, msb)
    polys = GetECHNPPolys(h0, hs, xqs, CurveParas)
    x, y = polys[0].parent().gens()
    lsb_polys = [poly(x*(2**k), y*(2**k)) for poly in polys]
    # full values
    oracle = ecdh.oracle_func(256, msb, True)
    H0, Hs, Xqs, pHs, mHs = GenECHNP(oracle, CurveParas, n, 256, msb)

    e0 = H0 >> k
    es = [ZZ(h1 >> k) + ZZ(h2 >> k) for h1, h2 in zip(pHs, mHs)]
    # check the polynomial small roots
    for i in range(n):
        assert lsb_polys[i](e0, es[i]) % p == 0,  "Invalid polynomial"
    print("[+] All LSB ECHNP Polynomials Checked")
    print("[+] ECHNP Polynomials test passed")
    
def ECHNP_Polys(k = 150, n = 10):
    E, G = secp256r1_curve()
    ecdh = ECDH(E, G)
    A, B = ecdh.get_pubs()
    p = ZZ(E.base_field().order())
    nbit = p.nbits()
    ebit = nbit - k
    msb = True
    CurveParas = (E, E.base_field().order(), A, B, G)
    oracle = ecdh.oracle_func(k, msb, True)
    h0, hs, xqs, phs, mhs = GenECHNP(oracle, CurveParas, n, k, msb)
    polys = GetECHNPPolys(h0, hs, xqs, CurveParas)
    # full values
    oracle = ecdh.oracle_func(256, msb, True)
    H0, Hs, Xqs, pHs, mHs = GenECHNP(oracle, CurveParas, n, 256, msb)
    assert H0 >> (nbit - k) == h0 >> (nbit - k), "Invalid h0"
    assert all([h >> (nbit - k) == hh  >> (nbit - k) for h, hh in zip(phs, pHs)]), "Invalid phs"
    assert all([h >> (nbit - k) == hh  >> (nbit - k) for h, hh in zip(mhs, mHs)]), "Invalid mhs"
    assert all([x == xx for x, xx in zip(xqs, Xqs)]), "Invalid xqs"

    e0 = H0 % 2**(ebit)
    es = [ZZ(h1 % 2**ebit) + ZZ(h2 % 2**ebit) for h1, h2 in zip(pHs, mHs)]
    # check the polynomial small roots
    for i in range(n):
        assert polys[i](e0, es[i]) % p == 0,  "Invalid polynomial"
    return polys, e0, es

    
if __name__ == "__main__":
    test_dh()
    test_ecdh()
    test_GetECHNPPolys()
    