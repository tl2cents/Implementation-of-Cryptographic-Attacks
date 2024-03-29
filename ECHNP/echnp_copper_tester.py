from echnp_oracle import GenECHNP, ECDH, secp256r1_curve, GetECHNPPolys, ECHNP_Polys
from XHS_Lattice import XHS20_Lattice_Polys, XHS22_Lattice_Polys
from coppersmith_reduce import coppersmith_multivarivate_reduced_poly
from sage.all import ZZ, binomial
from rootfind_ZZ import rootfind_ZZ
from findRootsZZ import find_roots_groebner, find_roots_variety, find_roots_gcd, find_roots_resultants


def ECDH_To_ECHNP_Polys(n, k):
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

def test_Coppersmallroot_With_Details_Secp256_XHS20(kbit = 185, n = 5, d = 2):
    E, G = secp256r1_curve()
    p = ZZ(E.base_field().order())
    pbit = p.nbits()
    ebit = pbit - kbit
    msb = True
    print(f"[+] The basic parameters {pbit = } {kbit = } {msb = }")
    print(f"[+] Try to use lattice parameters {n = } {d = } with XHS20 Lattice")
    polys, e, es = ECHNP_Polys(kbit, n)
    roots = [e] + es
    bounds = [2**ebit] + [2**(ebit+1)] * n
    count_poly = lambda d, n : (2*d + 1) * sum(binomial(n, l) for l in range(0,d+1))
    copper_polys = XHS20_Lattice_Polys(n, d, p, polys)
    total_poly_num = count_poly(d, n)
    print(f"[+] Shifted Polynomials test passed {len(copper_polys) = } {total_poly_num == len(copper_polys)}")
    res = []    
    tmp = []
    useful_copper_polys = []
    for poly in copper_polys:
        if poly == 0:
            continue
        useful_copper_polys.append(poly)
        res.append(poly(*roots) % p**d == 0)       
        tmp.append(poly(*roots) % p**(d+1) == 0)
        
    print(f"[+] Check Results (should be all true) {res.count(True) = } { res.count(False) = }")
    print(f"[+] Check tmp (should be all false) {tmp.count(True) = } { tmp.count(False) = }")
    polys = coppersmith_multivarivate_reduced_poly(useful_copper_polys, bounds)
    # check these polys
    good_polys_idx = []
    for i, poly in enumerate(polys):
        poly = poly.change_ring(ZZ)
        assert poly(*roots) % p**d == 0, "not copper poly"
        # check wether the poly has the original roots in ZZ
        if poly(*roots) == 0:
            good_polys_idx.append(i)
    print(f"[+] find {len(good_polys_idx)} polynomials with the desired roots in ZZ")
    print(f"[+] the good polynomials indexes : {good_polys_idx}")
    
    print("[+] try the groebner method to find roots")
    roots = find_roots_groebner(polys[0].parent(), polys)
    for root in roots:        
        print(f"[+] find : {root = }")
        
    print("[+] try the variety method to find roots")
    roots = find_roots_variety(polys[0].parent(), polys)
    for root in roots:
        print(f"[+] find : {root = }")
        
    print("[+] try the gcd method to find roots")
    roots = find_roots_gcd(polys[0].parent(), polys)
    for root in roots:
        print(f"[+] find : {root = }")
    print()
    # print("[+] try the resultant method to find roots")
    # roots = find_roots_resultants(polys[0].parent().gens(), polys)
    # for root in roots:
    #    print(f"[+] find : {root = }")
    
def test_Coppersmallroot_With_Details_Secp256_XHS22(kbit = 172, n = 5, d = 2, t = 1, first_n = 10):
    E, G = secp256r1_curve()
    p = ZZ(E.base_field().order())
    pbit = p.nbits()
    ebit = pbit - kbit

    msb = True
    print(f"[+] The basic parameters {pbit = } {kbit = } {msb = }")
    print(f"[+] Try to use lattice parameters {n = } {d = } {t = } with XHS22 Lattice")
    polys, e, es = ECHNP_Polys(kbit, n)
    roots = [e] + es
    bounds = [2**ebit] + [2**(ebit+1)] * n
    count_poly = lambda n, d, t : (2*d) * sum(binomial(n, l) for l in range(0,d+1)) + (t + 1) * binomial(n, d + 1)
    copper_polys = XHS22_Lattice_Polys(n, d, t, p, polys)
    total_poly_num = count_poly(n, d, t)
    print(f"[+] Shifted Polynomials test passed {len(copper_polys) = } {total_poly_num = }")
    res = []    
    tmp = []
    useful_copper_polys = []
    for poly in copper_polys:
        if poly == 0:
            continue
        useful_copper_polys.append(poly)
        res.append(poly(*roots) % p**d == 0)       
        tmp.append(poly(*roots) % p**(d+1) == 0)
        
    print(f"[+] Check Results (should be all true) {res.count(True) = } { res.count(False) = }")
    print(f"[+] Check tmp (should be all false) {tmp.count(True) = } { tmp.count(False) = }")
    polys = coppersmith_multivarivate_reduced_poly(useful_copper_polys, bounds)
    # check these polys
    good_polys_idx = []
    for i, poly in enumerate(polys):
        poly = poly.change_ring(ZZ)
        assert poly(*roots) % p**d == 0, "not copper poly"
        # check wether the poly has the original roots in ZZ
        if poly(*roots) == 0 and poly % p**d != 0:
            good_polys_idx.append(i)
    print(f"[+] find {len(good_polys_idx)} polynomials with the desired roots in ZZ")
    print(f"[+] the good polynomials indexes : {good_polys_idx}")
    if first_n is None:
        first_n = len(polys)
    final_polys = polys[:first_n]
        
    print(f"[+] select first {first_n} of {len(polys)} polynomials to find the roots in ZZ")
    # print("[+] try the groebner method to find roots")
    # roots = find_roots_groebner(polys[0].parent(), final_polys)
    # for root in roots:        
    #     print(f"[+] find : {root = }")
        
    print("[+] try the multiple method to find roots")
    roots = rootfind_ZZ(final_polys, bounds)
    print(f"[+] find : {roots = }")
    
    print("[+] try the variety method to find roots")
    roots = find_roots_variety(polys[0].parent(), final_polys)
    for root in roots:
        print(f"[+] find : {root = }")
        
    print("[+] try the gcd method to find roots")
    roots = find_roots_gcd(polys[0].parent(), final_polys)
    for root in roots:
        print(f"[+] find : {root = }")
    print()
    
if __name__ == "__main__":
    test_Coppersmallroot_With_Details_Secp256_XHS20()
    test_Coppersmallroot_With_Details_Secp256_XHS22(first_n = 20)
    test_Coppersmallroot_With_Details_Secp256_XHS22(kbit = 165, n = 5, d = 3, t = 2)
    test_Coppersmallroot_With_Details_Secp256_XHS22(kbit = 165, n = 5, d = 3, t = 2)
    