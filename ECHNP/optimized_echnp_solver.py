from echnp_oracle import GenECHNP, ECDH, secp256r1_curve, GetECHNPPolys, ECHNP_Polys
from XHS_Lattice import XHS20_Lattice_Polys, XHS22_Lattice_Polys
from coppersmith_reduce import coppersmith_multivarivate_reduced_poly
from sage.all import ZZ, binomial, inverse_mod
from rootfind_ZZ import rootfind_ZZ
from findRootsZZ import find_roots_groebner, find_roots_variety, find_roots_gcd, find_roots_resultants

def echnp_coppersmith_solver_groebner_op(Curve, G, pubA, pubB, kbit, H0, positiveH, negativeH, xQ, d, t=1, 
                             msb=True, lattice="XHS22", first_n=None):
    """solve the ECHNP problem using coppersmith method

    Args:
        Curve (EllipticCurve): the elliptic curve used in the ECDH oracle
        G (point in EllipticCurve): the generator point
        pubA (point in EllipticCurve): the public key of Alice
        pubB (point in EllipticCurve): the public key of Bob
        kbit (int or integer): the number of bits to leak
        H0 (int or integer): the leaked bits of the shared secret of the oracle at 0 
                            i.e. msb leak of ([ab]G + [0]pubA) or ([ab]G + [0]pubB)
        positiveH (list of ints or integers): the leak bits of oracle(m) with positive m = 1,2..., n 
                            i.e. msb leak of (abG + [m]pubA) or ([ab]G + [m]pubB)
        negativeH (list of ints or integers): the leak bits of oracle(-m) with negative m = 1,2..., n 
                            i.e. msb leak of (abG + [-m]pubA) or ([ab]G + [-m]pubB)
        xQ (list of ints or integers) : the x-coordinates of Q = [m]pubA for m = 1,2,...,n when oracle is based on privA i.e. (abG + [m]pubA) 
                                        the x-coordinates of Q = [m]pubB for m = 1,2,...,n when oracle is based on privB i.e. (abG + [m]pubB)
        d (int or integer): the `d` parameter for constructing the lattice
        t (int or integer): the `t` parameter for constructing the lattice (for XHS22). Default is 1.
        msb (bool, optional): the most significant bit leak or least significant bit leak. Defaults to True i.e MSB leak model.
            **Remarks** if msb is True, the value of H0, positiveH, negativeH should be the msb leak with bit in its original position : (x>> (pbit - kbit)) << (pbit - kbit)
        lattice (string) : the lattice model used : "XHS20" or "XHS22". Default is "XHS22"
        first_n (int or integer) : use the first n polynomials to generate groebner basis. Defaul is all.
    Returns:
        integer or None: the x-coordinate of the shared secret [ab]G if found, None otherwise
    """
    assert len(positiveH) == len(negativeH), "Invalid input"
    assert lattice in ["XHS20", "XHS22"], "unknown lattice"
    R = Curve.base_ring()
    p = R.order()
    pbit = p.nbits()
    ebit = pbit - kbit
    n = len(positiveH)
    CurveParas = (Curve, p, pubA, pubB, G)
    if msb:
        # check form (x>> (pbit - kbit)) << (pbit - kbit)
        if any(h % 2**ebit != 0 for h in [H0] + positiveH + negativeH):
            # transform the msb leak to the form (x>> (pbit - kbit)) << (pbit - kbit)
            H0 = H0 << ebit
            positiveH = [h << ebit for h in positiveH]
            negativeH = [h << ebit for h in negativeH]
    Hs = [ZZ(ph) + ZZ(nh)  for ph, nh in zip(positiveH, negativeH)]
    polys = GetECHNPPolys(H0, Hs, xQ, CurveParas)
    if not msb:
        # lsb case
        x, y = polys[0].parent().gens()
        shift = 2 ** kbit
        # make the leading coeffient of x^2y being zero
        inv_shift = ZZ(inverse_mod(shift ** 3, p))
        polys = [poly(shift *x, shift * y) * inv_shift for poly in polys]
        
    # recentering optimization  
    ebit = ebit - 1
    bounds = [2**ebit] + [2**(ebit+1)] * n
    new_poly = []
    for poly in polys:
        pr = poly.parent()
        xs = list(pr.gens())
        shift_xs = {x: x + b for x,b in zip(xs, bounds)}
        new_poly.append(pr(poly.subs(shift_xs)))
    polys = new_poly[:]
    
    # generate the lattice polynomials
    if lattice == "XHS22":
        copper_polys = XHS22_Lattice_Polys(n, d, t, p, polys)
    else:
        copper_polys = XHS20_Lattice_Polys(n, d, p, polys)
    
    print(f"[+] Generating {len(copper_polys)} polynomials for coppersmith ({lattice})")
    reduced_polys = coppersmith_multivarivate_reduced_poly(copper_polys, bounds)
    if first_n == None:
        first_n = len(reduced_polys)
    print(f"[+] try the groebner method to find roots using {first_n} of {len(reduced_polys)} polynomials")
    roots = find_roots_groebner(reduced_polys[0].parent(), reduced_polys[:first_n])
    for root in roots:
        print(f"[+] find : {root = }")
        # returns only one root
        if type(root) == list:
            e0 = root[0] + 2**ebit
        elif type(root) == dict:
            xs = reduced_polys[0].parent().gens()
            e0 = ZZ(root[xs[0]]) + 2**ebit
        else:
            assert False, "unknown root type"
        if msb:
            return ZZ(H0 + e0)
        else:
            return ZZ(H0 + (e0 << kbit))
    return None

def test_echnp_coppersmith_solver_groebner_secp256_op(kbit, n, d, t=1, lattice="XHS22", msb=True, first_n=None):
    E, G = secp256r1_curve()
    ecdh = ECDH(E, G)
    A, B = ecdh.get_pubs()
    shared_secret = ecdh.get_shared_secret()
    print(f"[+] The shared secret is {shared_secret}")
    p = ZZ(E.base_field().order())
    pbit = p.nbits()
    ebit = pbit - kbit
    print(f"[+] The basic parameters {pbit = } {kbit = } {msb = }")
    print(f"[+] Try to use lattice parameters {n = } {d = } {t = } with {lattice} Lattice")
    CurveParas = (E, p, A, B, G)
    oracle = ecdh.oracle_func(kbit, msb, True)
    h0, hs, xqs, phs, nhs = GenECHNP(oracle, CurveParas, n, kbit, msb)
    result = echnp_coppersmith_solver_groebner_op(E, G, A, B, kbit, h0, phs, nhs, xqs, 
                                               d, t, msb, lattice, first_n)
    if result is not None:
        print(f"[+] The shared secret : {result = }\n[+] check : {result == shared_secret}")
    else:
        print("[+] The shared secret is not found")
        
def echnp_coppersmith_solver_optimized(Curve, G, pubA, pubB, kbit, H0, positiveH, negativeH, xQ, d, t=1, 
                             msb=True, lattice="XHS22", first_n=None):
    """solve the ECHNP problem using coppersmith method with recentering optimization

    Args:
        Curve (EllipticCurve): the elliptic curve used in the ECDH oracle
        G (point in EllipticCurve): the generator point
        pubA (point in EllipticCurve): the public key of Alice
        pubB (point in EllipticCurve): the public key of Bob
        kbit (int or integer): the number of bits to leak
        H0 (int or integer): the leaked bits of the shared secret of the oracle at 0 i.e. msb leak of ([ab]G + [0]pubA) or ([ab]G + [0]pubB)
        positiveH (list of ints or integers): the leak bits of oracle(m) with positive m = 1,2..., n i.e. msb leak of (abG + [m]pubA) or ([ab]G + [m]pubB)
        negativeH (list of ints or integers): the leak bits of oracle(-m) with negative m = 1,2..., n i.e. msb leak of (abG + [-m]pubA) or ([ab]G + [-m]pubB)
        xQ (list of ints or integers) : the x-coordinates of Q = [m]pubA for m = 1,2,...,n when oracle is based on privA i.e. (abG + [m]pubA) 
                                        the x-coordinates of Q = [m]pubB for m = 1,2,...,n when oracle is based on privB i.e. (abG + [m]pubB)
        d (int or integer): the `d` parameter for constructing the lattice
        t (int or integer): the `t` parameter for constructing the lattice (for XHS22). Default is 1.
        msb (bool, optional): the most significant bit leak or least significant bit leak. Defaults to True i.e MSB leak model.
            **Remarks** if msb is True, the value of H0, positiveH, negativeH should be the msb leak with bit in its original position : (x>> (pbit - kbit)) << (pbit - kbit)
        lattice (string) : the lattice model used : "XHS20" or "XHS22". Default is "XHS22"
        first_n (int or integer) : use the first n polynomials to generate groebner basis. Defaul is all.
    Returns:
        integer or None: the x-coordinate of the shared secret [ab]G if found, None otherwise
    """
    assert len(positiveH) == len(negativeH), "Invalid input"
    assert lattice in ["XHS20", "XHS22"], "unknown lattice"
    R = Curve.base_ring()
    p = R.order()
    pbit = p.nbits()
    ebit = pbit - kbit
    n = len(positiveH)
    CurveParas = (Curve, p, pubA, pubB, G)
    if msb:
        # check form (x>> (pbit - kbit)) << (pbit - kbit)
        if any(h % 2**ebit != 0 for h in [H0] + positiveH + negativeH):
            # transform the msb leak to the form (x>> (pbit - kbit)) << (pbit - kbit)
            H0 = H0 << ebit
            positiveH = [h << ebit for h in positiveH]
            negativeH = [h << ebit for h in negativeH]
    Hs = [ZZ(ph) + ZZ(nh)  for ph, nh in zip(positiveH, negativeH)]
    polys = GetECHNPPolys(H0, Hs, xQ, CurveParas)
    if not msb:
        # lsb case
        x, y = polys[0].parent().gens()
        shift = 2 ** kbit
        # make the leading coeffient of x^2y being zero
        inv_shift = ZZ(inverse_mod(shift ** 3, p))
        polys = [poly(shift *x, shift * y) * inv_shift for poly in polys]
      
    # recentering optimization  
    ebit = ebit - 1
    bounds = [2**ebit] + [2**(ebit+1)] * n
    new_poly = []
    for poly in polys:
        pr = poly.parent()
        xs = list(pr.gens())
        shift_xs = {x: x + b for x,b in zip(xs, bounds)}
        new_poly.append(pr(poly.subs(shift_xs)))
    polys = new_poly[:]
    
    # generate the lattice polynomials
    if lattice == "XHS22":
        copper_polys = XHS22_Lattice_Polys(n, d, t, p, polys)
    else:
        copper_polys = XHS20_Lattice_Polys(n, d, p, polys)
    
    print(f"[+] Generating {len(copper_polys)} polynomials for coppersmith ({lattice})")
    reduced_polys = coppersmith_multivarivate_reduced_poly(copper_polys, bounds)
    # [TRIANGULATE, GROEBNER, JACOBIAN, HENSEL]
    if first_n == None:
        first_n = len(reduced_polys)
    print(f"[+] try the several methods to find roots using {first_n} of {len(reduced_polys)} polynomials")
    roots = rootfind_ZZ(reduced_polys[:first_n], bounds)
    if roots is None:
        return None
    for root in roots:
        print(f"[+] find : {root = }")
        # returns only one root
        if type(root) == list:
            e0 = ZZ(root[0]) + 2**ebit
        elif type(root) == dict:
            xs = reduced_polys[0].parent().gens()
            e0 = ZZ(root[xs[0]]) + 2**ebit
        else:
            assert False, "unknown root type"
        if msb:
            return ZZ(H0 + e0)
        else:
            return ZZ(H0 + (e0 << kbit))
    return None

def test_echnp_coppersmith_optimized_solver_secp256(kbit, n, d, t=1, lattice="XHS22", msb=True, first_n=None):
    E, G = secp256r1_curve()
    ecdh = ECDH(E, G)
    A, B = ecdh.get_pubs()
    shared_secret = ecdh.get_shared_secret()
    print(f"[+] The shared secret is {shared_secret}")
    p = ZZ(E.base_field().order())
    pbit = p.nbits()
    ebit = pbit - kbit
    print(f"[+] The basic parameters {pbit = } {kbit = } {msb = }")
    print(f"[+] Try to use lattice parameters {n = } {d = } {t = } with {lattice} Lattice")
    CurveParas = (E, p, A, B, G)
    oracle = ecdh.oracle_func(kbit, msb, True)
    h0, hs, xqs, phs, nhs = GenECHNP(oracle, CurveParas, n, kbit, msb)
    result = echnp_coppersmith_solver_optimized(E, G, A, B, kbit, h0, phs, nhs, xqs, 
                                                d, t, msb, lattice, first_n)
    if result is not None:
        print(f"[+] The shared secret : {result = }\n[+] check : {result == shared_secret}")        
    else:
        print("[+] The shared secret is not found")
        
def test_Coppersmallroot_With_Details_Secp256_XHS22_op(kbit = 172, n = 5, d = 2, t = 1, first_n = 10):
    E, G = secp256r1_curve()
    p = ZZ(E.base_field().order())
    pbit = p.nbits()
    ebit = pbit - kbit

    msb = True
    print(f"[+] The basic parameters {pbit = } {kbit = } {msb = }")
    print(f"[+] Try to use lattice parameters {n = } {d = } {t = } with XHS22 Lattice")
    polys, e, es = ECHNP_Polys(kbit, n)
    roots = [e] + es
    
    ebit = ebit - 1
    bounds = [2**ebit] + [2**(ebit+1)] * n
    new_poly = []
    for poly in polys:
        pr = poly.parent()
        xs = list(pr.gens())
        shift_xs = {x: x + b for x,b in zip(xs, bounds)}
        new_poly.append(pr(poly.subs(shift_xs)))
    # new roots
    for i, poly in enumerate(polys):
        assert poly(roots[0], roots[i+1]) % p == 0,  "Invalid polynomial for the original roots"
        
    roots = [root - bound for root, bound in zip(roots, bounds)]
    polys = new_poly[:]
    for i, poly in enumerate(polys):
        assert poly(roots[0], roots[i+1]) % p == 0,  "Invalid polynomial for the new roots"
    
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
    test_Coppersmallroot_With_Details_Secp256_XHS22_op(first_n = 40)
    # msb
    test_echnp_coppersmith_solver_groebner_secp256_op(kbit = 164, msb = True, n = 5, d = 3, t = 2, first_n = 80)
    # lsb
    test_echnp_coppersmith_solver_groebner_secp256_op(kbit = 164, msb = False, n = 5, d = 3, t = 2, first_n = 80)
    
    # test msb leak
    test_echnp_coppersmith_optimized_solver_secp256(kbit = 164, msb = True, n = 5, d = 3, t = 2, first_n = 20)
    # test lsb leak
    test_echnp_coppersmith_optimized_solver_secp256(kbit = 164, msb = True, n = 5, d = 3, t = 2, first_n = 20)