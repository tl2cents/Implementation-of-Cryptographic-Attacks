# from findRootsZZ import find_roots_gcd, find_roots_groebner, find_roots_resultants, find_roots_variety
# from rootfind_ZZ import rootfind_ZZ, JACOBIAN, HENSEL, TRIANGULATE, GROEBNER 
from sage.all import prod, Sequence, matrix, vector, PolynomialRing, GF, ZZ, QQ, PolynomialRing, Zmod

def XHS20_Lattice_l_i(d, l, i0, js, x0, ys, polys, W=None):
    # assert 0 <= l <= d, "Invalid l"
    # assert 0 <= i0 <= 2*d, "Invalid i0"
    # assert len(js) == l, "Invalid positions"
    # case a
    if l == 0:
        return x0**i0
    # case b, l = 1
    elif l == 1 and i0 <= 1:
        return x0**i0 * ys[js[0]]
    # case c
    elif l >= 1 and i0 >= 2*l:
        F = prod(polys[j] for j in js)
        return x0**(i0 - 2*l) * F
    # case d
    elif l >= 2 and i0 <= 2*l -1:
        res =  0
        F = prod(polys[j] for j in js)
        for u in range(l):
            # F = prod(polys[j] for j in js if j != js[u])
            P = F // polys[js[u]]
            for v in range(2):
                res += W[i0, u + l * v] * (x0**v) * (P) * ys[js[u]]
        return res
    else:
        assert False,  f"Unknown Case {d = } { l = } {i0 = }"
        
def XHS22_Lattice_l_i_t(d, l, i0, t, js, x0, ys, polys, W=None, C=None, B=None, E=None):
    if i0 <= 2*d - 1 and l <= d:
        # case a : the same as XHS20 Lattice
        return XHS20_Lattice_l_i(d, l, i0, js, x0, ys, polys, W)
    elif i0 <= t and l == d + 1:
        res =  0
        F = prod(polys[j] for j in js)
        for u in range(l):
            # F = prod(polys[j] for j in js if j != js[u])
            P = F // polys[js[u]]
            for v in range(2):
                # H poly items
                res += W[i0, u + l * v] * (x0**v) * (P) * ys[js[u]]
                # J poly items
                res += W[i0, u + l * v] * (x0**v) * (P) * C[js[u]]
            # K poly items
            res += W[i0, u + l] * (P) * (B[js[u]] - C[js[u]] * E[js[u]])
        return res
    else:
        print(f"[+] unhanled case {d = } {l = } {i0 = } {t = }")
        return None

    
def M_polys(xpolys, positions):
    S = prod(xpolys[pos] for pos in positions)
    x = xpolys[0].parent().gens()[0]
    polylist1 = []
    polylist2 = []
    for i in range(len(positions)):
        tmp = S // xpolys[positions[i]]
        assert tmp * xpolys[positions[i]] == S
        polylist1.append(tmp)
        polylist2.append(x * tmp)
    return polylist1 + polylist2
    
def fast_coef_mat(polys):
    # must define the order of the monomials
    # mat = matrix(GF(2), len(polys), len(monos))
    monos = set()
    for poly in polys:
        monos.update(poly.monomials())
    monos = sorted(list(monos))
    # print(monos)
    mono_to_index = {}
    for i, mono in enumerate(monos):
        mono_to_index[mono] = i
    mat = [[0] * len(monos) for i in range(len(polys))]
    
    for i, f in (list(enumerate(polys))):
        for coeff, mono in f:
            # mat[i,mono_to_index[mono]] = 1
            mat[i][mono_to_index[mono]] = coeff
    return mat, monos

def _coeff_mat_uni(mono_to_index, monos, polys):
    # mat = matrix(GF(2), len(polys), len(monos))
    mat = [[0] * len(monos) for i in range(len(polys))]
    pr = polys[0].parent()
    x = pr.gens()[0]
    for i, f in (list(enumerate(polys))):
        for coeff in f:
            # mat[i,mono_to_index[mono]] = 1
            mat[i][mono_to_index[pr(x**i)]] = coeff
    return mat

def coeff_mat_uni(polys, maxdeg=None):
    if maxdeg is None:
        maxdeg = max([poly.degree() for poly in polys])
    mat = [[0] * (maxdeg + 1) for i in range(len(polys))]
    pr = polys[0].parent()
    x = pr.gens()[0]
    seq = Sequence(polys, pr)
    for i, f in (list(enumerate(polys))):
        for j,coeff in enumerate(f):
            mat[i][j] = coeff
    return mat

def XHS20_Lattice_Polys(n, d, p, polys, remove_d=False):
    # remove_d : for XHS22 lattice
    assert d <= n, "Invalid d"
    from itertools import combinations
    poly_ring = PolynomialRing(Zmod(p**d), [f"x0"] + [f"y{i}" for i in range(1, n+1)], order='negdegrevlex')
    poly_ring_uni = PolynomialRing(Zmod(p**d),"x")
    # polys = [ poly.change_ring(Zmod(p**d))  for poly in poly]
    xs = poly_ring.gens()
    x0, ys = xs[0], xs[1:]
    x = poly_ring_uni.gens()[0]
    # change polys's variables to x, y1, y2, ..., y_n
    spolys = []
    for poly, y in zip(polys, ys):
        spolys.append(poly_ring(poly(x0, y)))
    xpolys = []
    # save only y_i(x**2 + Ex + D) from poly = A + Bx + Cx**2 + Dy + Exy + x**2y
    for poly in polys:
        xx, yy = poly.parent().gens()
        D = ZZ(poly.coefficient({xx: 0, yy: 1}))
        E = ZZ(poly.coefficient({xx: 1, yy: 1}))
        xpolys.append(poly_ring_uni(D + E*x + x**2))
    monos = [poly_ring_uni(x**i) for i in range(2*d)]    
    res  = []
    Ws = {}
    i0_bound = 2 * d if remove_d else 2*d + 1
    for i0 in range(0, i0_bound):
        for l in range(0, d + 1):
            for i_pos in combinations(list(range(n)), l):
                if 2 <= l and i0 <= 2*l -1:
                    if str(i_pos) not in Ws:
                        mpolys = M_polys(xpolys, i_pos)
                        M = matrix(Zmod(p**d), coeff_mat_uni(mpolys))
                        # M = fast_coef_mat_uni(mono_to_index, monos, mpolys)
                        xmonos = vector(monos[:2*l])
                        xxpolys = M * xmonos
                        for p1, p2 in zip(mpolys, xxpolys):
                            assert p1 == p2
                        Ws[str(i_pos)] = M**(-1)
                    tmp_poly = XHS20_Lattice_l_i(d, l, i0, i_pos, x0, ys, spolys, Ws[str(i_pos)])
                else:
                    tmp_poly = XHS20_Lattice_l_i(d, l, i0, i_pos, x0, ys, spolys)
                if 1 <= l <= d and  0 <= i0 <= 2*l - 1:
                    res.append(p**(d + 1 - l) * tmp_poly.change_ring(ZZ))
                elif i0 >= 2*l:
                    res.append(p**(d - l) * tmp_poly.change_ring(ZZ))
    return res

def _XHS22_Lattice_Polys(n, d, t, p, polys):
    assert d <= n, "Invalid d"
    assert t <= 2*d - 1
    from itertools import combinations
    poly_ring = PolynomialRing(Zmod(p**d), [f"x0"] + [f"y{i}" for i in range(1, n+1)], order='negdegrevlex')
    poly_ring_uni = PolynomialRing(Zmod(p**d),"x")
    xs = poly_ring.gens()
    x0, ys = xs[0], xs[1:]
    x = poly_ring_uni.gens()[0]
    
    # change polys's variables to x, y1, y2, ..., y_n
    spolys = []
    for poly, y in zip(polys, ys):
        spolys.append(poly_ring(poly(x0, y)))
        
    # save only y_i(x**2 + Ex + D) from poly = A + Bx + Cx**2 + Dy + Exy + x**2y
    xpolys = []
    Bs = []
    Cs = []
    Es = []
    for poly in polys:
        xx, yy = poly.parent().gens()
        B = ZZ(poly.coefficient({xx: 1, yy: 0}))
        C = ZZ(poly.coefficient({xx: 2, yy: 0}))
        D = ZZ(poly.coefficient({xx: 0, yy: 1}))
        E = ZZ(poly.coefficient({xx: 1, yy: 1}))
        Bs.append(B)
        Cs.append(C)
        Es.append(E)
        xpolys.append(poly_ring_uni(D + E*x + x**2))
    # monos = [poly_ring_uni(x**i) for i in range(2*d)]    
    res  = []
    Ws = {}
    for i0 in range(0, 2*d):
        for l in range(0, d + 2):
            for i_pos in combinations(list(range(n)), l):
                if (2 <= l <= d and i0 <= 2*l -1):
                    if str(i_pos) not in Ws:
                        mpolys = M_polys(xpolys, i_pos)
                        M = matrix(Zmod(p**d), coeff_mat_uni(mpolys))
                        Ws[str(i_pos)] = M**(-1)
                    tmp_poly = XHS22_Lattice_l_i_t(d, l, i0, t, i_pos, x0, ys, spolys, Ws[str(i_pos)])
                elif i0 <= t and l == d + 1:
                    if str(i_pos) not in Ws:
                        mpolys = M_polys(xpolys, i_pos)
                        M = matrix(Zmod(p**d), coeff_mat_uni(mpolys))
                        Ws[str(i_pos)] = M**(-1)
                    tmp_poly = XHS22_Lattice_l_i_t(d, l, i0, t, i_pos, x0, ys, spolys, Ws[str(i_pos)], Cs, Bs, Es)
                else:
                    tmp_poly = XHS22_Lattice_l_i_t(d, l, i0, t, i_pos, x0, ys, spolys)
                if tmp_poly == None:
                    break
                if l == d + 1:
                    res.append(tmp_poly.change_ring(ZZ))
                if 1 <= l <= d and  0 <= i0 <= 2*l - 1:
                    res.append(p**(d + 1 - l) * tmp_poly.change_ring(ZZ))
                elif i0 >= 2*l:
                    res.append(p**(d - l) * tmp_poly.change_ring(ZZ))
    return res    

def XHS22_Lattice_Polys(n, d, t, p, polys):
    assert d <= n, "Invalid d"
    assert t <= 2*d - 1
    from itertools import combinations
    poly_ring = PolynomialRing(Zmod(p**d), [f"x0"] + [f"y{i}" for i in range(1, n+1)], order='negdegrevlex')
    poly_ring_uni = PolynomialRing(Zmod(p**d),"x")
    xs = poly_ring.gens()
    x0, ys = xs[0], xs[1:]
    x = poly_ring_uni.gens()[0]
    
    # change polys's variables to x, y1, y2, ..., y_n
    spolys = []
    for poly, y in zip(polys, ys):
        spolys.append(poly_ring(poly(x0, y)))
        
    # save only y_i(x**2 + Ex + D) from poly = A + Bx + Cx**2 + Dy + Exy + x**2y
    xpolys = []
    Bs = []
    Cs = []
    Es = []
    for poly in polys:
        xx, yy = poly.parent().gens()
        B = ZZ(poly.coefficient({xx: 1, yy: 0}))
        C = ZZ(poly.coefficient({xx: 2, yy: 0}))
        D = ZZ(poly.coefficient({xx: 0, yy: 1}))
        E = ZZ(poly.coefficient({xx: 1, yy: 1}))
        Bs.append(B)
        Cs.append(C)
        Es.append(E)
        xpolys.append(poly_ring_uni(D + E*x + x**2))
    # monos = [poly_ring_uni(x**i) for i in range(2*d)]    
    res  = []
    Ws = {}
    res += XHS20_Lattice_Polys(n, d, p, polys, remove_d=True)
    for i0 in range(0, t + 1):
        l = d + 1
        for i_pos in combinations(list(range(n)), l):
            if str(i_pos) not in Ws:
                mpolys = M_polys(xpolys, i_pos)
                M = matrix(Zmod(p**d), coeff_mat_uni(mpolys))
                Ws[str(i_pos)] = M**(-1)
            tmp_poly = XHS22_Lattice_l_i_t(d, l, i0, t, i_pos, x0, ys, spolys, Ws[str(i_pos)], Cs, Bs, Es)
            if tmp_poly == None:
                break
            res.append(tmp_poly.change_ring(ZZ))
    return res
    