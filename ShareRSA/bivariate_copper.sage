from math import ceil
from Crypto.Util.number import isPrime
from random import randrange, getrandbits

def gen_rsa_chall(delta, beta, nbits=1000):
    # alpha = 2, typically
    kbit = ceil(nbits * beta)
    ubit = (nbits - 2*kbit)//2
    while True:
        k = getrandbits(kbit)
        pp = getrandbits(ubit) << kbit
        qq = getrandbits(ubit) << kbit
        p = pp + k
        q = qq + k
        if isPrime(p) and isPrime(q):
            break
    if q > p:
        p, q = q, p

    n = p * q
    phi = (p ** 2 - 1) * (q ** 2 - 1)
    ub = int(n ** delta)
    lb = int(n ** (delta - 0.02))
    while True:
        d = randrange(lb, ub)
        if gcd(d, phi) == 1:
            break
    e = int(inverse_mod(d, phi))
    sk = (p, q, d)
    pk = (n, e)
    return pk, sk

# pk, sk = gen_rsa_chall(0.7, 0.1)

def delat_upper_bound(beta, alpha=2):
    # alpha = 2, typically
    return RR(7/3 - (4/3 * beta) - 2/3 * sqrt((1 - 4 * beta)*(1 + 3 * alpha - 4 * beta)))

def derive_2r_bits_leak(N, r):
    Zr = Zmod(2^r)
    u0_list = Zr(N).nth_root(2, all = True)
    vs = []
    for u0 in u0_list:
        u0 = ZZ(u0)
        v0 = 2*u0 + (N - u0^2) * ZZ(inverse_mod(u0, 2^(2 * r)))
        vs.append(v0 % 2^(2*r))
    return vs

def construct_polynomials(N, e, beta):
    nbits = ZZ(N).nbits()
    r = ceil(nbits * beta)
    v0_list = derive_2r_bits_leak(N, r)
    PR = PolynomialRing(Zmod(e), ["x", "y"],  order = "lex")
    x, y = PR.gens()
    polys = []
    for v0 in v0_list:
        a1 = v0 * ZZ(inverse_mod(2^(2 * r -1), e)) % e
        a2 = (v0^2 - (N+1)^2) * ZZ(inverse_mod(2^(4 * r), e)) % e
        a3 = ZZ(inverse_mod(-2^(4 * r), e))
        poly = x * y^2 + a1 * x * y + a2 * x + a3
        polys.append(poly)
    return v0_list, polys

def flatter(M):
    from subprocess import check_output
    from re import findall
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), map(int, findall(b"-?\\d+", ret)))

def G_sij(f, s, i, j, m, e):
    x, y = f.parent().gens()
    return x^(i - s) * y^(j - 2*s) * f^s * e^(m - s)

def gen_copper_polys(f, m, t, e):
    gpolys = []
    x, y = f.parent().gens()
    for s in range(0, m+1):
        for i in range(s, m+1):
            for j in range(2*s, 2*s+2):
                gpolys.append(G_sij(f(x, y), s, i, j, m, e))
    hpolys = []
    for s in range(0, m+1):
        for i in range(s, s + 1):
            for j in range(2*s + 2, 2*s + t + 1):
                hpolys.append(G_sij(f(x, y), s, i, j, m, e))
    return gpolys, hpolys

def bivarivate_small_roots(f, X, Y, m, t, poly_num=3):
    R = f.base_ring()
    e = R.cardinality()
    f = f.change_ring(ZZ)
    
    g_polys, hpolys = gen_copper_polys(f, m, t, e)
    
    G = Sequence(g_polys + hpolys, f.parent())
    B, monomials = G.coefficient_matrix()
    monomials = vector(monomials)
    factors = [monomial(X, Y) for monomial in monomials]
    for i, factor in enumerate(factors):
        B.rescale_col(i, factor)
    
    print(f"[+] start LLL with {B.dimensions() = }")
    B = flatter(B.dense_matrix())
    print("[+] LLL done")

    B = B.change_ring(QQ)
    for i, factor in enumerate(factors):
        B.rescale_col(i, 1/factor)
    polys = B * monomials
    selected_polys = polys[:poly_num]
    
    x, y = polys[0].parent().gens()
    roots = []
    for poly1 in selected_polys:
        for poly2 in selected_polys:
            if poly1 == poly2:
                continue
            poly_res = poly1.resultant(poly2, x)
            if poly_res.is_constant():
                continue
            poly_univar_y = poly_res.univariate_polynomial()
            y_roots = poly_univar_y.roots(ring=ZZ, multiplicities=False)
            if len(y_roots)!= 0:
                for y_root in y_roots:
                    if abs(y_root) >= Y:
                        print(f"[+] unbounded root find {y_root = }, please check")
                        continue
                    poly_univar_x = poly1(x, y_root).univariate_polynomial()
                    x_roots = poly_univar_x.roots(ring=ZZ, multiplicities=False)
                    for x_root in x_roots:
                        if abs(x_root) >= X:
                            print(f"[+] unbounded root find {x_root = }, please check")
                            continue
                        roots.append((x_root, y_root))
                # we will not check other polynomials
                return roots
    return roots
            
def check_paper_samples():
    N = 611402847577596838649117628567007514815745193613363898749361
    e = 256620536587836325389289742308993660982466972819568896656661249105081887104266414292000389611061562475252796784804753727
    beta = 0.1
    delta = 0.7
    r = 20

    alpha = RR(log(e, N))
    X = int(2 * N^(alpha + delta - 2))
    Y = int(3 * N^(0.5 - 2*beta))
    m = 4
    t = 4

    k = 17387477862024536259218971283032599828907
    v = 1433911212640302358

    polys = construct_polynomials(N, e, beta)
    for poly in polys:
        roots = bivarivate_small_roots(poly, X, Y, m, t)
        if len(roots) != 0:
            print(f"[+] recovered roots {roots = }")

if __name__ == "__main__":
    print(f"[+] beta = 0.0, delta = {delat_upper_bound(0)}")
    print(f"[+] beta = 0.1, delta = {delat_upper_bound(0.1)}")

    nbits = 1000
    beta =  0.1
    delta = 0.7
    r = ceil(nbits * beta)
    print("[*] generating samples")
    (N,e), (p,q,d) = gen_rsa_chall(delta, beta, nbits)
    print("[+] generation done")
    
    v0_list, polys = construct_polynomials(N, e, beta)
    
    alpha = RR(log(e, N))
    X = int(2 * N^(alpha + delta - 2))
    Y = int(3 * N^(0.5 - 2*beta))
    m = 4
    t = 4
    
    # k = (e*d - 1)//((p^2 -1)*(q^2-1))
    # v = ((p+q) - (p+q)%(2^(2*r)))>>(2*r)

    for poly, v0 in zip(polys, v0_list):
        roots = bivarivate_small_roots(poly, X, Y, m, t)
        if len(roots) != 0:
            print(f"[+] recovered roots {roots = }")
            k, v = roots[0]
            if v < 0:
                p_plus_q = (- v - 1) * 2^(2*r) + (2^(2*r) - v0)
            else:
                p_plus_q = v * 2^(2*r) + v0
            
            p_minus_q = ZZ(sqrt(p_plus_q**2 - 4 * N))
            rp = (p_plus_q + p_minus_q) // 2
            rq = (p_plus_q - p_minus_q) // 2
            assert rp * rq == N
            print(f"[+] successfully factored {N} = {p} * {q}")