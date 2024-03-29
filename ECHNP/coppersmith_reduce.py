from sage.all import matrix, ZZ, QQ, vector
from XHS_Lattice import fast_coef_mat

def flatter(M):
    from subprocess import check_output
    from re import findall
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), map(int, findall(b"-?\\d+", ret)))

def matrix_overview(mat):
    # 0 -> space, non-zero -> x
    for row in mat:
        for num in row:
            if num == 0:
                print("0", end="")
            else:
                print("x", end="")
        print()

def coppersmith_multivarivate_reduced_poly(polys, bounds):
    qq_poly_ring = polys[0].parent().change_ring(QQ)
    polys = [poly.change_ring(QQ) for poly in polys]
    mat, monomials = fast_coef_mat(polys)
    B = matrix(ZZ, mat)
    # matrix_overview(B)
    # print(f"[+] {B.is_triangular('lower') = }")
    # print(f"[+] {B.is_triangular('upper') = }")
    print(f"[+] {B.dimensions() = }. LLLing...")
    factors = [monomial(*bounds) for monomial in monomials]
    for i, factor in enumerate(factors):
        B.rescale_col(i, factor)
    try:
        # doing flatter, much faster than LLL
        B = flatter(B.dense_matrix())
        print("[+] Flatter-LLL Done")
    except:
        # native LLL reduction
        B = B.LLL()
        print("[+] Native-LLL Done")
    B = B.change_ring(QQ)
    for i, factor in enumerate(factors):
        B.rescale_col(i, 1/factor)
    monomials = vector(monomials)
    return list(filter(None, B*monomials))