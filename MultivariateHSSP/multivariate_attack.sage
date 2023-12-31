from tqdm import tqdm
import time 

# https://github.com/Neobeo/HackTM2023/blob/main/solve420.sage
def flatter(M):
    from subprocess import check_output
    from re import findall
    M = matrix(ZZ,M)
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = '[[' + ']\n['.join(' '.join(map(str,row)) for row in M) + ']]'
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), map(int,findall(b'-?\\d+', ret)))

def gen_hssp(n = 10, m = 20, Mbits = 100):
    M = random_prime(2^Mbits, proof = False, lbound = 2^(Mbits - 1)) 
    # alpha vectors
    a = vector(ZZ, n)
    for i in range(n):
        a[i] = ZZ.random_element(M)
  
    # The matrix X has m rows and must be of rank n
    while True:
        X = random_matrix(GF(2), m, n).change_ring(ZZ)
        if X.rank() == n: break
    
    # Generate an instance of the HSSP: h = X*a % M
    h = X * a % M
    return M, a, X, h

def derive_mod_bits(n):
    iota=0.035
    mod_bits=int(2 * iota * n^2 + n * log(n,2))
    return mod_bits

# compute kernel space and obtain an LLL-reduced basis
def kernel_LLL(vec_mat, mod=None, verbose=True, w=None):
    """
    Input :
    vec_mat : m * n matrix, the m vectors are : v1,..., vm in Z^n
    mod     : if mod is not None, we find kernel in Zmod(mod) else in ZZ
    Output :
    B : matrix, an LLL-reduced basis b1,b2,...,bk such that bi * vj for all i in [1,k], j in [1,m]   
    """
    m, n = vec_mat.dimensions()
    if mod is None:
        # if n < 2*m : return vec_mat.right_kernel().matrix()
        if w is None : w=2^(n//2)* vec_mat.height()
        M = block_matrix(ZZ, [w * vec_mat.T, 1], ncols = 2).dense_matrix()
    else:
        if w is None : w = 2^(ZZ(mod).nbits())
        M = block_matrix(ZZ, [
            [w * vec_mat.T, 1],
            [w * mod, 0]
        ]).dense_matrix()
    if verbose: print(f"    [*] start to LLL reduction with dim {M.dimensions()}")
    # L = M.LLL()
    t0 = time.time()
    L = flatter(M)
    t1 = time.time()
    if verbose: print(f"    [+] LLL reduction done in {t1 -t0}")
    # filter orthogonal vectors
    basis = []
    for row in L:
        if row[:m] == 0:
            # orthogonal vector
            basis.append(row[m:m+n])
    if verbose: print(f"    [+] find {len(basis)} orthogonal vectors for {m} {n}-dim vectors")
    return matrix(ZZ, basis)


# implement the new multivariate attack with m//n LLL reductions
# algorithm 2 introduced in section 4.1 of https://eprint.iacr.org/2020/461.pdf 
# the paper implementations uses optimization mentioned in appendix F
def multivariate_ol_attack_step1(h, m, n, M, verbose=True):
    # return basis that can generate X, i.e. return basis_matrix of L(X)
    # m = d * m
    d = m // n
    h = list(h)
    hs = [h[i*n : (i+1)*n] for i in range(d)]
    Cs = []
    if verbose:
        iter_item = tqdm(range(d - 1))
    else:
        iter_item  = range(d - 1)
    for i in iter_item:
        yi = vector(ZZ, hs[0] + hs[i+1])
        C = ol_attack(yi, 2*n, n, M)[:n]
        if i == 0:
            C0 = C[:, :n]
            C1 = C[:, n:]
            C0_inv = C0
            Cs += [C0, C1]
        else:
            C0_i = C[:, :n]
            Ci = C[:, n:]
            Ci = C0 * C0_i^(-1) * Ci
            Cs.append(Ci)
    return block_matrix(ZZ, Cs, ncols = len(Cs))

def sub_kernel_LLL(h, M, verbose=False):
    # h = [h1 * h0_inv, h2 * h0_inv, ... , h2n* h0_inv]
    # return kernel space orthogonal to [h0, h1, ..., h2n] mod M
    t0 = time.time()
    n = len(h)
    L = identity_matrix(ZZ, n + 1)
    L[0,0] = M
    for i in range(n):
        L[i+1, 0] = h[i]
    L = flatter(L)
    # L = L.LLL()
    t1 = time.time()
    if verbose: print(f"    [+] sub kernel LLL with dim ({n+1, n+1}) costs {t1 - t0} seconds")
    return L

# implement kernel space computation with only one LLL reduction when m >> n
# algorithm 5 introduced in section 5.2 of https://eprint.iacr.org/2020/461.pdf 
def optimized_lll_reduced_orthogonal_vectors(h, m, n, M, verbose=False):
    t0 = time.time()
    if verbose: print(f"    [+] start to compute basis orthogonal to h")
    h0 = h[0]
    H1 = vector(ZZ, h[1:]) * ZZ(inverse_mod(-h0, M)) % M
    M_2n_basis = sub_kernel_LLL(H1.list()[:2 * n - 1], M, verbose)
    # THIS COSTS A LOT OF TIME !!!
    RF = RealField(2 * ZZ(M).nbits())
    M_2n_INV = matrix(RF, M_2n_basis).inverse()
    ortho_vecs = [list(vec) + [0] * (m - 2 * n) for vec in M_2n_basis[:n]]
    for i in range(m - 2 * n):
        idx = 2*n + i
        t = vector(ZZ, [H1[idx - 1]] + (2 * n -1) * [0])
        # Babai rounding process
        u = list(t[0] * M_2n_INV[0])
        u_vec = vector(ZZ, [round(num) for num in u])
        v_vec = vector(ZZ, u_vec * M_2n_basis)
        ai_2n = (t  - v_vec)
        ai_prime = list(ai_2n) + [1 if j==i else 0 for j in range(m - 2 * n)]
        ortho_vecs.append(ai_prime)
    t1 = time.time()
    if verbose: print(f"    [+] basis orthogonal to h computed in {t1 - t0} seconds")
    return matrix(ZZ, ortho_vecs)

# implement the new multivariate attack step 1 with only one LLL reduction
def optimized_multivariate_ol_attack_step1(h, m, n, M, verbose=False):
    # return basis that can generate X^T with size n * m, i.e. return basis_matrix of L(X^T)
    # return a matrix M with size n * m such that L(M) =  L(X^T)
    orthogonal_h_basis = optimized_lll_reduced_orthogonal_vectors(h, m, n, M, verbose=verbose)
    # orthogonal_h_basis is almost in Hermite Normal Form (after the first 2n components)
    # therefore we will use a straightforward and fast way to compute its kernel
    U = orthogonal_h_basis[:n, :2*n]
    V = orthogonal_h_basis[n:, :2*n]
    P = kernel_LLL(U, verbose=verbose).T
    C = block_matrix(ZZ, [P, - V * P], nrows = 2)
    assert orthogonal_h_basis * C == 0, "not orthogonal"
    return C.T

# implement the new multivariate attack step 1 with only one LLL reduction
def optimized_multivariate_ol_attack_step1_GF3(h, m, n, M, verbose=False):
    # return basis that can generate X^T with size n * m, i.e. return basis_matrix of L(X^T)
    # return a matrix M with size n * m such that L(M) =  L(X^T)
    orthogonal_h_basis = optimized_lll_reduced_orthogonal_vectors(h, m, n, M, verbose=verbose)
    # orthogonal_h_basis is almost in Hermite Normal Form (after the first 2n components)
    # therefore we will use a straightforward and fast way to compute its kernel
    U = orthogonal_h_basis[:n, :2*n]
    V = orthogonal_h_basis[n:, :2*n]
    P = kernel_LLL(U, verbose=verbose).T
    P3 = P.change_ring(GF(3))
    V3 = V.change_ring(GF(3))
    C = block_matrix(GF(3), [P3, - V3 * P3], nrows = 2)
    return C.T

# algorithm 3 introduced in section 4.2 of https://eprint.iacr.org/2020/461.pdf 
def _multivariate_attack(C, m, n, M, verbose=False):
    """
    input
    `C`: a given basis with size n * m of x1,...,xn
    return : x1,x2,...,xn in {0,1}^m  , such that L(x1,...,xn) = L(C) i.e. wi C = xi for i = 1,..,n
    """
    # turn matrix into GF(3)
    t1 = time.time()
    C = matrix(GF(3), C)
    mod = M
    # linearization matrix R: rj = ((2-delta_{i,k}) Cij Ckj) 1<=i<=k<=n
    R = matrix(GF(3), m, ZZ((n**2 + n)/2))
    look_up_dict = {}
    for j in range(m):
        idx = 0
        for i in range(0, n):
            for k in range(i,n):
                delta_ik = int(i==k)
                R[j, idx] = (2 - delta_ik) * C[i,j] * C[k, j]
                look_up_dict[tuple((i,k))] = idx 
                idx += 1
    E = block_matrix(GF(3), [R.T, - C], nrows = 2)
    ker = E.left_kernel().matrix()
    t2 = time.time()
    if verbose: print(f"    [+] {ker.dimensions() = } , R matrix construction and kernel computation cost {t2-t1} seconds")
    assert ker.dimensions() == (n, ZZ((n**2 + 3 * n)/2)), "weird dimensions"
    # turn kernel matrix into systematic form
    K1 = ker[:, :-n]
    K2 = ker[:, -n:]
    M = (K2^-1) * K1
    # duplicating columns of M
    M_prime =  matrix(GF(3),n, n**2)
    for i in range(n):
        for k in range(n):
            if i <= k:
                idx = look_up_dict[tuple((i,k))] 
            else:
                idx = look_up_dict[tuple((k,i))]
            M_prime[:,i*n + k] = M[:,idx]
    t3 = time.time()
    if verbose: print(f"    [+] systematic matrix and M_prime matrix constructed in {t3 -t2} seconds.")
    # We will compute the list of eigenvectors
    # We start with the full space.
    # We loop over the coordinates. This will split the eigenspaces.
    L1 = [Matrix(GF(3),identity_matrix(n))]    
    for j in range(n):               # We loop over the coordinates of the wi vectors.
        M = M_prime[:,j*n:(j+1)*n]   # We select the submatrix corresponding to coordinate j
        L2 = []                      # We initialize the next list
        for v in L1:
            if v.nrows()==1:         # We are done with this eigenvector 
                L2.append(v)
            else:                    # eigenspace of dimension >1
                A = v.solve_left(v * M)
                # v*M=A*v. When we apply M on the right, this is equivalent to applying the matrix A.
                # The eigenvalues of matrix A correspond to the jth coordinates of the wi vectors in that eigenspace
                for e, v2 in A.eigenspaces_left():    # We split the eigenspace according to the eigenvalues of A.
                    v2_mat = v2.matrix()
                    L2.append(v2_mat * v)                   # The new eigenspaces 
        L1 = L2
    t4 = time.time()
    if verbose: print(f"    [+] eigenvectors fully recovered in {t4 -t3} seconds.")
    XT = matrix([v[0] for v in L1]) * C
    for i in range(n):
        if any(c==2 for c in XT[i]): XT[i] = -XT[i]
    if verbose: print(f"    [+] Number of recovered xi vectors {XT.nrows()}")
    if XT.nrows() < n:
        return False, L
    X =  XT.T
    # h = X * a
    a = matrix(Zmod(mod) ,X[:n]).inverse() * vector(Zmod(mod),h[:n])
    return True, a

def multivariate_attack(h, m, n, M, verbose=True):
    # C = optimized_multivariate_ol_attack_step1(h, m, n, M, verbose)
    if verbose: print("[+] Multivariate attack for HSSP.")
    if verbose: print("[*] Step 1, recovering X basis...")
    t0 = time.time()
    C = optimized_multivariate_ol_attack_step1_GF3(h, m, n, M, verbose)
    t1 = time.time()
    if verbose: print(f"    [+] Total Step 1 : {t1-t0} seconds")
    if verbose: print(f"[*] Step 2 multivariate equation solving...")
    a = _multivariate_attack(C, m, n, M, verbose)
    t2 = time.time()
    if verbose: print(f"    [+] Total Step 2 : {t2-t1} seconds")
    if verbose: print(f"[+] Total Attack {t2-t0} seconds")
    return a

if __name__ == "__main__":
    n = 190

    if n % 2==1:
        m = n * (n + 3) // 2 # n is odd
    else:
        m = n * (n + 4) // 2 # n is even

    Mbits = derive_mod_bits(n)
    M, a, X, h  = gen_hssp(n, m, Mbits)

    if n <= 32:
        C = multivariate_ol_attack_step1(h, m, n, M)
        XT_lattice = IntegerLattice(X.T)
        for row in C:
            assert (row in XT_lattice)
        print("ok with multivariate_ol_attack_step1")
        
        # C1 = _optimized_multivariate_ol_attack_step1(h, m, n, M)
        C1 = optimized_multivariate_ol_attack_step1(h, m, n, M)

        XT_lattice = IntegerLattice(X.T)
        for row in C1:
            assert (row in XT_lattice)
        print("ok with optimized_multivariate_ol_attack_step1")
        
    find, recovered_a = multivariate_attack(h, m, n, M)
    if find:
        print(f"[+] check {sorted(a) == sorted(recovered_a)}")