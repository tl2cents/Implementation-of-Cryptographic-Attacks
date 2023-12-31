from collections.abc import Iterable
from sage.modules.free_module_integer import IntegerLattice
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

# compute kernel space and obtain an LLL-reduced basis
def kernel_LLL(vec_mat, mod=None, verbose=False, w=None):
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

def derive_mod_bits(n):
    iota=0.035
    mod_bits=int(2 * iota * n^2 + n * log(n,2))
    return mod_bits


def ol_attack(h, m, n , M, verbose=False):
    """
    HSSP : h = X * a  % M
    Input :
    h : hssp result m-dim vector
    m,n : X.dimensions()
    M : the mod
    Output:
    the basis b1, ..., bn generating x1,...,xn (column vectors of X)
    """
    H = matrix(ZZ,h)
    # we only need m - n generating basis
    basis = kernel_LLL(H, mod = M, verbose=verbose)[:m - n]
    # check
    assert basis * H.T % M == 0, "not kernel, do check"
    # the basis is orthogonal to x_i in ZZ by assumption
    # try to recover the basis of xi
    xi_basis = kernel_LLL(basis, verbose=verbose)
    return xi_basis


def check_matrix(M, white_list = [1,0,-1]):
    # check wheter all values in M fall into white_list
    for row in M:
        for num in row:
            if num not in white_list:
                return False
    return True

def all_ones(v):
    if len([vj for vj in v if vj in [0,1]])==len(v):
        return v
    if len([vj for vj in v if vj in [0,-1]])==len(v):
        return -v
    return None
    
def recover_binary_basis(basis):
    lv = [all_ones(vi) for vi in basis if all_ones(vi)]
    n = basis.nrows()
    for v in lv:
        for i in range(n):
            nv = all_ones(basis[i] - v)
            if nv and nv not in lv:
                lv.append(nv)
            nv = all_ones(basis[i] + v)
            if nv and nv not in lv:
                lv.append(nv)
    return matrix(lv)

def find_original_basis(basis, new_basis):
    n, m = basis.dimensions()
    origin_lattice = IntegerLattice(basis)
    origin_basis = []
    for row in new_basis:
        if sum(row) == m:
            continue
        # seems like we cannot determine wether 1,-1 represents 1 or 0 in this lattcie
        # therefore we do some checking in the original lattice
        v = vector(ZZ, [1 if num == 1 else 0 for num in row])
        if v in origin_lattice:
            origin_basis.append(v)
        else:
            v = vector(ZZ, [0 if num == 1 else 1 for num in row])
            assert v in origin_lattice, "oops, something wrong"
            origin_basis.append(v)
    return matrix(ZZ, origin_basis)


def recover_binary_basis_by_lattice(basis, blocksize = None):
    new_lattice = 2 * basis
    n, m = basis.dimensions()
    new_lattice = new_lattice.insert_row(0, [1] * m)
    if blocksize is None: 
        # new_basis = new_lattice.LLL()
        new_basis = flatter(new_lattice)
    else:
        new_basis = new_lattice.BKZ(block_size = blocksize)
    
    if not check_matrix(new_basis, [1,-1]):
        print("[+] fails to recover basis")
        return None
        
    origin_lattice = IntegerLattice(basis)
    origin_basis = []
    for row in new_basis:
        if sum(row) == m:
            continue
        # seems like we cannot determine wether 1,-1 represents 1 or 0 in this lattcie
        # therefore we do some checking in the original lattice
        v = vector(ZZ, [1 if num == 1 else 0 for num in row])
        if v in origin_lattice:
            origin_basis.append(v)
        else:
            v = vector(ZZ, [0 if num == 1 else 1 for num in row])
            assert v in origin_lattice, "oops, something wrong"
            origin_basis.append(v)
    return matrix(ZZ, origin_basis)

# Nguyen-Stern attack using greedy method mentioned in appendix D of https://eprint.iacr.org/2020/461.pdf
def ns_attack_greedy(h, m, n, M, bkz = range(2,12,2), verbose=True):
    t0 = time.time()
    if verbose: print(f"[*] start to ns attack with greedy method")
    xi_basis = ol_attack(h, m, n, M)
    assert isinstance(bkz, Iterable), "give a list or iterable object as block_size para"
    L = xi_basis
    if verbose: print(f"    [+] basis dimensions : {L.dimensions()}")
    assert L.dimensions() == (n,m) , "basis generating xi's is not fully recovered"
    for bs in bkz:
        if verbose: print(f"    [*] start to BKZ reduction with block_size {bs}")
        L = L.BKZ(block_size = bs)
        if verbose: print(f"    [+] BKZ reduction with block_size {bs} done")
        if check_matrix(L,[-1,1,0]):
            if verbose: print("    [+] find valid basis")
            break
            
    XT = recover_binary_basis(L)
    if verbose: print(f"    [+] Number of recovered xi vectors {XT.nrows()}")
    if XT.nrows() < n:
        print(f"    [+] not enough xi vectors recovered, {XT.nrows()} out of {n}")
        print(f"    [*] trying new lattice recovery...")
        XT = recover_binary_basis_by_lattice(L)
        if XT.nrows() < n:
            print(f"[+] failed.")
            return False, L
    X =  XT.T
    # h = X * a
    a = matrix(Zmod(M) ,X[:n]).inverse() * vector(Zmod(M),h[:n])
    t1 = time.time()
    if verbose : print(f"    [+] total time cost in {t1 - t0}")
    return True, a

# Nguyen-Stern attack using 2*Lx + E lattice method
def ns_attack_2Lx(h, m, n, M, bkz = range(2,12,2), verbose=True):
    t0 = time.time()
    if verbose: print(f"[*] start to ns attack with 2*Lx + E method")
    xi_basis = ol_attack(h, m, n, M)
    assert isinstance(bkz, Iterable), "give a list or iterable object as block_size para"
    # we use the new lattice : 2 * basis + [1, ..., 1], the final vectors all fall into [-1, 1]
    L = 2 * xi_basis
    L = L.insert_row(0, [1] * m)
    if verbose: print(f"    [+] basis dimensions : {L.dimensions()}")
    assert L.dimensions() == (n + 1, m) , "basis generating xi's is not fully recovered"
    for bs in bkz:
        if verbose: print(f"    [*] start to BKZ reduction with block_size {bs}")
        L = L.BKZ(block_size = bs)
        if verbose: print(f"    [+] BKZ reduction with block_size {bs} done")
        if check_matrix(L,[-1,1]):
            if verbose: print("    [+] find valid basis")
            break
            
    XT = find_original_basis(xi_basis, L)
    if verbose: print(f"    [+] Number of recovered xi vectors {XT.nrows()}")
    if XT.nrows() < n:
        return False, L
    X =  XT.T
    # h = X * a
    a = matrix(Zmod(M) ,X[:n]).inverse() * vector(Zmod(M),h[:n])
    t1 = time.time()
    if verbose : print(f"    [+] total time cost in {t1 - t0}")
    return True, a

if __name__ == "__main__":
    n = 128
    m = 256
    Mbits = derive_mod_bits(n)
    M, a, X, h = gen_hssp(n, m, Mbits)
    find, recovered_a = ns_attack_greedy(h, m, n, M, range(2,32,4))
    if find:
        print(f"[+] check {sorted(a) == sorted(recovered_a)}")
        print()
    
    find, recovered_a = ns_attack_2Lx(h, m, n, M, range(2,32,4))
    if find:
        print(f"[+] check {sorted(a) == sorted(recovered_a)}")    
        print()