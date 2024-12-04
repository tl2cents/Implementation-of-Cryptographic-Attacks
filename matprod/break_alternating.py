from chall import AlternatingMatrixProductCryptosystem
from random import SystemRandom
from hashlib import sha256
from itertools import product
from Crypto.Cipher import AES
from sage.all import matrix, ZZ, Zmod, randint, prod, log, GF, block_matrix, vector
from tqdm import tqdm
import time

# Break AlternatingMatrixProductCryptosystem of https://eprint.iacr.org/2023/1745.pdf

# https://github.com/Neobeo/HackTM2023/blob/main/solve420.sage
def flatter(M):
    from subprocess import check_output
    from re import findall
    M = matrix(ZZ,M)
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = '[[' + ']\n['.join(' '.join(map(str,row)) for row in M) + ']]'
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), map(int,findall(b'-?\\d+', ret)))

def gen_alternating_challenge_local(flag=None, paras=(10, 64, 2, 2**553 + 549)):
    rand = SystemRandom()
    H = sha256()
    cry = AlternatingMatrixProductCryptosystem(*paras)
    priv, pub = cry.keygen(rand)
    msg = cry.randmsg(rand)
    M = cry.encrypt(pub, msg)
    if cry.decrypt(priv, M) != msg:
        raise ValueError("Decryption failed")
    H.update(str(msg).encode())
    challenge = (pub, M)

    if flag is None:
        flag = b"flag{" + rand.randbytes(16).hex().encode() + b"}"
    cipher = AES.new(H.digest(), AES.MODE_CTR)
    enc_flag = cipher.encrypt(flag)

    return {"challenge": challenge, "enc_flag": enc_flag, "nonce": cipher.nonce}, msg

def estimate_trace_bound(n, a, k):
    # https://eprint.iacr.org/2023/1745.pdf
    # The strict upper bound of the trace is too large in Corollary 1.
    # I don't find good formula to estimate the general trace bound, 
    # so I just use the average trace value of the random product matrix as the estimation.
    if k == 0:
        return n
    Ms = []
    for i in range(k):
        M = matrix(
            ZZ,
            n,
            n,
            [randint(0, a) for _ in range(n * n)],
        )
        Ms.append(M)
        
    pM = prod(Ms)
    # avg_pm = sum(pM.list()) // (n * n)
    # upper_bound = n * (a*n)**(k-1)
    return pM.trace() * 2


def gen_partial_ciphertext(pubkey, i, j, n_samples):
    """ Generate `n_samples` partial ciphertexts from index i to j including i.

    Args:
        pubkey (list): [A^0, A^1], two matrix sets.
        i (int): the start index
        j (int): the end index (exclusive)
        n_samples (int): the number of partial ciphertexts to generate
        
    Returns:
        list: the list of partial ciphertexts
    """
    # print(f"Generating partial ciphertexts from {i} to {j} with {n_samples} samples")
    Cs = []
    if 2**(j - i + 1) < n_samples:
        for num in range(2**(j - i)):
            bits = [int(bit) for bit in bin(num)[2:].zfill(j-i)]
            C = prod([pubkey[idx][bit] for idx, bit in enumerate(bits, start=i)]).list()
            Cs.append(C)
        return Cs
         
    while len(Cs) != n_samples:
        randbits = [randint(0, 1) for _ in range(j-i)]
        C = prod([pubkey[idx][bit] for idx, bit in enumerate(randbits, start=i)]).list()
        if C not in Cs:
            Cs.append(C)
    return Cs


def modulo_reduction(M, p, verbose=False):
    """ Perform LLL reduction on the matrix M with modulo p.

    Args:
        M (matrix): the matrix to reduce
        p (int): the modulo
        verbose (bool, optional): whether to print the debug information. Defaults to False.

    Returns:
        matrix: The reduced matrix
    """
    n, m = M.nrows(), M.ncols()
    if n < m:
        Me = M.change_ring(GF(p)).echelon_form()
        delta = Me.ncols() - n
        zero_mat = matrix.zero(delta, n)
        pI = matrix.identity(delta) * p
        L = block_matrix(ZZ, [[Me], [zero_mat.augment(pI)]])
        if L.rank() != L.nrows():
            L = L.echelon_form(algorithm="pari0", include_zero_rows=False, include_zero_columns=False)
        L = L.change_ring(ZZ)
    else:
        pI = matrix.identity(m) * p
        L = block_matrix(ZZ, [[M], [pI]])
        
    if verbose:
        st = time.time()
        print(f"Starts to do LLL(flatter) reduction with dimensions {L.dimensions()}")
    try:
        L = flatter(L)
    except Exception as e:
        print(f"Failed to use flatter: {e}")
        print(f"Starts to do sage built-in LLL reduction")
        L = L.LLL()
    if verbose: print(f"Ends LLL reduction in {time.time() - st:.2f}s")
    return L

def recover_Eij(paras, pubkey, i, j, n_sample = None):
    """ Recover the secret Ei*Ej^{-1} from `AlternatingMatrixProductCryptosystem`'s public key.

    Args:
        paras (tuple): the parameters of the cryptosystem i.e. (n, k, a, p)
        pubkey (list): [A^0, A^1], two matrix sets.
        n_sample (int): the number of partial ciphertexts to generate, if None, it will be set as estimated value.
    """
    assert i > j >= 0, f"Invalid {i = }, {j = }"
    (n, k, a, p) = paras
    # estimate the trace bound of Aj * A_{j+1} ... A_{i} wher A_i <= alpha
    trace_bound = int(estimate_trace_bound(n, a, i-j))
    # estimated_t = int(((ZZ(p).nbits() - 1) * n ** 2 / ZZ(p // trace_bound).nbits()))
    estimated_t = int(n**2 * log((p/2), p/trace_bound))
    
    # to make LLL algorithm work, we need n_samples > estimated_t
    if n_sample is None:
        n_sample = min(estimated_t + 32, estimated_t * 2)
    print(f"Number of samples used: {n_sample} (also the dimension of lattice)")
    print(f"Estimated bit-length of trace bound: {trace_bound.bit_length()}")
    partial_ciphertexts = gen_partial_ciphertext(pubkey, j, i, n_sample)
    
    # build lattice
    M = matrix(ZZ, partial_ciphertexts).T
    L = modulo_reduction(M, p, verbose=True)
    L = [v for v in L if v!=0]
    traces_vector= L[0]
    checks = [num < trace_bound for num in traces_vector]
    traces_vector_bits = [int(num).bit_length() for num in traces_vector]
    print(f"Average bit-length of recovered traces: {sum(traces_vector_bits)//len(traces_vector_bits)}")
    assert all(checks), "Failed to recover Eij"
    sol = M.change_ring(GF(p)).solve_left(traces_vector)
    return sol, trace_bound

def balanced_mod(x, p):
    x = ZZ(x) % p 
    return x if x <= p // 2 else x - p

def break_alternating_cryptosystem(paras, pubkey, C, step_size=16):
    """ Break the direct cryptosystem, decrypting the ciphertext C.

    Args:
        paras (tuple): the parameters of the cryptosystem i.e. (n, k, a, p)
        pubkey (list): list of matrix: bar A
        C (matrix): the ciphertext: C = prod(sigma, A)
        step_size (int): the step size of recovering bits, default is 16

    Returns:
        list: the decrypted message bits
    """
    assert step_size <= 24, "The step size is too large and may use a lot of memory and time"
    (n, k, a, p) = paras
    recovered_bits = []
    # pubkey_inv = [(A0 ** -1, A1 ** -1) for A0, A1 in pubkey]
    
    for i in range(0, k, step_size):
        print(f"Recovering bits from {i} to {i + step_size}")
        # the special case of the last bits
        if i + step_size >= k:
            # using direct brute-force to find the last bits
            step_size = k - i
            table1 = [(prod([pubkey[idx][bit] for idx, bit in enumerate(partial_sol, start=i)]) ** -1, partial_sol)
                  for partial_sol in product([0, 1], repeat=step_size//2)]
            # the second half of prod
            table2 = [(prod([pubkey[idx][bit] for idx, bit in enumerate(partial_sol, start=i + step_size//2)]) ** -1, partial_sol)
                    for partial_sol in product([0, 1], repeat=step_size - step_size//2)]
            for C1, partial_sol1 in tqdm(table1, desc="Trying to brute-force the valid partial solution"):
                for C2, partial_sol2 in table2:
                    partial_sol = partial_sol1 + partial_sol2
                    partial_ciphertext = C2 * C1 * C
                    if partial_ciphertext == 1:
                        print(f"Found partial solution of {i} to {i + step_size}")
                        recovered_bits.extend(partial_sol)
                        print(f"The recovered bits: {recovered_bits}")
                        return recovered_bits
            assert False, "Failed to find the last bits"
            
        # The general case
        Eki, trace_bound = recover_Eij(paras, pubkey, k, i + step_size)
        # find the solutions
        # bf_space = product([0, 1], repeat=step_size)
        # two step prod to speed up
        # the first half of prod
        table1 = [(prod([pubkey[idx][bit] for idx, bit in enumerate(partial_sol, start=i)]) ** -1, partial_sol)
                  for partial_sol in product([0, 1], repeat=step_size//2)]
        # the second half of prod
        table2 = [(prod([pubkey[idx][bit] for idx, bit in enumerate(partial_sol, start=i + step_size//2)]) ** -1, partial_sol)
                  for partial_sol in product([0, 1], repeat=step_size - step_size//2)]
        find_sol = False
        for C1, partial_sol1 in tqdm(table1, desc="Trying to brute-force the valid partial solution"):
            for C2, partial_sol2 in table2:
                partial_sol = partial_sol1 + partial_sol2
                partial_ciphertext = C2 * C1 * C
                trace_mA = abs(balanced_mod(Eki * vector(partial_ciphertext.list()), p))
                if trace_mA <= trace_bound:
                    print(f"Trace_mA : {int(trace_mA).bit_length()} bits, estimated trace_bound: {trace_bound.bit_length()} bits")
                    print(f"Found partial solution of {i} to {i + step_size}")
                    # right solution
                    C = partial_ciphertext
                    recovered_bits.extend(partial_sol)
                    print(f"Current recovered bits: {recovered_bits}")
                    find_sol = True
                    break
            if find_sol:
                break
        print()
    return recovered_bits

def test_break_alternating_cryptosystem():
    # a small size one
    # (n, k, a, p) = 10, 64, 2, 2**553 + 549
    (n, k, a, p) = 10, 128, 2, 2**553 + 549
    paras = (n, k, a, p)
    data, msg = gen_alternating_challenge_local(None, paras)
    pubkey, C = data["challenge"]
    # Eij = recover_Eij(paras, pubkey, 48, 0)
    m_bits = break_alternating_cryptosystem(paras, pubkey, C)
    msg_ = int("".join(map(str, m_bits[::-1])), 2)
    print(f"Original message: {msg}")
    print(f"Recovered message: {msg_}")
    assert msg == msg_, "Failed to recover the message"

if __name__ == "__main__":
    test_break_alternating_cryptosystem()