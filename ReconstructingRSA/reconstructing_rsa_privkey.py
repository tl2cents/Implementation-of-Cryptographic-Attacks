from Crypto.Util.number import getPrime, inverse
from random import sample
from itertools import product
from sage.all import PolynomialRing, ZZ, GF, Zmod, var, solve_mod
import sys

sys.setrecursionlimit(3000)
print(f"[+] current recursion limit: {sys.getrecursionlimit()}")

def gen_rsa_priv_key(nbits, ebits=None):
    n = 1
    while n.bit_length() != nbits:
        p = getPrime(nbits//2)
        q = getPrime(nbits//2)
        n = p * q
    if ebits is None:
        e = 65537
    else:
        e = getPrime(ebits)
    d = inverse(e, (p-1)*(q-1))
    dp = d % (p-1)
    dq = d % (q-1)
    qinv = inverse(q, p)
    return (n, e, d, p, q, dp, dq, qinv)

def rsa_random_leak_model(priv_key, p_per, q_per, d_per, dp_per, dq_per):
    n, e, d, p, q, dp, dq, qinv = priv_key
    nbits = int(n).bit_length()
    L = list(range(nbits))
    L2 = list(range(nbits//2))
    leak_d_bits = sample(L, int(nbits * d_per))
    leak_p_bits = sample(L2, int(nbits//2 * p_per))
    leak_q_bits = sample(L2, int(nbits//2 * q_per))
    leak_dp_bits = sample(L2, int(nbits//2 * dp_per))
    leak_dq_bits = sample(L2, int(nbits//2 * dq_per))
    leak_d = {}
    for pos in leak_d_bits:
        leak_d[pos] = (d >> pos) & 1
    leak_p = {}
    for pos in leak_p_bits:
        leak_p[pos] = (p >> pos) & 1
    leak_q = {}
    for pos in leak_q_bits:
        leak_q[pos] = (q >> pos) & 1
    leak_dp = {}
    for pos in leak_dp_bits:
        leak_dp[pos] = (dp >> pos) & 1
    leak_dq = {}
    for pos in leak_dq_bits:
        leak_dq[pos] = (dq >> pos) & 1
    return (leak_d, leak_p, leak_q, leak_dp, leak_dq)


def generate_k_d_pair(n,e):
    # e * d = k(p-1)(q-1) + 1 = k(n - p - q + 1) + 1
    # Given n and small e, compute all k candidates and its corresponding d
    assert e <= 2**20, "e is too large, consider using C implementation"
    for k in range(1, e):
        d_ = (k * n + 1) // e
        yield k, d_
        
def compute_k_d_pair(n, e, leak_d, start_pos=None):
    # p, q must be balanced : p < q < 2p
    # then d_ = (k * n + 1) // e agrees with d on their floor[nbits/2] - 2 msbs, omit extra 2 bits
    if start_pos == None:
        start_pos = int(n).bit_length()//2 + 2
    leak_d_highbits = {pos for pos in leak_d if pos >= start_pos}
    for k, d in generate_k_d_pair(n, e):
        if all((d >> pos) & 1 == leak_d[pos] for pos in leak_d_highbits):
            yield k, d
    
def test_recovering_k(leak_percent = 0.3):
    n, e, d, p, q, dp, dq, qinv = gen_rsa_priv_key(1024)
    k = (e * d - 1) // ((p-1)*(q-1))
    leak_d, leak_p, leak_q, leak_dp, leak_dq = rsa_random_leak_model((n, e, d, p, q, dp, dq, qinv), leak_percent, leak_percent, leak_percent, leak_percent, leak_percent)
    pos = int(n).bit_length()//2 + 32
    L = list(compute_k_d_pair(n, e, leak_d))
    print(f"Recovered {len(L) = }")
    if len(L) == 0:
        print(f"Failed to recover k, try with pos = {pos}")
        L = list(compute_k_d_pair(n, e, leak_d, pos))
        print(f"Recovered {len(L) = }, {L = }")

        
def compute_kp_kq(n, e, k):
    # Given n, e and k, compute kp and kq such that
    # e * d = 1 + k(p-1)(q-1) = 1 + k * (n - p - q + 1)
    # e * dp = kp * (p-1) + 1
    # e * dq = kq * (q-1) + 1
    # n = p * q
    # solve x for eq: kp^2 - [k*(n-1) + 1] * kp - k = 0 \mod e
    # e should be small and prime (or product of several primes)
    x = var('x')
    sols = solve_mod([x**2 - (k*(n-1) + 1)*x - k == 0], e)
    return [int(sol[0]) for sol in sols]

def tau(x):
    r = 0
    while x % 2 == 0:
        x = x // 2
        r += 1
    return int(r)

def recover_info_from_rsa_priv_leaks(n, e, leaks):
    assert e % 2 == 1 and e > 2, "e must be odd and greater than 2"
    d_leak, p_leak, q_leak, dp_leak, dq_leak = leaks
    # copy leaks to avoid modifying the original leaks
    d_bits, p_bits, q_bits, dp_bits, dq_bits = d_leak.copy(), p_leak.copy(), q_leak.copy(), dp_leak.copy(), dq_leak.copy()
    nbits = int(n).bit_length()
    # Recover k
    if d_leak == None:
        return None
    
    pos = nbits // 2 + 2
    k_d = list(compute_k_d_pair(n, e, d_leak, pos))
    while len(k_d) == 0:
        pos += 5
        k_d = list(compute_k_d_pair(n, e, d_leak, pos))
        
    if len(k_d) > 1:
        print(f"[+] Warning: Multiple k candidates found: {len(k_d) = }")
    k, d_ = k_d[0]
    print(f"[+] Recovered k = {k}")
    kps = compute_kp_kq(n, e, k)
    # init all information we can derive from the known leaks
    for i in range(pos, nbits):
        if i in d_leak:
            assert (d_ >> i) & 1 == d_leak[i], "Inconsistent d leaks"
        else:
            d_bits[i] = (d_ >> i) & 1
    
    # lsb leaks of p, q, dp, dq, d
    p_bits[0] = 1
    q_bits[0] = 1
    kp_kq_s = product(kps, repeat=2)
    
    # lsb of d 
    tau_k = tau(k)
    d_lsb = inverse(e, 2**(2 + tau_k))
    for i in range(2 + tau_k):
        if i in d_leak:
            assert (d_lsb >> i) & 1 == d_leak[i], "Inconsistent d leaks"
        else:
            d_bits[i] = (d_lsb >> i) & 1
    
    for kp, kq in kp_kq_s:
        dp_bits = dp_leak.copy()
        dq_bits = dq_leak.copy()
        if kp == kq:
            # not considering kp == kq case
            continue
        tau_kp = tau(kp)
        tau_kq = tau(kq)
        dp_lsb = inverse(e, 2**(1 + tau_kp))
        dq_lsb = inverse(e, 2**(1 + tau_kq))
        good_guess = True
        for i in range(1 + tau_kp):
            if i in dp_leak:
                if (dp_lsb >> i) & 1 != dp_leak[i]:
                    # "Inconsistent dp leaks"
                    good_guess = False
                    break
            else:
                dp_bits[i] = (dp_lsb >> i) & 1
        if not good_guess:
            continue
        for i in range(1 + tau_kq):
            if i in dq_leak:
                if (dq_lsb >> i) & 1 != dq_leak[i]:
                    # "Inconsistent dq leaks"
                    good_guess = False
                    break
            else:
                dq_bits[i] = (dq_lsb >> i) & 1
        if good_guess:
            yield (n, e, k, kp, kq), (d_bits, p_bits, q_bits, dp_bits, dq_bits), (tau_k, tau_kp, tau_kq)
            
            
def bit_i(x, i):
    return (x >> i) & 1
            
def reconstructing_rsa_priv_key_v0(n, e, leaks):
    iter_table = list(product([0, 1], repeat=5))
    
    def sub_process(pk, known_bits, taus, p, q, d, dp, dq, pos):
        (n, e, k, kp, kq), (d_bits, p_bits, q_bits, dp_bits, dq_bits), (tau_k, tau_kp, tau_kq) = pk, known_bits, taus
        if pos == nbits//2 and (n % p == 0 or n % q == 0):
            print("[+] RSA private key recovered !")
            print(f"[+] {p = }")
            print(f"[+] {q = }")
            return (p, q, d, dp, dq)
        i = pos
        for pi, qi, dpi, dqi, di in iter_table:
            # the known bits            
            if i in p_bits and pi != p_bits[i]:
                continue
            if i in q_bits and qi != q_bits[i]:
                continue
            if (i + tau_kp) in dp_bits and dpi != dp_bits[i + tau_kp]:
                continue
            if (i + tau_kq) in dq_bits and dqi != dq_bits[i + tau_kq]:
                continue
            if (i + tau_k) in d_bits and di != d_bits[i + tau_k]:
                continue
            # the 4 equations
            if (pi + qi) % 2 != bit_i(n- p*q, i):
                continue
            if (di + pi + qi) % 2 != bit_i(k * (n + 1) + 1 - k *(p + q) - e * d, i + tau_k):
                continue
            if (dpi + pi) % 2 != bit_i(kp * (p - 1) + 1 - e * dp, i + tau_kp):
                continue
            if (dqi + qi) % 2 != bit_i(kq * (q - 1) + 1 - e * dq, i + tau_kq):
                continue
            # yield from sub_process(n,e,p,q,d,dp,dq,i+1)
            new_p = p | (pi << i)
            new_q = q | (qi << i)
            new_d = d | (di << (i + tau_k))
            new_dp = dp | (dpi << (i + tau_kp))
            new_dq = dq | (dqi << (i + tau_kq))
            result = sub_process(pk, known_bits, taus, new_p, new_q, new_d, new_dp, new_dq, i+1)
            if result is not None:
                return result  
        return None
                
    for pk, known_bits, taus in recover_info_from_rsa_priv_leaks(n, e, leaks):
        (n, e, k, kp, kq), (d_bits, p_bits, q_bits, dp_bits, dq_bits), (tau_k, tau_kp, tau_kq) = pk, known_bits, taus
        nbits = int(n).bit_length()
        print(f"[+] Try with {kp = }, {kq = }")
        p = 1
        q = 1
        d = ZZ([d_bits[i] for i in range(1 + tau_k)], base=2)
        dp = ZZ([dp_bits[i] for i in range(1 + tau_kp)], base=2)
        dq = ZZ([dq_bits[i] for i in range(1 + tau_kq)], base=2)
        res = sub_process(pk, known_bits, taus, p, q, d, dp, dq, 1)
        if res is not None:
            return res
        
def _reconstructing_rsa_priv_key(pk, known_bits, taus):
    iter_table = list(product([0, 1], repeat=5))
    def core_process(p, q, d, dp, dq, pos):
        if pos == max_pos and (n % p == 0 or n % q == 0):
            print("[+] RSA private key recovered !")
            print(f"[+] {p = }")
            print(f"[+] {q = }")
            return (p, q, d, dp, dq)
        i = pos
        for pi, qi, dpi, dqi, di in iter_table:
            # the known bits            
            if i in p_bits and pi != p_bits[i]:
                continue
            if i in q_bits and qi != q_bits[i]:
                continue
            if (i + tau_kp) in dp_bits and dpi != dp_bits[i + tau_kp]:
                continue
            if (i + tau_kq) in dq_bits and dqi != dq_bits[i + tau_kq]:
                continue
            if (i + tau_k) in d_bits and di != d_bits[i + tau_k]:
                continue
            # the 4 equations
            if (pi + qi) % 2 != bit_i(n- p*q, i):
                continue
            if (di + pi + qi) % 2 != bit_i(k * (n + 1) + 1 - k *(p + q) - e * d, i + tau_k):
                continue
            if (dpi + pi) % 2 != bit_i(kp * (p - 1) + 1 - e * dp, i + tau_kp):
                continue
            if (dqi + qi) % 2 != bit_i(kq * (q - 1) + 1 - e * dq, i + tau_kq):
                continue
            # yield from sub_process(n,e,p,q,d,dp,dq,i+1)
            new_p = p | (pi << i)
            new_q = q | (qi << i)
            new_d = d | (di << (i + tau_k))
            new_dp = dp | (dpi << (i + tau_kp))
            new_dq = dq | (dqi << (i + tau_kq))
            result = core_process(new_p, new_q, new_d, new_dp, new_dq, i+1)
            if result is not None:
                return result
        return None
    
    (n, e, k, kp, kq), (d_bits, p_bits, q_bits, dp_bits, dq_bits), (tau_k, tau_kp, tau_kq) = pk, known_bits, taus
    nbits = int(n).bit_length()
    max_pos = nbits//2
    print(f"[+] Try with {kp = }, {kq = }")
    p = 1
    q = 1
    d = int(ZZ([d_bits[i] for i in range(1 + tau_k)], base=2))
    dp = int(ZZ([dp_bits[i] for i in range(1 + tau_kp)], base=2))
    dq = int(ZZ([dq_bits[i] for i in range(1 + tau_kq)], base=2))
    return core_process(p, q, d, dp, dq, 1)

def _reconstructing_rsa_priv_key_iter(pk, known_bits, taus):
    from collections import deque 
    iter_table = list(product([0, 1], repeat=5))
    (n, e, k, kp, kq), (d_bits, p_bits, q_bits, dp_bits, dq_bits), (tau_k, tau_kp, tau_kq) = pk, known_bits, taus
    nbits = int(n).bit_length()
    max_pos = nbits//2
    print(f"[+] Try with {kp = }, {kq = }")
    p = 1
    q = 1
    d = int(ZZ([d_bits[i] for i in range(1 + tau_k)], base=2))
    dp = int(ZZ([dp_bits[i] for i in range(1 + tau_kp)], base=2))
    dq = int(ZZ([dq_bits[i] for i in range(1 + tau_kq)], base=2))
    stack = deque()
    stack.append((p, q, d, dp, dq, 1))
    # not empty
    while stack:
        p, q, d, dp, dq, pos = stack.pop()
        if pos == max_pos:
            if (n % p == 0 or n % q == 0):
                print("[+] RSA private key recovered !")
                print(f"[+] {p = }")
                print(f"[+] {q = }")
                return (p, q, d, dp, dq)
            else:
                continue
        i = pos
        for pi, qi, dpi, dqi, di in iter_table:
            # the known bits            
            if i in p_bits and pi != p_bits[i]:
                continue
            if i in q_bits and qi != q_bits[i]:
                continue
            if (i + tau_kp) in dp_bits and dpi != dp_bits[i + tau_kp]:
                continue
            if (i + tau_kq) in dq_bits and dqi != dq_bits[i + tau_kq]:
                continue
            if (i + tau_k) in d_bits and di != d_bits[i + tau_k]:
                continue
            # the 4 equations
            if (pi + qi) % 2 != bit_i(n- p*q, i):
                continue
            if (di + pi + qi) % 2 != bit_i(k * (n + 1) + 1 - k *(p + q) - e * d, i + tau_k):
                continue
            if (dpi + pi) % 2 != bit_i(kp * (p - 1) + 1 - e * dp, i + tau_kp):
                continue
            if (dqi + qi) % 2 != bit_i(kq * (q - 1) + 1 - e * dq, i + tau_kq):
                continue
            # yield from sub_process(n,e,p,q,d,dp,dq,i+1)
            new_p = p | (pi << i)
            new_q = q | (qi << i)
            new_d = d | (di << (i + tau_k))
            new_dp = dp | (dpi << (i + tau_kp))
            new_dq = dq | (dqi << (i + tau_kq))
            stack.append((new_p, new_q, new_d, new_dp, new_dq, i+1))
    return None

def _reconstructing_rsa_priv_key_iter_pq(pk, known_bits):
    # with only leaks of p and q
    from collections import deque 
    iter_table = list(product([0, 1], repeat=2))
    (n, e), (p_bits, q_bits) = pk, known_bits
    nbits = int(n).bit_length()
    max_pos = nbits//2
    p = 1
    q = 1
    stack = deque()
    stack.append((p, q, 1))
    # not empty
    while stack:
        p, q, pos = stack.pop()
        if pos == max_pos:
            if (n % p == 0 or n % q == 0):
                print("[+] RSA private key recovered !")
                print(f"[+] {p = }")
                print(f"[+] {q = }")
                return (p, q)
            else:
                continue
        i = pos
        for pi, qi in iter_table:
            # the known bits            
            if i in p_bits and pi != p_bits[i]:
                continue
            if i in q_bits and qi != q_bits[i]:
                continue
            # the equations
            if (pi + qi) % 2 != bit_i(n- p*q, i):
                continue
            new_p = p | (pi << i)
            new_q = q | (qi << i)
            stack.append((new_p, new_q, i+1))
    return None
        
def reconstructing_rsa_priv_key(n, e, leaks):            
    for pk, known_bits, taus in recover_info_from_rsa_priv_leaks(n, e, leaks):
        res = _reconstructing_rsa_priv_key(pk, known_bits, taus)
        if res is not None:
            return res
    return None

def reconstructing_rsa_priv_key_iter(n, e, leaks):        
    for pk, known_bits, taus in recover_info_from_rsa_priv_leaks(n, e, leaks):
        res = _reconstructing_rsa_priv_key_iter(pk, known_bits, taus)
        if res is not None:
            return res
    return None

def reconstructing_rsa_priv_key_iter_pq(n, e, leaks):        
    return _reconstructing_rsa_priv_key_iter_pq((n,e), leaks)

def test_reconstructing_rsa_priv_key():
    n, e, d, p, q, dp, dq, qinv = gen_rsa_priv_key(2048)
    print("[+] RSA Key information")
    print(f"[+] {p = }")
    print(f"[+] {q = }")
    k = (e * d - 1) // ((p-1)*(q-1))
    kp = (e * dp - 1) // (p-1)
    kq = (e * dq - 1) // (q-1)
    print(f"[+] {k = }")
    print(f"[+] {kp = }")
    print(f"[+] {kq = }")
    print()
    leaks = rsa_random_leak_model((n, e, d, p, q, dp, dq, qinv), 0.3, 0.3, 0.3, 0.3, 0.3)
    print("[+] Try iter version p q d dp dq leak model")
    reconstructing_rsa_priv_key_iter(n, e, leaks)
    print()

    print("[+] Try iter version p q d leak model")
    leaks = rsa_random_leak_model((n, e, d, p, q, dp, dq, qinv), 0.43, 0.43, 0.43, 0.0, 0.0)
    reconstructing_rsa_priv_key_iter(n, e, leaks)
    print()
    
    print("[+] Try iter version p q leak model")
    leaks = rsa_random_leak_model((n, e, d, p, q, dp, dq, qinv), 0.51, 0.51, 0.0, 0.0, 0.0)
    (leak_d, leak_p, leak_q, leak_dp, leak_dq) = leaks
    reconstructing_rsa_priv_key_iter_pq(n, e, (leak_p, leak_q))
    print()
    # print("Try recursing version")
    # priv = reconstructing_rsa_priv_key(n, e, leaks)

if __name__ == "__main__":
    test_reconstructing_rsa_priv_key()