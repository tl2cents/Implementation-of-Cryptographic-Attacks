from chall import alternating, direct, DirectMatrixProductCryptosystem
from random import Random, SystemRandom
from hashlib import sha256
from Crypto.Cipher import AES
from sage.all import load, save, matrix, ZZ, Zmod


# Break DirectMatrixProductCryptosystem of https://eprint.iacr.org/2023/1745.pdf


def gen_direct_challenge_local(paras=(10, 35, 2, 2**302 + 307)):
    rand = SystemRandom()
    cry = DirectMatrixProductCryptosystem(*paras)
    priv, pub = cry.keygen(rand)
    msg = cry.randmsg(rand)
    M = cry.encrypt(pub, msg)
    if cry.decrypt(priv, M) != msg:
        raise ValueError("Decryption failed")
    challenge = (pub, M)
    return challenge, msg

def dfs_search_message(C, pubkey):
    pubkey_inv = [A ** (-1) for A in pubkey]
    def dfs_search(current_c, current_path):
        # print(f"current_path: {current_path}")
        if len(current_path) == len(pubkey_inv) - 1:
            if current_c in pubkey:
                # print(f"[+] possible perm: {current_path}")
                yield current_path + [pubkey.index(current_c)]
        for i in range(len(pubkey_inv)):
            if i not in current_path:
                try_mat = pubkey_inv[i]
                if (try_mat * current_c).trace() <= current_c.trace():
                    yield from dfs_search(try_mat * current_c, current_path + [i])
    yield from dfs_search(C, [])

def break_direct_cryptosystem(paras, pubkey, C):
    """ Break the direct cryptosystem, decrypting the ciphertext C.

    Args:
        paras (tuple): the parameters of the cryptosystem i.e. (n, k, a, p)
        pubkey (list): list of matrix: bar A
        C (matrix): the ciphertext: C = prod(sigma, A)

    Returns:
        list: the decrypted permutation sigma
    """
    (n, k, a, p) = paras
    for perm in dfs_search_message(C, pubkey):
        if perm:
            print(f"Found perm {perm}")
            direct_cipher = DirectMatrixProductCryptosystem(n, k, a, p)
            return direct_cipher.decode(perm)
        
def test_break_direct_cryptosystem():
    (n, k, a, p) = 10, 35, 2, 2**302 + 307
    direct_paras = (n, k, a, p)
    (pubkey, C), msg = gen_direct_challenge_local(direct_paras)
    sol = break_direct_cryptosystem(direct_paras, pubkey, C)
    assert sol == msg, f"Failed: {sol} != {msg}"
    print("Passed test_break_direct_cryptosystem")

if __name__ == "__main__":
    test_break_direct_cryptosystem()