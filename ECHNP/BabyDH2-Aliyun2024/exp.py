import sys
# add .. path
sys.path.append('../')
from pwn import remote, process, context, log
from pwnlib.util.iters import mbruteforce
from Crypto.Cipher import AES
from ecdsa import NIST256p
from ecdsa.ecdsa import Public_key
from ecdsa.ellipticcurve import Point
from hashlib import sha256
from os import urandom
from random import seed, choice, randrange
from string import ascii_letters, ascii_lowercase, digits
from optimized_echnp_solver import echnp_coppersmith_solver_optimized
from sage.all import EllipticCurve, GF, ZZ

io = process(["python3", "task.py"])

def solve_pow(io:remote):
    # print('sha256(XXXX + {}) == {}'.format(proof[4: ], digest))
    # print('Give me XXXX:')
    io.recvuntil(b"sha256(XXXX + ")
    suffix = io.recvuntil(b") == ", drop=True)
    target_digest = io.recvuntil(b"\n", drop=True).decode()
    log.info(f"suffix: {suffix}, target_digest: {target_digest}")
    sol = mbruteforce(lambda x: sha256(x.encode() + suffix).hexdigest() == target_digest, ascii_letters + digits, 4)
    log.info(f"solved pow: {sol}")
    io.sendlineafter(b"Give me XXXX:\n", sol.encode())
    return True

def submit_pk(Qx, Qy):
    io.sendline(hex(Qx)[2:].zfill(64).encode())
    io.sendline(hex(Qy)[2:].zfill(64).encode())
    io.recvuntil(b'Leak: ')
    return int(io.recvline().strip().decode(), 16)

def ECHNP_Sample_OracleA(n, A, B, G):
    # oracleA returns x-coordinate of a * Q
    Q = 0 * G  + B
    H0 = ZZ(submit_pk(Q[0], Q[1]))
    positiveH = []
    negativeH = []
    xQ = []
    for i in range(1, n + 1):
        Q = i * G + B
        positiveH.append(ZZ(submit_pk(Q[0], Q[1])))
        Q = -i * G + B
        negativeH.append(ZZ(submit_pk(Q[0], Q[1])))
        xQ.append(ZZ((i * A)[0]))
    return H0, positiveH, negativeH, xQ
    
E = NIST256p.curve
p = E.p()
a = E.a()
b = E.b()
Sage_Curve = EllipticCurve(GF(p), [a, b])
log.info(f"Curve: {Sage_Curve}")
G = NIST256p.generator
m = 2**164
n = G.order()
R = GF(p)

solve_pow(io)

io.recvuntil(b'Hi, Bob! Here is my public key: ')
Ax = R(int(io.recvline().strip().decode(), 16))
io.recvuntil(b'Hi, Alice! Here is my public key: ')
Bx = R(int(io.recvline().strip().decode(), 16))
signed_secret = bytes.fromhex(io.recvline().strip().decode())

A = Sage_Curve.lift_x(Ax)
B = Sage_Curve.lift_x(Bx)
G = Sage_Curve.lift_x(R(G.x()))
log.info(f"A: {A}")
log.info(f"B: {B}")
log.info(f"G: {G}")
n = 5 
d = 3
t = 2
kbit = 164
# def echnp_coppersmith_solver_optimized(Curve, G, pubA, pubB, kbit, H0, positiveH, negativeH, xQ, d, t=1, 
                            #  msb=True, lattice="XHS22", first_n=None):
H0, positiveH, negativeH, xQ = ECHNP_Sample_OracleA(n, A, B, G)
shared_AB = echnp_coppersmith_solver_optimized(Sage_Curve, G, A, B, kbit, H0, positiveH, negativeH, xQ, d, t, False, "XHS22", 20)
log.info(f"shared_AB: {shared_AB}")

key = sha256(str(shared_AB).encode()).digest()
msg = AES.new(key, AES.MODE_ECB).decrypt(signed_secret)
# pr(f'Secret sign is {secret_sign}.', shared_AB)
real_secret = msg.decode()[len("Secret sign is "):-1]
assert len(real_secret) == 64, "Invalid secret length"
log.info(f"real_secret: {real_secret}")
key = sha256(str(shared_AB).encode()).digest()
payload = AES.new(key, AES.MODE_ECB).encrypt(real_secret.encode()).hex()
io.sendline(payload.encode())
print(io.recvline())
io.interactive()