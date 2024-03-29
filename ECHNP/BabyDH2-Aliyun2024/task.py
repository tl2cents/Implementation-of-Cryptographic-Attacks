from Crypto.Cipher import AES
from ecdsa import NIST256p
from ecdsa.ecdsa import Public_key
from ecdsa.ellipticcurve import Point
from hashlib import sha256
from os import urandom
from random import seed, choice, randrange
# from secret import FLAG
FLAG = 'flag{this_is_a_fake_flag_12abcd}'
from signal import signal, alarm, SIGALRM
from string import ascii_letters, ascii_lowercase, digits
from sys import stdin

def proof_of_work():
    seed(urandom(8))
    proof = ''.join([choice(ascii_letters + digits) for _ in range(20)])
    digest = sha256(proof.encode('latin-1')).hexdigest()
    print('sha256(XXXX + {}) == {}'.format(proof[4: ], digest))
    print('Give me XXXX:')
    x = read(4)
    if x != proof[: 4]:
        return False
    return True

def handler(signum, frame):
    print('Time out!')
    raise exit()

def read(l):
    return stdin.read(l + 1).strip()

def pr(msg, key=None):
    if not key:
        print(msg)
    else:
        key = sha256(str(key).encode()).digest()
        print(AES.new(key, AES.MODE_ECB).encrypt(msg.encode()).hex())

def inp():
    try:
        return Point(E, int(read(64), 16), int(read(64), 16))
    except:
        pass
    return None

def DH(priv, pub):
    shared = priv * pub
    return shared.x()

signal(SIGALRM, handler)
alarm(60)

if not proof_of_work():
    exit()

E = NIST256p.curve
G = NIST256p.generator
m = 2**164
n = G.order()

Alice_sk = randrange(n)
Alice_pk = Public_key(G, Alice_sk * G).point
pr(f'Hi, Bob! Here is my public key: {Alice_pk.x() :x}')

Bob_sk = randrange(n)
Bob_pk = Public_key(G, Bob_sk * G).point
pr(f'Hi, Alice! Here is my public key: {Bob_pk.x() :x}')

shared_AB = DH(Alice_sk, Bob_pk)
shared_BA = DH(Bob_sk, Alice_pk)
assert shared_AB == shared_BA

secret_sign = ''.join([choice(ascii_letters + digits) for _ in range(64)])
pr(f'Secret sign is {secret_sign}.', shared_AB)

alarm(20)
pr('Now, it is your turn:')
for _ in range(11):
    Mallory_pk = inp()
    if not Mallory_pk:
        pr('Invalid pk!')
        exit()
    shared_AC = DH(Alice_sk, Mallory_pk)
    pr(f'Leak: {shared_AC % m :x}')

ct = bytes.fromhex(read(128))
key = sha256(str(shared_AB).encode()).digest()
if AES.new(key, AES.MODE_ECB).decrypt(ct).decode() == secret_sign:
    pr(FLAG)
else:
    pr('Wrong!')
