# Implementation of Fast Syndrome Based Hash
from sage.all import matrix, random_matrix, GF, ZZ, vector, block_matrix
from secrets import randbits

from tqdm import trange, tqdm
secure_paras = {
    # security level : [(r,w,n,n/w), ...]
    "64": [(512, 512, 131072, 256),
           (512, 450, 230400, 512),
           (1024, 2**17, 2**25, 256)],
    "80": [(512, 170, 43520, 256),
           (512, 144, 73728, 512)],
    "128": [(1024, 1024, 262144, 256),
            (1024, 904, 462848, 512),
            (1024, 816, 835584, 1024)]
}

class Fast_Syndrome_Based_Hash:
    
    def __init__(self, r, w, n, word_size=None, random_init=True) -> None:
        # FSB H matrix dimension : r x n
        # n//w is usually 2^8 
        # the hash result is r bits, single message size is n bits
        # total w sub-matrixes
        assert n % w == 0, "n must be divisible by w"
        self.r = r
        self.w = w
        self.n = n
        if word_size != None:
            assert word_size == n//w, "Word size does not match"
        self.word_size = n//w
        if random_init:
            self.init()
        
    def __generate_matrix(self):
        self.Hcols = [randbits(self.r) for i in range(self.n)]
        # self.H = matrix(GF(2), [ZZ(num).digits(base = 2, padto = self.r) for num in self.Hcols]).T
        self.H = None
    
    def __set_matrix(self, mat: matrix):
        assert mat.dimensions() == (self.r, self.n), "Matrix dimensions do not match"
        if mat.base_ring() != GF(2):
            mat = mat.change_ring(GF(2))
        self.H = mat
        
    def init(self, H=None):
        if H == None:
            self.__generate_matrix()
        else:
            self.__set_matrix(H)
            # fast computation
            self.Hcols = [ZZ(self.H.column(i).list(), base=2) for i in trange(self.n)]
        self.initialized = True
    
    def msg2list(self, m:bytes):
        assert len(m) * 8 <= self.n, "Message length does not match"
        num = int.from_bytes(m, byteorder='big')
        vec = ZZ(num).digits(self.word_size)
        padded = [0]*(self.w - len(vec))
        return vec + padded
    
    def __hash_core(self, m):
        # deprecated
        assert self.initialized, "FSB not initialized, call init() first"
        ss = self.msg2list(m)
        val = sum([self.H.column(i*self.word_size + s) for i,s in enumerate(ss)])
        return val.list()
    
    def fast_hash_core(self, m):
        assert self.initialized, "FSB not initialized, call init() first"
        ss = self.msg2list(m)
        val = 0
        for i,s in enumerate(ss):
            val ^= self.Hcols[i*self.word_size + s]
        return val
    
    def __digest(self, m:bytes):
        # deprecated
        res = ZZ(self.__hash_core(m), base=2)
        return int(res).to_bytes(self.r//8, byteorder='big')
    
    def digest(self, m:bytes):
        res = self.fast_hash_core(m)
        return int(res).to_bytes(self.r//8, byteorder='big')
    
    def __hexdigest(self, m:bytes):
        # deprecated
        return self.__digest(m).hex()
    
    def hexdigest(self, m:bytes):
        return self.digest(m).hex()
    
def fsb_linearization(Hcols, r, w, n, candidates=[0,1]):
    # Return FSB(x) = Ax + c in GF(2)^r where x_i in candidates
    word_size = n//w
    c = 0
    As = []
    for i in range(w):
        c ^= Hcols[i*word_size + candidates[0]]
        As.append(ZZ(Hcols[i*word_size + candidates[1]] ^ Hcols[i*word_size + candidates[0]]).digits(base=2, padto=r))
    A = matrix(GF(2), As).T
    return A, c
    
def test_FSB():
    r = 256
    w = 64
    n = 256 * w
    fsb = Fast_Syndrome_Based_Hash(r, w, n)
    fsb.init()
    m = b"Hello, World!"
    print(fsb.digest(m))
    print(fsb.hexdigest(m))

def test_fsb_pre_image():
    print("[+] Pre-image attack on FSB64")
    r, w, n, word_size = secure_paras["64"][0]
    FSB64 = Fast_Syndrome_Based_Hash(*secure_paras["64"][0])
    ascii_table = [i for i in range(65, 0x7f)]
    full_rank = False
    for c1 in ascii_table:
        for c2 in ascii_table:
            cs=[c1,c2]
            A, c = fsb_linearization(FSB64.Hcols, r, w, n, candidates=cs)
            if A.rank() == r:
                full_rank = True
                break
        if full_rank:
            break
    hsh_val = int.from_bytes(b"[FSB64] pre-image attack. Implemented by tl2cents in 2024.03.13!", byteorder='big')
    C = vector(GF(2), ZZ(c).digits(base=2, padto=r))
    H = vector(GF(2), ZZ(hsh_val).digits(base=2, padto=r))
    # solve x for Ax + C = H
    x = A.solve_right(H - C).list()
    x_num = ZZ([cs[idx] for idx in x], base=word_size)
    x_bytes = int(x_num).to_bytes(ZZ(word_size**w).nbits()//8, byteorder='big')    
    print("[+] Pre-image found: ", x_bytes)
    print(f"[+] Hash result {FSB64.digest(x_bytes) = }")
    print("-" * 80)

def test_fsb_collision_attack():
    # case r = 2w
    from itertools import combinations
    print("[+] Collision attack on FSB")
    r = 256
    w = 128
    word_size = 256
    n = word_size * w
    FSB = Fast_Syndrome_Based_Hash(r, w, n)
    letter_table = [i for i in range(ord("a"), ord("z") + 1)]
    table_iter1 = combinations(letter_table, 2)
    table_iter2 = combinations(letter_table, 2)
    for candidates1 in table_iter1:
        for candidates2 in table_iter2:
            if candidates1 == candidates2:
                continue
            A1, c1 = fsb_linearization(FSB.Hcols, r, w, n, candidates=candidates1)
            A2, c2 = fsb_linearization(FSB.Hcols, r, w, n, candidates=candidates2)
            A = block_matrix([A1, A2], ncols=2)
            c = c1^c2
            try:
                sol = A.solve_right(vector(GF(2), ZZ(c).digits(base=2, padto=r))).list()
                b, b_ = sol[:w], sol[w:]
                s1 = int(ZZ([candidates1[idx] for idx in b], base=word_size))
                s2 = int(ZZ([candidates2[idx] for idx in b_], base=word_size))
                if s1 == s2:
                    continue
                s1_bytes = s1.to_bytes(ZZ(word_size**w).nbits()//8, byteorder='big')
                s2_bytes = s2.to_bytes(ZZ(word_size**w).nbits()//8, byteorder='big')
                print(f"[+] collision found")
                print(f"[+] {s1_bytes = }\n[+] {s2_bytes = }")
                print(f"[+] Hash result {FSB.digest(s1_bytes) = }")
                print(f"[+] Hash result {FSB.digest(s2_bytes) = }")
                return
            # print error info
            except Exception as e:
                continue

if __name__ == "__main__":
    # test_FSB()
    test_fsb_pre_image()
    test_fsb_collision_attack()