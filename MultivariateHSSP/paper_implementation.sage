# -*- mode: python;-*-

# This is an implementation of the Nguyen-Stern algorithm, and of our new multivariate attack

# To run the Nguyen-Stern algorithm, run the function NSattack(). 
# To run all the experiments with the Nguyen-Stern algorithm, run the function statNS().
# We provide the experimental results below.

# To run our algorithm, run the function multiAttack().
# To run all the experiments with our algorithm, run the function statMulti().
# We provide the experimental results below.

from time import time

# https://github.com/Neobeo/HackTM2023/blob/main/solve420.sage
def flatter(M):
    from subprocess import check_output
    from re import findall
    M = matrix(ZZ,M)
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = '[[' + ']\n['.join(' '.join(map(str,row)) for row in M) + ']]'
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), map(int,findall(b'-?\\d+', ret)))

def genpseudoprime(eta,etamin=211):
    if eta<=(2*etamin):
        return random_prime(2**eta,False,2**(eta-1))
    else:
        return random_prime(2**etamin,False,2**(etamin-1))*genpseudoprime(eta-etamin)

def genParams(n=10,m=20,nx0=100):
    #print "Generation of x0",
    t=time()
    x0=genpseudoprime(nx0)
    #print time()-t

    # We generate the alpha_i's
    a=vector(ZZ,n)
    for i in range(n):
        a[i]=mod(ZZ.random_element(x0),x0)

    # The matrix X has m rows and must be of rank n
    while True:
        X=Matrix(ZZ,m,n)
        for i in range(m):
            for j in range(n):
                X[i,j]=ZZ.random_element(2)
        if X.rank()==n: break

    # We generate an instance of the HSSP: b=X*a
    c=vector(ZZ,[0 for i in range(m)])
    s=ZZ.random_element(x0)
    b=X*a
    for i in range(m):
        b[i]=mod(b[i],x0)

    return x0,a,X,b

# We generate the lattice of vectors orthogonal to b modulo x0
def orthoLattice(b,x0):
    m=b.length()
    M=Matrix(ZZ,m,m)

    for i in range(1,m):
        M[i,i]=1
    M[1:m,0]=-b[1:m]*inverse_mod(b[0],x0)
    M[0,0]=x0

    for i in range(1,m):
        M[i,0]=mod(M[i,0],x0)

    return M

def allones(v):
    if len([vj for vj in v if vj in [0,1]])==len(v):
        return v
    if len([vj for vj in v if vj in [0,-1]])==len(v):
        return -v
    return None

def recoverBinary(M5):
    lv=[allones(vi) for vi in M5 if allones(vi)]
    n=M5.nrows()
    for v in lv:
        for i in range(n):
            nv=allones(M5[i]-v)
            if nv and nv not in lv:
                lv.append(nv)
            nv=allones(M5[i]+v)
            if nv and nv not in lv:
                lv.append(nv)
    return Matrix(lv)

def allpmones(v):
    return len([vj for vj in v if vj in [-1,0,1]])==len(v)

# Computes the right kernel of M using LLL.
# We assume that m>=2*n. This is only to take K proportional to M.height()
# We follow the approach from https://hal.archives-ouvertes.fr/hal-01921335/document
def kernelLLL(M):
    n=M.nrows()
    m=M.ncols()
    # if m<2*n: return M.right_kernel().matrix()
    K=2**(m//2)*M.height()
  
    MB=Matrix(ZZ,m+n,m)
    MB[:n]=K*M
    MB[n:]=identity_matrix(m)
  
    # MB2=MB.T.LLL().T
    MB2=flatter(MB.T).T
    
    assert MB2[:n,:m-n]==0
    Ke=MB2[n:,:m-n].T

    return Ke
        
# This is the Nguyen-Stern attack, based on BKZ in the second step
def hssp_ns_attack(h, m, n, x0):
    b = h
    M=orthoLattice(b,x0)

    t=walltime()
    # M2=M.LLL()
    M2=flatter(M)
    print("LLL step1: %.1f" % walltime(t))

    MOrtho=M2[:m-n]

    print("  log(Height,2)=",int(log(MOrtho.height(),2)))

    t2=walltime()
    ke=kernelLLL(MOrtho)

    print("  Kernel: %.1f" % walltime(t2))
    print("  Total step1: %.1f" % walltime(t))

    if n>170: return

    beta=2
    tbk=walltime()
    while beta<n:
        if beta==2:
            M5=flatter(ke)
            # M5=ke.LLL()
        else:
            M5=M5.BKZ(block_size=beta)

        # we break when we only get vectors with {-1,0,1} components
        if len([True for v in M5 if allpmones(v)])==n: break

        if beta==2:
            beta=10
        else:
            beta+=10

    print("BKZ beta=%d: %.1f" % (beta,walltime(tbk)))
    t2=walltime()
    MB=recoverBinary(M5)
    print("  Recovery: %.1f" % walltime(t2))
    print("  Number of recovered vector=",MB.nrows())

    NS=MB.T
    invNSn=matrix(Integers(x0),NS[:n]).inverse()
    ra=invNSn*b[:n]
    print("  Total step2: %.1f" % walltime(tbk))
    print("  Total time: %.1f" % walltime(t))
    return ra
    
def test_hssp_ns_attack(n):
    print("NS Attack")
    m=int(max(2*n,16*log(n,2)))
    print("n=",n,"m=",m)
    iota=0.035
    nx0=int(2*iota*n**2+n*log(n,2))
    print("nx0=",nx0)

    x0,a,X,b=genParams(n,m,nx0)
    ra = hssp_ns_attack(b, m, n, x0)
    check = sorted(list(ra)) == sorted(list(a))
    print(f"{check = }")
    
# This is the Nguyen-Stern attack, based on BKZ in the second step
def NSattack(n=60):
    print("NS Attack")
    m=int(max(2*n,16*log(n,2)))
    print("n=",n,"m=",m)
    iota=0.035
    nx0=int(2*iota*n**2+n*log(n,2))
    print("nx0=",nx0)

    x0,a,X,b=genParams(n,m,nx0)
    ra = hssp_ns_attack(b, m, n, x0)
    check = sorted(list(ra)) == sorted(list(a))
    print(f"{check = }")

def statNS():
    for n in list(range(70,190,20)) + list(range(190,280,30)):
        NSattack(n)
        print()
    
def matNbits(M):
    return max([M[i,j].nbits() for i in range(M.nrows()) for j in range(M.ncols())])

# Matrix rounding to integers
def roundM(M):
    M2=Matrix(ZZ,M.nrows(),M.ncols())
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            M2[i,j]=round(M[i,j])
    return M2

def orthoLatticeMod(b,n,x0):
    m=b.length()
    assert m>=3*n
    assert m % n==0
    M=Matrix(ZZ,m,3*n)
    M[:2*n,:2*n]=identity_matrix(2*n)
    for i in range(2,m//n):
        M[i*n:(i+1)*n,2*n:3*n]=identity_matrix(n)

    M[1:,0]=-b[1:]*inverse_mod(b[0],x0)
    M[0,0]=x0

    for i in range(1,m):
        M[i,0]=mod(M[i,0],x0)
    return M

def NZeroVectors(M):
    return sum([vi==0 and 1 or 0 for vi in M])

# This is our new multivariate attack        
def hssp_multiAttack(h, m, n, x0, k=4):

    print("n=",n,"m=",m,"k=",k)
    b = h

    M=orthoLatticeMod(b,n,x0)

    print("Step 1")
    t=walltime()

    # M[:n//k,:n//k]=M[:n//k,:n//k].LLL()
    M[:n//k,:n//k] = flatter(M[:n//k,:n//k])
    

    # M2 = M[:2*n,:2*n].LLL()
    M2 = flatter(M[:2*n,:2*n])
    
    tprecomp=walltime(t)
    print("  LLL:%.1f" % tprecomp)

    RF=RealField(matNbits(M))

    M4i=Matrix(RF,M[:n//k,:n//k]).inverse()
    M2i=Matrix(RDF,M2).inverse()

    ts1=walltime()
    while True:
        flag=True
        for i in range((m/n-2)*k):
            indf=2*n+n//k*(i+1)
            if i==(m//n-2)*k-1:
                indf=m

            mv=roundM(M[2*n+n//k*i:indf,:n//k]*M4i)
            if mv==0: 
                continue
            flag=False
            M[2*n+n//k*i:indf,:]-=mv*M[:n//k,:]
        if flag: break
    print("  Sred1:%.1f" % walltime(ts1))

    M[:2*n,:2*n]=M2

    ts2=walltime()
    while True:
        #print("  matNBits(M)=",matNbits(M[2*n:]))
        mv=roundM(M[2*n:,:2*n]*M2i)
        if mv==0: break
        M[2*n:,:]-=mv*M[:2*n,:]
    print("  Sred2:%.1f" % walltime(ts2))

    # Orthogonal of the orthogonal vectors
    # We compute modulo 3
    MO=Matrix(GF(3),n,m)

    tk=walltime()
    MO[:,:2*n]=kernelLLL(M[:n,:2*n])
    print("  Kernel LLL: %.1f" % walltime(tk))

    for i in range(2,m/n):
        MO[:,i*n:(i+1)*n]=-(M[i*n:(i+1)*n,:2*n]*MO[:,:2*n].T).T
    #print("Total kernel computation",walltime(tk))
    print("  Total Step 1: %.1f" % walltime(t))

    print("Step 2")
    t2=walltime()
    xt23=Matrix(GF(3),[(-x).list()+[x[i]*x[j]*((i==j) and 1 or 2) for i in range(n) for j in range(i,n)] for x in MO.T])
    ke3=xt23.right_kernel().matrix()
    print("  Kernel: %.1f" % walltime(t2))

    assert xt23.nrows()==m
    assert xt23.ncols()==n*(n+1)//2+n

    ke23=Matrix(GF(3),n,n*n)
    ind=n
    for i in range(n):
        for j in range(i,n):
            ke23[:,i*n+j]=ke3[:,ind]
            ke23[:,j*n+i]=ke3[:,ind]
            ind+=1

    tei=walltime()
    # We will compute the list of eigenvectors
    # We start with the full space.
    # We loop over the coordinates. This will split the eigenspaces.
    li=[Matrix(GF(3),identity_matrix(n))]    
    for j in range(n):       # We loop over the coordinates of the wi vectors.
        #print("j=",j)
        M=ke23[:,j*n:(j+1)*n]   # We select the submatrix corresponding to coordinate j
        li2=[]                 # We initialize the next list
        for v in li:
            if v.nrows()==1:     # We are done with this eigenvector 
                li2.append(v)
            else:     # eigenspace of dimension >1
                #print("eigenspace of dim:",v.nrows())
                A=v.solve_left(v*M)  # v*M=A*v. When we apply M on the right, this is equivalent to applying the matrix A.
                                      # The eigenvalues of matrix A correspond to the jth coordinates of the wi vectors in that
                                      # eigenspace
                for e,v2 in A.eigenspaces_left():    # We split the eigenspace according to the eigenvalues of A.
                    vv2=v2.matrix()
                    #print("  eigenspace of dim:",(vv2*v).nrows())
                    li2.append(vv2*v)                   # The new eigenspaces 

        li=li2

    #print("Eigenvectors computation",walltime(tei))

    NS=Matrix([v[0] for v in li])*MO
    for i in range(n):
        if any(c==2 for c in NS[i]): NS[i]=-NS[i]

    print("  Number of recovered vectors:",NS.nrows())

    NS=NS.T

    # b=X*a=NS*ra
    invNSn=matrix(Integers(x0),NS[:n]).inverse()
    ra=invNSn*b[:n]
    print("  Total step2: %.1f" % walltime(t2))
    print("  Total runtime: %.1f" % walltime(t))
    return ra
    
def test_multi_attack(n):
    if n % 2==1:
        m=n*(n+3)//2 # n is odd
    else:
        m=n*(n+4)//2 # n is even
    k=4
    iota=0.035
    nx0=int(2*iota*n**2+n*log(n,2))
    x0,a,X,b=genParams(n,m,nx0)
    print("multiAttack")
    ra = hssp_multiAttack(b, m, n, x0)
    check = sorted(list(ra)) == sorted(list(a))
    print(f"{check = }")

def multiAttack(n=16):
    if n % 2==1:
        m=n*(n+3)//2 # n is odd
    else:
        m=n*(n+4)//2 # n is even
    k=4
    print("multiAttack")

    print("n=",n,"m=",m,"k=",k)

    iota=0.035
    nx0=int(2*iota*n**2+n*log(n,2))
    print("nx0=",nx0)

    x0,a,X,b=genParams(n,m,nx0)
    ra = hssp_multiAttack(b, m, n, x0)
    check = sorted(list(ra)) == sorted(list(a))
    print(f"{check = }")

def statMulti():
    for n in list(range(70,210,20)) + [220,250]:
        multiAttack(n)
        print()

def _statMulti():
    for n in list(range(70,120,20)):
        multiAttack(n)
        print()
        
def _statNS():
    for n in list(range(70,120,20)):
        NSattack(n)
        print()
        

if __name__ == "__main__":
    # statNS()
    # statMulti()
    test_hssp_ns_attack(64)
    print()
    test_multi_attack(64)
    print()
    test_multi_attack(190)
    print()