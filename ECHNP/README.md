## Remarks

Implementations of paper  [Improving Bounds on Elliptic Curve Hidden Number Problem for ECDH Key Exchange](https://eprint.iacr.org/2022/1239.pdf).

- [echnp_solver.ipynb](./echnp_solver.ipynb): A step by step implementation for solving ECHNP with details.
- [lll.py](./lll.py), [logger.py](./logger.py), [rootfind_ZZ.py](./rootfind_ZZ.py) : a toolkit for finding the roots of reduced polynomials in $\mathbb{Z}$ from [kiona's coppersmith repository](https://github.com/kionactf/coppersmith).
- [findRootsZZ.py](./findRootsZZ.py) :  a toolkit for finding the roots of reduced polynomials in $\mathbb{Z}$ from [jvdsn-crypto-attacks](https://github.com/jvdsn/crypto-attacks/blob/master/shared/small_roots/__init__.py).
- [XHS_Lattice](./XHS_Lattice.py), [coppersmith_reduce.py](./coppersmith_reduce.py), [echnp_oracle.py](./echnp_oracle.py) : some basic functions for solving ECHNP.
- [echnp_solver](./echnp_solver.ipynb), [optimized_echnp_solver.py](./optimized_echnp_solver.py), [echnp_copper_tester.py](./echnp_copper_tester.py) : the main implementation of ECHNP solver and a detailed tester for coppersmith's method.


**If you want to learn about details, please directly refer to the [notebook](./echnp_solver.ipynb).**

**If you want to use it directly, please check the sage scripts under this directory.**

## CTF Challenges

There is a subdirectory [BabyDH2-Aliyun2024](./BabyDH2-Aliyun2024/) which contains a CTF challenge named `BabyDH2` which is completely based on the ECHNP problem.

- [task.py](./BabyDH2-Aliyun2024/task.py) : the task file for the challenge.
- [exp.py](./BabyDH2-Aliyun2024/exp.py) : a python (sage) script for solving the challenge locally.