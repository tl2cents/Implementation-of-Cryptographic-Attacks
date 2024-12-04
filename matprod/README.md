
## Remarks

Implementations of attacks on paper [New Public-Key Cryptosystem Blueprints Using Matrix Products in Fp](https://eprint.iacr.org/2023/1745.pdf). The scripts are my exploits from [HITCON CTF 2024 MatProd](https://github.com/maple3142/My-CTF-Challenges/tree/master/HITCON%20CTF%202024/MatProd) and the attack ideas originate from [@maple3142](https://github.com/maple3142/). **For more details, please refer to my blog post [here](https://blog.tanglee.top/2024/07/15/HITCON-CTF-2024-Qual-Crypto-Writeup.html#matprod).** There is also a paper proposing a smiliar attack [Attacking trapdoors from matrix products](https://eprint.iacr.org/2024/1332).

- [Attack.ipynb](./attack.ipynb): detailed output of the attack.
- [break_alternating.py](./break_alternating.py) : script for breaking the alternating matrix product cryptosystem.
- [break_direct.py](./break_direct.py) : script for breaking the direct matrix product cryptosystem.
- [chall.py](./chall.py) : the challenge script from HITCON CTF 2024 MatProd.
- [2024 hitcon matprod writeup](https://blog.tanglee.top/2024/07/15/HITCON-CTF-2024-Qual-Crypto-Writeup.html#matprod): my writeup.

**If you want to learn about details, please directly refer to the [notebook](./attack.ipynb).**

**If you want to use it directly, please check the python/sage scripts under this directory.**