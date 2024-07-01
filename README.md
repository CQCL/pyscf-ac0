# CAS-AC0 module for PySCF
This is a [PySCF](https://pyscf.org/) extension implementing [CAS-AC0](https://doi.org/10.1021/acs.jctc.8b00213), using refactored 
portions of [GAMMCOR](https://github.com/pernalk/GAMMCOR). GAMMCOR and the ACn theory have been developed by the Pernal 
group at Łódź University of Technology

# Getting started

`pyscf-ac0` is available on pypi for Python 3.10, 3.11, and 3.12 on Linux and MacOS. Install with:
```
pip install pyscf-ac0
```
A command line executable `rdm_ac0` is provided for interfacing with the
[`inquanto-pyscf`](https://inquanto.quantinuum.com/) extension to compute the AC0 correlation energy. Use this 
executable with: 
```
rdm_ac0 file.h5
```
where `file.h5` is a data file produced by `inquanto-pyscf`.

# Installing from source

`pyscf-ac0` uses the [scikit-build-core](https://github.com/scikit-build/scikit-build-core) build system. This requires 
a C and Fortran compiler. Try installing from source with:
``` 
pip install .
```
If this runs into issues, you may need to be more deliberate with compilers. Specify C and Fortran compilers with:
```
pip install -v . -Ccmake.args="-DCMAKE_C_COMPILER=/my/C/compiler; -DCMAKE_Fortran_COMPILER=/my/fortran/compiler"
```
On linux, we test with gcc-10 for C and Fortran. And on MacOS, we test with AppleClang 15 for C, and gcc-11 for Fortran.

# Development

Install packages for development with
```shell
pip install -r tests/test-requirements.txt
```
This repository comes with two simple tests. Once installed from source, run tests from the project root with
``` 
pytest
```


