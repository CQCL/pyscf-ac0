# CAS-AC0 module for PySCF
This is a PySCF extension implementing CAS-AC0 (see https://doi.org/10.1021/acs.jctc.8b00213), using refactored portions of GAMMCOR https://github.com/pernalk/GAMMCOR .

An executable script is provided, to enable calling cas-ac0 from inquanto-pyscf extension.

Note: install with `pip install -e .`. If installed as non-editable, the f2py compilation will fail. 

There are 2 tests included, executed with `python pyscf/cas_ac0/accas.py`.
