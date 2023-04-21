# CAS-AC0 module for PySCF
This is a PySCF extension implementing CAS-AC0, using refactored portions of GAMMCOR https://github.com/pernalk/GAMMCOR .
It has to remain internal or be released under GPL, because GAMMCOR itself is released under GPL.
An executable script is (well, will be) provided, to enable calling cas-ac0 from inquanto-pyscf extension.
Note: install with `pip install -e .`. If installed as non-editable, the f2py compilation will fail.
