# CAS-AC0 module for PySCF
This is a PySCF extension implementing CAS-AC0 (see https://doi.org/10.1021/acs.jctc.8b00213), using refactored 
portions of GAMMCOR https://github.com/pernalk/GAMMCOR . GAMMCOR and the ACn theory have been developed by the Pernal 
group at Łódź University of Technology

An executable script is provided, to enable calling cas-ac0 from inquanto-pyscf extension. Tested with PySCF 2.4.0 and 2.5.0.

Install with `pip install .`

There are 2 tests included, executed with `python pyscf/cas_ac0/accas.py`.

Use the command-line executable with `rdm_ac0 file.h5`.
