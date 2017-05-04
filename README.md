# Atom
The program is for calculation of atomic energy using Gaussian basis set.
The program reads a basis set file (obtained from https://bse.pnl.gov/bse/portal in Gaussian94 format) for maximum 7 atoms.
The basis set file should be specified as the first program argument. 
I just started this program, hence so far it can calculate only hydrogen-like atoms :) For them, the energy is calculated analytically,
no numerical integration is used so far. Currently, I am adding coulomb and exchange integrals. 
