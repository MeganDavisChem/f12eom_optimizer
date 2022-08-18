# hoo rad ROHF-EOM-CCSD/aug-cc-pVTZ App
memory 7950 mb

molecule HNO{
  0 1
  H      0.000000000    1.442306288   -2.754574215
  C      0.000000000   -0.186314094   -1.423042607
  F      0.000000000    0.041269805    1.045804233

 units bohr
 }

set globals {
 g_convergence gau_verytight
 max_force_g_convergence 1.0e-8
 max_energy_g_convergence 1.0e-12
 max_disp_g_convergence 1.0e-8
 reference rohf
 basis aug-cc-pVDZ
 freeze_core true
 cachelevel=0
# print=2
 maxiter=500
 dertype none
 ints_tolerance 20
 roots_per_irrep = [0, 1]
}

set scf d_convergence 13
set ccenergy r_convergence 15
set cceom r_convergence 8

optimize('scf', dertype='energy')

