"""kinetics of microbially induced calcium carbonate precipitation(MICP)."""

import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

phreeqc = phreeqc_mod.IPhreeqc('/usr/local/lib/libiphreeqc.so')

initial_conditions = """
#-------------------------------------------
SOLUTION_MASTER_SPECIES
    Urea Urea  0 CO(NH2)2 60.06 
    
SOLUTION_SPECIES 
Urea = Urea
        log_k  0 
        Vm  14.5
Cl- + H+ = HCl
        log_k   0.7
        delta_h 0       kcal
        -gamma  0       0.06
Cl- + Ca+2 = CaCl+
        log_k   0.3
        delta_h 0       kcal
        -gamma  0       0.06
2Cl- + Ca+2 = CaCl2
        log_k   .66
        delta_h 0       kcal
        -gamma  0       0.06
#---------------------------------------------

SOLUTION 1
      Urea 1
      Ca 1
      Cl 2
      -units 	mol/kgw
      -pH 	7	charge_balance

Rates
    Ureolysis
    -START
    10 REM PARM(1) = biomass density in CFU/l
    20 if (M <= 0) then goto 40                                              
    30 rate = PARM(1)*6.4e-12 /3600 * mol("Urea") / (mol("Urea") + .3)
    40 moles = rate * TIME                                                   
    50 save moles                                                            
    -END

    Calcite
    -START
    10  REM PARM(1) = specific surface area of calcite, cm^2/mol
    20  REM PARM(2) = exponent for M/M0
    30  si_cc = SI("Calcite")
    40  IF (M <= 0  and si_cc < 0) THEN GOTO 100
    50  k1 = 1.55e-6 #mol/m2/s
    60  k2 = 5.01E-01 #mol/m2/s
    70  IF M0 > 0 THEN area = PARM(1)*M0*(M/M0)^PARM(2) ELSE area = PARM(1)*M #area(m²)
    80  rate = area * (k1 + ACT("H+")* k2) * (1 - 10^(si_cc)) #mol/s
    90  moles = rate * TIME 
    100 SAVE moles
    -END

Kinetics 1 #Ureolysis 
    Ureolysis
        -formula AmmH+ 2 CO3-2 1 Urea -1 H2O -2
        -parms 1.58e9 
    Calcite
        -m0   .04 # mol/L # same as lb
        -parms 1.1e1   0.6
       -steps 900000 s in 100 steps 
 
USER_PUNCH
    -heading  soln_vol 
    -START
    10 punch soln_vol
    -END

SELECTED_OUTPUT
    -step    true
    -distance fasle
    -state  false
    -sim    false
    -soln   false
    -pe     false
    -time   true
    -pH 
    -tot    Ca Urea 
    -kinetic_reactants Calcite Ureolysis
    -si calcite

END
"""
np.set_printoptions(precision=2)
pd.set_option('display.precision', 2)

phreeqc.load_database(r"/home/amin/phreeqc_databases/Amm_biozement.dat")
phreeqc.run_string(initial_conditions)
output = phreeqc.get_selected_output_array()
table = pd.DataFrame(output[1:], columns=output[0])
print("\n")
print(table)

dk_Calcite = table['dk_Calcite'][1:]
k_Calcite = table['k_Calcite'][1:]
dk_Ureolysis = table['dk_Ureolysis'][1:]
time = table['time'][1:]
step = table['step'][1:]
R_p = dk_Calcite/time/step
R_u = -dk_Ureolysis/time/step
sv_calcite = 36.9 # specific volume of calcite (cm³/mol)
tot_v = table['soln_vol'][1:] + k_Calcite*sv_calcite/1000 # total volume in litre

fig, ax = plt.subplots()
ax.plot(time / 3600, R_p / R_u)
ax.set(xlabel='time(h)', ylabel='$R_p/R_u$')

plt.show()
