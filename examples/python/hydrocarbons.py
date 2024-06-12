# Fellipe Carvalho de Oliveira - 23/05/2023

from reaktoro import *

#import reaktoro
#print(reaktoro.__path__)
#from reaktplot import *
#import numpy as np

import os

homedir = os.path.expanduser("~")

kij = [[0.,0.,0.], 
      [0.,0.,0.],
      [0.,0.,0.]]

db = SupcrtDatabase.fromFile(homedir+"/miniconda3/reaktoro/embedded/databases/reaktoro/custom-supcrtbl.yaml")

#liquid = LiquidPhase("CO2(l) C7H16(l) C12H26(l)")
#liquid = ("C7H16(l) C12H26(l)")
#liquid.setActivityModel(ActivityModelPengRobinson())

liquid = LiquidPhase("CO2(l) C7H16(l) C12H26(l)")
liquid.setActivityModel(ActivityModelPcsaft(homedir+
  "/miniconda3/reaktoro/embedded/databases/pcsaft/params-pcsaft.yml",kij))

gas = GaseousPhase("CO2(g) C7H16(g) C12H26(g)")
# GaseousPhase gas("C7H16(g) C12H26(g)")
gas.setActivityModel(ActivityModelPcsaft(homedir+
  "/miniconda3/reaktoro/embedded/databases/pcsaft/params-pcsaft.yml",kij))
# gas.setActivityModel(ActivityModelPengRobinson())    

system = ChemicalSystem(db,liquid,gas)
# ChemicalSystem system(db, liquid)

state = ChemicalState(system)
state.temperature(293.15, "kelvin")
state.pressure(70, "MPa")
state.set("CO2(l)", 0.5, "mol")
state.set("C7H16(l)", 0.25, "mol")
state.set("C12H26(l)", 0.25, "mol")
state.set("CO2(g)", 0.5, "mol")
state.set("C7H16(g)", 0.25, "mol")
state.set("C12H26(g)", 0.25, "mol")


solver = EquilibriumSolver(system)
solver.solve(state)

# Output the chemical state to a text file
state.output("state.txt")    

props = ChemicalProps(state)
props.output("props.txt")

print("Success! Check outputted files `props.txt`.")
print("The molar volume of Liquid Phase is: " + props.molarVolume() + " m³/mol")
print("The molar density of Liquid Phase is: " + 1./props.molarVolume() + " mol/m³")
print("The density of Liquid Phase is: " + props.density() + " kg/m³")


