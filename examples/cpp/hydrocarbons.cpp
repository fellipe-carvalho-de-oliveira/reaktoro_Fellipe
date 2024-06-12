// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Fellipe Carvalho de Oliveira (05 May 2023)
//
// and since revised by:
//   •
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

#include <iomanip>

std::string const homedir = getenv("HOME");

MatrixXr kij{{0.,0.,0.}, 
             {0.,0.,0.},
             {0.,0.,0.}}; 

int main()
{   
    std::cout << std::fixed << std::setprecision(20);

    auto db = SupcrtDatabase::fromFile(homedir+
     "/miniconda3/reaktoro/embedded/databases/reaktoro/custom-supcrtbl.yaml"); 

    //LiquidPhase liquid("CO2(l) C7H16(l) C12H26(l)");
    //LiquidPhase liquid("C7H16(l) C12H26(l)");
    //liquid.setActivityModel(ActivityModelPengRobinson());

    LiquidPhase liquid("CO2(l) C7H16(l) C12H26(l)");
    liquid.setActivityModel(ActivityModelPcsaft(homedir+
     "/miniconda3/reaktoro/embedded/databases/pcsaft/params-pcsaft.yml",kij)); 

    GaseousPhase gas("CO2(g) C7H16(g) C12H26(g)");
    //GaseousPhase gas("C7H16(g) C12H26(g)");
    gas.setActivityModel(ActivityModelPcsaft(homedir+
    "/miniconda3/reaktoro/embedded/databases/pcsaft/params-pcsaft.yml",kij));
    //gas.setActivityModel(ActivityModelPengRobinson());    

    ChemicalSystem system(db, liquid, gas);
    //ChemicalSystem system(db, liquid);

    ChemicalState state(system);
    state.temperature(293.15, "kelvin");
    state.pressure(70, "MPa");
    state.set("CO2(l)", 0.5, "mol");
    state.set("C7H16(l)", 0.25, "mol");
    state.set("C12H26(l)", 0.25, "mol");
    state.set("CO2(g)", 0.5, "mol");
    state.set("C7H16(g)", 0.25, "mol");
    state.set("C12H26(g)", 0.25, "mol");


    EquilibriumSolver solver(system);
    solver.solve(state);

    // Output the chemical state to a text file
    state.output("state.txt");    

    ChemicalProps props(state);
    props.output("props.txt");

    std::cout << "Success! Check outputted files `props.txt`." << std::endl;
    std::cout << "The molar volume of Liquid Phase is: " << props.molarVolume() 
      << " m³/mol" << std::endl;
    std::cout << "The molar density of Liquid Phase is: " << 1./props.molarVolume() 
      << " mol/m³" << std::endl;       
    std::cout << "The density of Liquid Phase is: " << props.density() 
      << " kg/m³" << std::endl;      
 
    return 0;
}