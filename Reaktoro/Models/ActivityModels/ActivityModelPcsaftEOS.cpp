// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "ActivityModelPcsaftEOS.hpp"

// C++ includes
#include <fstream>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Models/ActivityModels/Support/PcsaftEOS.hpp>

namespace Reaktoro {

using std::log;

auto activityModelPcsaftEOS(const SpeciesList& species, String file, MatrixXr kij) -> ActivityModel
{
    // The number of species
    const auto nspecies = species.size();

    // Enum and Map for pcsaft types
    Map<String,PcsaftEOSModel> types{{"PCSAFT", PCSAFT},{"A_PCSAFT", A_PCSAFT}};

    // Get the parameters of the Pcsaft EOS
    Vec<String> pcsaft_type(nspecies);
    ArrayXr m(nspecies), sigma(nspecies), epsilon(nspecies), kappa(nspecies),
      eassoc(nspecies);
    Vec<unsigned> n_donor_sites(nspecies), n_acceptor_sites(nspecies);  

    read_pcsaft_database(species, file, types, pcsaft_type, m, sigma, 
      epsilon, kappa, eassoc, n_donor_sites, n_acceptor_sites);

    /*for (int j=0;j<nspecies;j++) std::cout << types[pcsaft_type[j]] << " " << m(j) << " " 
    << sigma(j) << " " << epsilon(j) << " " << kappa(j) << " "
    << eassoc(j) << " " << n_donor_sites[j] << " "
    << n_acceptor_sites[j] << " " << std::endl;*/

    // Transform kij(matrix) into k_ij(vector)
    auto count = 0;
    ArrayXr k_ij(nspecies*nspecies);
    for(auto i = 0; i < nspecies; i++){
      for(auto j = 0; j < nspecies; j++){
        k_ij[count] =  kij(i,j);
        count++;
      }
    }        

    // Initialize the PcsaftEOS instance
    PcsaftEOS eos({k_ij,types,pcsaft_type,nspecies, m, sigma, epsilon, kappa, eassoc, n_donor_sites, n_acceptor_sites});

     /// The thermodynamic properties calculated with PcsaftEOS
    PcsaftEOSProps res;
    res.ln_phi.resize(nspecies);

    // Define the activity model function of the phase
    ActivityModel model = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        const auto Pbar = P * 1.0e-5; // convert from Pa to bar

        eos.compute(res, T, P, x);

        props.Vx   = res.V;
        props.VxT  = res.VT;
        props.VxP  = res.VP;
        props.Gx   = res.Gres;
        props.Hx   = res.Hres;
        props.Cpx  = res.Cpres;
        props.ln_g = res.ln_phi;
        props.ln_a = res.ln_phi + log(x) + log(Pbar);
        props.som  = res.som;
    }; 

    return model;
}

auto ActivityModelPcsaftEOS(String file, MatrixXr kij) -> ActivityModelGenerator
{
    return [=](const SpeciesList& species)
    {
        return activityModelPcsaftEOS(species, file, kij);
    };
}

auto ActivityModelPcsaft(String file, MatrixXr kij) -> ActivityModelGenerator
{ 
    return ActivityModelPcsaftEOS(file, kij);
}


/// Return a substance name corrected for checking in the ChemicalProps database.
/// This function removes suffix from given substance name. For
/// example, `H2O(aq)` becomes `H2O`. It also replaces spaces by dashes.
auto correctName(String name) -> String
{
    const auto [name0, _] = splitSpeciesNameSuffix(name);
    name = name0;
    name = replace(name, " ", "-");
    return name;
}

auto split_string (const String &s) -> Vec<String>
{
  Vec<String> result;
  Vec<char> line (s.begin(),s.end());
  int n = -1;
  for (int i = 0; i < line.size(); i++)
  {
    if (std::isgraph(line[i])) //(line[i] != ' ') 
    {
      if ((i!=0) && (line[i-1]==' '))
      {
        result.push_back(String{line[i]});
        n++;
      }
      else
      {
        result[n].append(String{line[i]});
      } 
    } 
  }
  return result;
}

auto get_param(std::ifstream &file) -> double
{
  String line;
  std::getline (file,line);
  auto parameters = split_string(line);
  return std::stod(parameters[1]);
}

auto read_pcsaft_database(const SpeciesList& species, String file,
      Map<String,PcsaftEOSModel> types,Vec<String> &pcsaft_type, ArrayXr &m, ArrayXr &sigma, 
      ArrayXr &epsilon, ArrayXr &kappa, ArrayXr &eassoc, Vec<unsigned> &n_donor_sites, 
      Vec<unsigned> &n_acceptor_sites) -> void
{
  // The number of species
  const auto nspecies = species.size();

//for (int j = 0; j < nspecies; j++) 
//        std::cout << species[j].name() <<  " " << correctName(species[j].name()) << std::endl;

  // The number of parameters of pcsaft
  int n_param;  
  
  std::ifstream input_file (file);
  if (!input_file.is_open()) {
    std::cout << "Error opening Pcsaft database file: `" << file << "`." << std::endl; 
    exit(0);
  }

  auto i = 0;
  bool test;
  String line;
  for (int i=0;i<2;i++) input_file.ignore(100,'\n'); // skipping 2 first lines of Pcsaft database file    
  while(getline (input_file,line))
  { 
    // getting species name
    //std::cout << line << std::endl;
    auto substring = split_string(line);
    //for (int j = 0; j < substring.size(); j++) std::cout << substring[j] << std::endl;
    
    // getting species pcsaft type
    getline (input_file,line);
    pcsaft_type[i] = split_string(line)[1];
    //std::cout << pcsaft_type[i] << std::endl;      

    switch (types[pcsaft_type[i]])
    {
    case (PCSAFT):
      n_param = 3;
      test = false;
      for (int j = 0; j < nspecies; j++)
      {
        if (substring[2]==correctName(species[j].name()))
        { 
          m(j) = get_param(input_file); 
          sigma(j) = get_param(input_file);
          epsilon(j) = get_param(input_file);
          //std::cout << m(j) << " " << sigma(j) << " " << epsilon(j) << " " << std::endl;
          i++; 
          test = true;
          break;
        } 
      }
      break;
    case (A_PCSAFT):
      n_param = 7;
      test = false;
      for (int j = 0; j < nspecies; j++)
      {
        if (substring[2]==correctName(species[j].name()))
        { 
          m(j) = get_param(input_file); 
          sigma(j) = get_param(input_file);
          epsilon(j) = get_param(input_file);
          kappa(j) = get_param(input_file);
          eassoc(j) = get_param(input_file);
          n_donor_sites[j] = get_param(input_file);
          n_acceptor_sites[j] = get_param(input_file);
          //std::cout << m(j) << " " << sigma(j) << " " << epsilon(j) << " " 
          //  << eassoc(j) << " " << n_acceptor_sites[j] << " " << std::endl;
          i++; 
          test = true;
          break;
        } 
      }      
      break;       
    default:
      //not working
      std::cout << "Not available pcsaft model named " << pcsaft_type[i] << std::endl;
      break;
    }
    if (i==nspecies) break;
    if (!test) {for (int k=0;k<n_param;k++) input_file.ignore(100,'\n');}
  }
  input_file.close();
}


} // namespace Reaktoro
