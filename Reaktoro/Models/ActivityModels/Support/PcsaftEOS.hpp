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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/StateOfMatter.hpp>

namespace Reaktoro {

// Tolerance for numeric calculations
const auto tol = 1.0E-9;

// Universal Constants
const auto R = universalGasConstant;
const auto n_av = avogadroNumber;
const auto kb = boltzmannConstant; 
const auto pi = piNumber; 

// Number of different association sites (acceptor or donor)
const unsigned s = 2;

// Constants for dispersion term
const static real a0[7] = { 0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084 };
const static real a1[7] = { -0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930};
const static real a2[7] = { -0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368};
const static real b0[7] = { 0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612};
const static real b1[7] = { -0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346};
const static real b2[7] = { 0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585};

/// Pcsaft types
enum PcsaftEOSModel {PCSAFT,A_PCSAFT};  
const static int n_pcsaftTypes = A_PCSAFT + 1; 

// Macros to auxiliary vectors
#define S(i,j) S[(j)+n_comp*(i)] 
#define delta(i,j,k,l) delta[(l)+n_comp*((k)+s*((j)+n_comp*(i)))]
#define deltaD(i,j,k,l) deltaD[(l)+n_comp*((k)+s*((j)+n_comp*(i)))]
#define deltaDD(i,j,k,l) deltaDD[(l)+n_comp*((k)+s*((j)+n_comp*(i)))]
#define deltaDDD(i,j,k,l) deltaDDD[(l)+n_comp*((k)+s*((j)+n_comp*(i)))]
#define deltaT(i,j,k,l) deltaT[(l)+n_comp*((k)+s*((j)+n_comp*(i)))]
#define deltaTT(i,j,k,l) deltaTT[(l)+n_comp*((k)+s*((j)+n_comp*(i)))]
#define delta_xm(i,j,k,l,m) delta_xm[(m)+n_comp*((l)+n_comp*((k)+s*((j)+n_comp*(i))))]
#define delta_xm_rho(i,j,k,l,m) delta_xm_rho[(m)+n_comp*((l)+n_comp*((k)+s*((j)+n_comp*(i))))]
#define delta_xm_rho_rho(i,j,k,l,m) delta_xm_rho_rho[(m)+n_comp*((l)+n_comp*((k)+s*((j)+n_comp*(i))))]
#define delta_xm_rho_rho_rho(i,j,k,l,m) delta_xm_rho_rho_rho[(m)+n_comp*((l)+n_comp*((k)+s*((j)+n_comp*(i))))]
#define delta_xm_T(i,j,k,l,m) delta_xm_T[(m)+n_comp*((l)+n_comp*((k)+s*((j)+n_comp*(i))))]
#define delta_xm_T_T(i,j,k,l,m) delta_xm_T_T[(m)+n_comp*((l)+n_comp*((k)+s*((j)+n_comp*(i))))]
#define XA(i,j) XA[(j)+n_comp*(i)] 
#define XA_old(i,j) XA_old[(j)+n_comp*(i)] 
#define lambda(i,j) lambda[(j)+order*(i)]
#define inv_lambda(i,j) inv_lambda[(j)+order*(i)]
#define psi2D(i,j) psi2D[(j)+n_comp*(i)]
#define XAT(i,j) XAT[(j)+n_comp*(i)]
#define XATT(i,j) XATT[(j)+n_comp*(i)]
#define XAD(i,j) XAD[(j)+n_comp*(i)]
#define XADD(i,j) XADD[(j)+n_comp*(i)]
#define XADDD(i,j) XADDD[(j)+n_comp*(i)]
#define XA_xm(i,j,k) XA_xm[(k)+n_comp*((j)+n_comp*(i))]
#define XA_xm_T(i,j,k) XA_xm_T[(k)+n_comp*((j)+n_comp*(i))]
#define XA_xm_T_T(i,j,k) XA_xm_T_T[(k)+n_comp*((j)+n_comp*(i))]
#define XA_xm_rho(i,j,k) XA_xm_rho[(k)+n_comp*((j)+n_comp*(i))]
#define XA_xm_rho_rho(i,j,k) XA_xm_rho_rho[(k)+n_comp*((j)+n_comp*(i))]
#define XA_xm_rho_rho_rho(i,j,k) XA_xm_rho_rho_rho[(k)+n_comp*((j)+n_comp*(i))]
#define Fprime(i,j) Fprime[(j)+order*(i)]


struct result_pcsaft_calculation {
    real P; // Pressure (Pa)
    ArrayXr ln_fugacity_coefficients;  // ln of fugacity coefficients of components in one phase (Dimensionless)
    real PD; // Derivative of pressure wrt density (Pa/(mol/m³))
    real PDD; // Second derivative of pressure wrt density (Pa / (mol/m³)²)
    real PDDD; // Second derivative of pressure wrt density (Pa / (mol/m³)³)
    real PT; // Derivative of pressure wrt temperature (Pa/K)
    real PTT; // Second Derivative of pressure wrt temperature (Pa/K/K)
    real Ares; // Residual Helmholtz Free energy (Dimensionless)
    real AresT; // Derivative of Residual Helmholtz Free energy wrt Temperature (1/K)
    real AresTT; // second derivative of Residual Helmholtz Free energy wrt Temperature (1/K/K)
    real Hres; // Residual Enthalpy (J/mol)
    real Gres; // Residual Gibbs Free energy (J/mol)
    real Z;  // Compressibility Factor (Dimensionless)
    real ZD; // Derivative of Compressibility Factor wrt density (m³/mol)
    real ZDD; // Second derivative of Compressibility Factor wrt density (m⁶/mol²)
    real ZDDD; // Third derivative of Compressibility Factor wrt density (m⁹/mol³)
    real ZT; // Derivative of Compressibility Factor wrt temperature (1/K)
    real ZTT; // Second Derivative of Compressibility Factor wrt temperature (1/K/K)
    real Cvres; // Isochoric heat capacity (J/mol/K)
    real kappaT; // Isothermal compressibility  (m²/N)
    real alpha; // Isobaric expansivity (1/K)
    real Cpres; // Isobaric heat capacity (J/mol/K)
};


/// The properties calculated by evaluating Pcsaft EOS.
struct PcsaftEOSProps
{
    /// The molar volume of the phase (in m3/mol).
    real V = {};

    /// The temperature derivative of the molar volume at constant pressure (in m3/(mol*K)).
    real VT = {};

    /// The pressure derivative of the molar volume constant temperature (in m3/(mol*Pa)).
    real VP = {};

    /// The residual molar Gibbs energy of the phase (in J/mol).
    real Gres = {};

    /// The residual molar enthalpy of the phase (in J/mol).
    real Hres = {};

    /// The residual molar heat capacity at constant pressure of the phase (in J/(mol*K)).
    real Cpres = {};

    /// The residual molar heat capacity at constant volume of the phase (in J/(mol*K)).
    real Cvres = {};

    /// The ln fugacity coefficients of the species in the phase.
    ArrayXr ln_phi;

    /// The state of matter of the fluid phase
    StateOfMatter som;
};

/// Calculates thermodynamic properties of fluid phases based on Pcsaft equation of state model.
class PcsaftEOS
{
public:
    /// The arguments needed to construct a PcsaftEOS object.
    struct Args
    {
        /// Matrix of BIP's
        ArrayXrConstRef kij;

        /// Map for Pcsaft types
        const Map<String,PcsaftEOSModel> types;
        
        /// Pcsaft types
        const Vec<String> pcsaft_type;

        /// The number of species in the fluid phase.
        const Index nspecies = {};

        /// The number of segments in a chain of the species
        ArrayXrConstRef m;

        /// The temperature-independent segment diameter of the species (in angstrom)
        ArrayXrConstRef sigma;

        /// The depth of the potential well or dispersion energy, epsilon/kb, (in K)
        ArrayXrConstRef epsilon;

        /// The effective association volume (dimensionless)
        ArrayXrConstRef kappa;

        /// The association well depth energy, energy/kb, (in K)
        ArrayXrConstRef eassoc;

        /// The number of electron donor sites of the species
        const Vec<unsigned> n_donor_sites;

        /// The number of electron acceptor sites of the species
        const Vec<unsigned> n_acceptor_sites;                               

    };

    /// Construct a PcsaftEOS instance.
    explicit PcsaftEOS(const Args& args);

    /// Construct a copy of a PcsaftEOS instance
    PcsaftEOS(const PcsaftEOS& other);

    /// Destroy this PcsaftEOS instance
    virtual ~PcsaftEOS();

    /// Assign a PcsaftEOS instance to this
    auto operator=(PcsaftEOS other) -> PcsaftEOS&;

    /// Compute the thermodynamic properties of the phase.
    /// @param[in] props The evaluated thermodynamic properties of the phase.
    /// @param T The temperature of the phase (in K)
    /// @param P The pressure of the phase (in Pa)
    /// @param x The mole fractions of the species in the phase (in mol/mol)
    auto compute(PcsaftEOSProps& props, real T, real P, ArrayXrConstRef x) -> void;

    /// Compute the thermodynamic properties of the phase based on T, rho, x.
    /// @param T The temperature of the phase (in K)
    /// @param rho The density of the phase (in mol/m³)
    /// @param x The mole fractions of the species in the phase (in mol/mol)
    auto pcsaft_calculation(result_pcsaft_calculation &pcsaftCalcProps, real T, real rho, ArrayXrConstRef x) -> void;

    /// Determine the intervals for calculating the roots of PD (derivative of pressure wrt density)
    auto zbrakPD(const real rho1, const real rho2, const int n,
      ArrayXr &rhob1, ArrayXr &rhob2, int &nroot,
      real T, ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void;

    /// determine the intervals for calculating the roots of P (derivative of pressure wrt density)
    auto zbrakP(const real rho1, const real rho2, const int n,
      ArrayXr &rhob1, ArrayXr &rhob2, int &nroot,
      real P, real T, ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void;
    
    /// determine the roots of P 
    auto newton_bissection(real &rho_c, real rho_a, real rho_b, real P, real T, 
      ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void;

    /// determine the roots of PD (derivative of pressure wrt density)  
    auto newton_bissectionD(real &rho_c, real rho_a, real rho_b, real T, 
      ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void;      

    // Calculate XA using Newton's Method for a system of non-linear equations
    auto calculate_XA(int n_comp, ArrayXrConstRef x,real num_den,Vec<int> &indx) -> void;

    auto resizing_vectors(result_pcsaft_calculation &pcsaftCalcProps, 
      int n_comp_New, int &n_comp_Old) -> void;

private:
    struct Impl;
    Ptr<Impl> pimpl;
    
    // Struct of thermodynamic Properties, private to the instance
    result_pcsaft_calculation pcsaftCalcProps;

    // Stores the number of components of the system in the PcsaftEOS instance between runs
    int n_comp0;

    // Auxiliary vectors
    ArrayXr d,dDT,dDTT;
    ArrayXr p_dj;
    ArrayXr ghs; // It is the radial pair distribution function for segments of components in the hard sphere system
    ArrayXr ghsD; // derivative of ghs wrt rho
    ArrayXr ghsDD; // second derivative of ghs wrt rho
    ArrayXr ghsDDD; // third derivative of ghs wrt rho
    ArrayXr ghsT; // derivative of ghs wrt T  
    ArrayXr ghsTT; // second derivative of ghs wrt T                                      
    ArrayXr Dghs; //rho*Dghs_Drho
    ArrayXr DghsD; // derivative of Dghs wrt rho
    ArrayXr DghsDD; // second derivative of Dghs wrt rho
    ArrayXr DghsDDD; // third derivative of Dghs wrt rho
    ArrayXr DghsT; // derivative of Dghs wrt T
    ArrayXr DghsTT; // second derivative of Dghs wrt T
    ArrayXr epsilon_ij;
    ArrayXr sigma_ij;
    ArrayXr mu_assoc;
    ArrayXr col;  
    ArrayXr psi;  
    ArrayXr dXA;  // Auxiliary vector for calculating all XA derivatives
    ArrayXr Dghsii_Dx;
    ArrayXr Dahs_Dx;
    ArrayXr Dadisp_Dx;
    ArrayXr Dahc_Dx;  
    ArrayXr mu_hc;
    ArrayXr mu_disp;  
    ArrayXr mu;
    ArrayXr F;        
    ArrayXr dif;  
    Vec<int> indx;    
    Vec<int> S;
    ArrayXr delta;
    ArrayXr deltaD;
    ArrayXr deltaDD;
    ArrayXr deltaDDD;
    ArrayXr deltaT;
    ArrayXr deltaTT;
    ArrayXr delta_xm;
    ArrayXr delta_xm_rho;
    ArrayXr delta_xm_rho_rho;
    ArrayXr delta_xm_rho_rho_rho;
    ArrayXr delta_xm_T;
    ArrayXr delta_xm_T_T;
    ArrayXr XA;
    ArrayXr XA_old; 
    ArrayXr lambda;
    ArrayXr inv_lambda;
    ArrayXr psi2D;
    ArrayXr XAT;
    ArrayXr XATT;
    ArrayXr XAD;
    ArrayXr XADD;
    ArrayXr XADDD;
    ArrayXr XA_xm;
    ArrayXr XA_xm_T;
    ArrayXr XA_xm_T_T;
    ArrayXr XA_xm_rho;
    ArrayXr XA_xm_rho_rho;
    ArrayXr XA_xm_rho_rho_rho;
    ArrayXr Fprime;
};

} // namespace Reaktoro
