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

#include "PcsaftEOS.hpp"

// C++ includes
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <thread>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <string>
#include <chrono>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/StateOfMatter.hpp>

//=================================================================================================
// == REFERENCES ==
//=================================================================================================
// The implementation of this module is based on the following articles and codes:
//
//     Gross, J., & Sadowski, G. (2001). Perturbed-chain SAFT: An equation of state based 
//     on a perturbation theory for chain molecules. Industrial & engineering chemistry 
//     research, 40(4), 1244-1260.
//     
//     Privat, R., Gani, R., & Jaubert, J. N. (2010). Are safe results obtained when the 
//     PC-SAFT equation of state is applied to ordinary pure chemicals?. Fluid Phase Equilibria, 295(1), 76-92.
//
//     Gross, J., & Sadowski, G. (2002). Application of the perturbed-chain SAFT equation of state to 
//     associating systems. Industrial & engineering chemistry research, 41(22), 5510-5515.
//
//     Huang, S. H., & Radosz, M. (1991). Equation of state for small, large, polydisperse, and associating
//     molecules: extension to fluid mixtures. Industrial & Engineering Chemistry Research, 30(8), 1994-2005.
//
//     Huang, S. H., & Radosz, M. (1990). Equation of state for small, large, polydisperse, and associating 
//     molecules. Industrial & Engineering Chemistry Research, 29(11), 2284-2294.
//
//     Chapman, W. G., Gubbins, K. E., Jackson, G., & Radosz, M. (1990). New reference equation of state 
//     for associating liquids. Industrial & engineering chemistry research, 29(8), 1709-1721.
//
//     Tan, S. P., Adidharma, H., & Radosz, M. (2004). Generalized procedure for estimating the fractions of 
//     nonbonded associating molecules and their derivatives in thermodynamic perturbation theory. Industrial 
//     & engineering chemistry research, 43(1), 203-208.
//       
//     Michelsen, M. L., & Mollerup, J. (2004). Thermodynamic Modelling: Fundamentals and 
//     Computational Aspects. Tie-Line Publications.
//
//
//     https://pcsaft.readthedocs.io/en/latest/index.html
//
//     https://www.th.bci.tu-dortmund.de/cms/de/Forschung/PC-SAFT/Download/index.html
//
//     http://hpp.uva.es/open-source-software-eos/
//
//=================================================================================================

namespace Reaktoro {

using std::abs;
using std::log;
using std::sqrt;
using std::pow;

// Auxiliary functions
namespace detail 
{

  auto ludcmp(ArrayXr &a, int n, Vec<int> &indx, real &d) -> void
  {
    const real TINY=1.0e-20;
    int i,imax,j,k;
    real big,dum,sum,temp;

    static ArrayXr vv(n);
    d=1.0;
    for (i=0;i<n;i++) 
    {
      big=0.0;
      for (j=0;j<n;j++) if ((temp=abs(a[j+n*i])) > big) big=temp;
      if (big == 0.0) std::cout << "Singular matrix in routine ludcmp" << std::endl; 
      vv[i]=1.0/big;
    }
    for (j=0;j<n;j++) 
    {
      for (i=0;i<j;i++) 
      {
        sum=a[j+n*i];
        for (k=0;k<i;k++) sum -= a[k+n*i]*a[j+n*k];
        a[j+n*i]=sum;
      }
      big=0.0;
      for (i=j;i<n;i++) 
      {
        sum=a[j+n*i];
        for (k=0;k<j;k++) sum -= a[k+n*i]*a[j+n*k];
        a[j+n*i]=sum;
        if ((dum=vv[i]*abs(sum)) >= big) 
        {
          big=dum;
          imax=i;
        }
      }
      if (j != imax) 
      {
        for (k=0;k<n;k++) 
        {
          dum=a[k+n*imax];
          a[k+n*imax]=a[k+n*j];
          a[k+n*j]=dum;
        }
        d = -d;
        vv[imax]=vv[j];
      }
      indx[j]=imax;
      if (a[j+n*j] == 0.0) a[j+n*j]=TINY;
      if (j != n-1) 
      {
        dum=1.0/(a[j+n*j]);
        for (i=j+1;i<n;i++) a[j+n*i] *= dum;
      }
    }
  
  }

  auto lubksb(ArrayXr &a, int n, Vec<int> &indx, ArrayXr &b) -> void
  {
    int i,ii=0,ip,j;
    real sum;

    for (i=0;i<n;i++) 
    {
      ip=indx[i];
      sum=b[ip];
	  b[ip]=b[i];
	  if (ii != 0)  for (j=ii-1;j<i;j++) sum -= a[j+n*i]*b[j];
	  else if (sum != 0.0) ii=i+1;
	  b[i]=sum;
    }
    for (i=n-1;i>=0;i--) 
    {
	  sum=b[i];
	  for (j=i+1;j<n;j++) sum -= a[j+n*i]*b[j];
	  b[i]=sum/a[i+n*i];
    }
  
  } 

}

struct PcsaftEOS::Impl
{

  /// The number of species in the phase.
  unsigned nspecies;

  ArrayXr kij;
  Map<String,PcsaftEOSModel> types;
  Vec<String> pcsaft_type;  
  ArrayXr m;
  ArrayXr sigma;
  ArrayXr epsilon;
  ArrayXr kappa;
  ArrayXr eassoc;
  Vec<unsigned> n_donor_sites;      
  Vec<unsigned> n_acceptor_sites;

  result_pcsaft_calculation pcsaftCalcProps; 
  int n_comp0; 

  //Auxiliary vectors;
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

  /// Construct a PcsaftEOS::Impl instance.
  Impl(const Args& args)
  : kij(args.kij),
    types(args.types),
    pcsaft_type(args.pcsaft_type),
    nspecies(args.nspecies),
    m(args.m),
    sigma(args.sigma),
    epsilon(args.epsilon),
    kappa(args.kappa),
    eassoc(args.eassoc),
    n_donor_sites(args.n_donor_sites),
    n_acceptor_sites(args.n_acceptor_sites)
  {
    assert(kij.size() == nspecies*nspecies);  
    assert(types.size() == n_pcsaftTypes);
    assert(pcsaft_type.size() == nspecies);  
    assert(m.size() == nspecies);
    assert(sigma.size() == nspecies);
    assert(epsilon.size() == nspecies);
    assert(kappa.size() == nspecies);
    assert(eassoc.size() == nspecies);
    assert(n_donor_sites.size() == nspecies);
    assert(n_acceptor_sites.size() == nspecies);
  }

  auto compute(PcsaftEOSProps& props, real T, real P, ArrayXrConstRef x) -> void
  { 
    
    auto n_comp = x.size();

    // Check if the mole fractions are zero or non-initialized  
    if(n_comp == 0 || x.maxCoeff() <= 0.0)  return;  
  
    // Resizing vectors if needed
    if (pcsaftCalcProps.ln_fugacity_coefficients.size() == 0) n_comp0 = 0;
    resizing_vectors(pcsaftCalcProps,n_comp,n_comp0);

    static ArrayXr rho_range(2); // Reduced densities or packing fractions, eta, gas and liquid limits

    // Calculating the factor to be multiplied to reduced densites to get densities in mol/m³
    real factor,di,summation = 0.;
    for (int i = 0; i < n_comp; i++) {
      di = sigma(i)*(1.-0.12*exp(-3.*epsilon(i)/T)); // Angs³ EQ. 3
      summation += x[i]*m(i)*(di*di*di); // Angs³/parti
    }
    factor = 6./pi/summation*1.0e30/n_av; // mol/m³, Based on EQ. 9
    rho_range[0] = 1.e-10*factor;   // low boundary for finding densities (in mol/m³), 1e-10 is the reduced low boundary density 
    rho_range[1] = 0.7405*factor; // upper boundary for finding densities (in mol/m³), 0,7405 is the reduced upper boundary density (max packing)
    
    // Zbrak method for determining the nº of roots of PD (Pressure derivative wrt density) 
    // within the reduced density range 1.e-10 - 0.7405
    auto n_intervals = 15;
    auto max_roots = 8;  // 6 roots + P(0) + P(inf)
    int n_rootsPD; // nº of roots of the Pressure derivative wrt density 
    
    // The PD (pressure derivative) density roots will be between rhob1[0] and rhob2[0], rhob1[1] and rhob2[1] ... and so on 
    static ArrayXr rhob1(max_roots);  // the nº of elements of this vector is the max number of roots to be found which should be the number of phases
    static ArrayXr rhob2(max_roots);
    std::fill(rhob1.begin(),rhob1.end(),0.);
    std::fill(rhob2.begin(),rhob2.end(),0.);
    zbrakPD(rho_range[0],rho_range[1],n_intervals,rhob1,rhob2,n_rootsPD,T,x,pcsaftCalcProps);

    auto n_rootsP = 0; // nº of roots of the Pressure , must be initialized as 0
    ArrayXr rho_extr; // Vector that contains density for extreme values
    rho_extr.resize(n_rootsPD+2);
    if(n_rootsPD>0)
    {
      // In this case the pressure curve has some extreme values
      // It can be a gas, liquid or a gas-liquid equilibrium depending on the
      // system pressure

      std::fill(rho_extr.begin(),rho_extr.end(),0.);
      for (int i = 0; i < n_rootsPD; i++)
      newton_bissectionD(rho_extr[i],rhob1[i],rhob2[i],T,x,pcsaftCalcProps);
        
        
      // Finding rho for which P is infinity (1e9 Pa)
      auto rho_aux = rho_extr[n_rootsPD-1];
      auto incre = rho_aux*0.1;
      pcsaft_calculation(pcsaftCalcProps,T,rho_aux,x);
      auto P_aux = pcsaftCalcProps.P;
      for(P_aux;P_aux<1.0e9;rho_aux+=incre)
      {
        pcsaft_calculation(pcsaftCalcProps,T,rho_aux,x);
        P_aux = pcsaftCalcProps.P;
      }
      rho_extr[n_rootsPD + 1] = rho_aux-incre; // It is the rho for which P > 1e9 Pa (Pinf)

      // Ordering rho_extr
      for (int i = n_rootsPD; i > 0  ; i--) rho_extr[i] = rho_extr[i-1];
      rho_extr[0] = 1.e-20;
    
      // rhob1 and rhob2 should be reinitialized as 0 
      std::fill(rhob1.begin(),rhob1.end(),0.);
      std::fill(rhob2.begin(),rhob2.end(),0.);
      
      // Obtaining rhob1, rhob2, and the number of roots for Pressure
      // Now, the Pressure density roots will be between rhob1[0] and rhob2[0], rhob1[1] and rhob2[1] ... and so on    
      real P1,P2;
      for (int i = 0; i < n_rootsPD + 1; i++)
      {
        pcsaft_calculation(pcsaftCalcProps,T,rho_extr[i],x);
        P1 = pcsaftCalcProps.P;
        pcsaft_calculation(pcsaftCalcProps,T,rho_extr[i+1],x);
        P2 = pcsaftCalcProps.P;      
        if (((P1 < P) && (P < P2)) || ((P1 > P) && (P > P2)))
        {
          rhob1[n_rootsP] = rho_extr[i];
          rhob2[n_rootsP] = rho_extr[i+1];
          n_rootsP += 1;
        }
      } 

    } else
    {
      // The Pressure density roots will be between rhob1[0] and rhob2[0]. 
      // In this case the pressure curve has no extreme values
      // It is a supercritical fluid
      zbrakP(rho_range[0],rho_range[1],n_intervals,rhob1,rhob2,n_rootsP,P,T,x,pcsaftCalcProps);
    }
    
    assert(("The number of density roots for pressure should be greater than 0 ",
      n_rootsP > 0));
    
    ArrayXr rho(n_rootsP);
    for (auto i = 0; i < n_rootsP; i++)
      newton_bissection(rho[i],rhob1[i],rhob2[i],P,T,x,pcsaftCalcProps);
    

    ArrayXr Gres(n_rootsP);
    for (auto i = 0 ; i < n_rootsP ; i++)
    {
      pcsaft_calculation(pcsaftCalcProps,T,rho[i],x);
      Gres[i] = pcsaftCalcProps.Gres; 
    }

    // get the element(position) of the Gres vector with minimum Gibbs Free Energy
    auto minElementIndex = std::min_element(Gres.begin(),Gres.end()) - Gres.begin();
    
    pcsaft_calculation(pcsaftCalcProps,T,rho[minElementIndex],x);    

    props.V = 1./rho[minElementIndex];
    auto rhoT = -pcsaftCalcProps.PT/pcsaftCalcProps.PD; // derivative of rho wrt T
    props.VT = -rhoT/(rho[minElementIndex]*rho[minElementIndex]);
    props.VP = 1./(-rho[minElementIndex]*rho[minElementIndex]*pcsaftCalcProps.PD);
    props.Gres = pcsaftCalcProps.Gres;
    props.Hres = pcsaftCalcProps.Hres;
    props.Cpres = pcsaftCalcProps.Cpres;
    props.Cvres = pcsaftCalcProps.Cvres;
    props.ln_phi = pcsaftCalcProps.ln_fugacity_coefficients;
    props.som = (n_rootsPD > 0) ? ((pcsaftCalcProps.PDD > 0.0) ? StateOfMatter::Liquid : 
      StateOfMatter::Gas) : StateOfMatter::Supercritical;

  } // end compute() calculation  

  auto zbrakPD(const real rho1, const real rho2, const int n,
    ArrayXr &rhob1, ArrayXr &rhob2, int &nroot,
    real T, ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void
  {  
    unsigned i;
    real rho,fp,fc,drho;
    auto nb=rhob1.size(); // maximum number of roots to be found for DP/Drho
    nroot=0;
    drho=(rho2-rho1)/n;
    pcsaft_calculation(pcsaftCalcProps,T,rho=rho1,x);
    fp = pcsaftCalcProps.PD;
    for (i=0;i<n;i++) 
    {
      pcsaft_calculation(pcsaftCalcProps,T,rho+=drho,x);    
      fc = pcsaftCalcProps.PD;
      if (fc*fp <= 0.0) 
      {
        rhob1[nroot]=rho-drho;
        rhob2[nroot++]=rho;
        if(nroot == nb) return;
      }
      fp=fc;
    }
  }

  auto zbrakP(const real rho1, const real rho2, const int n,
    ArrayXr &rhob1, ArrayXr &rhob2, int &nroot,
    real P, real T, ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void
  {  
    unsigned i;
    real rho,fp,fc,drho;
    int nb=rhob1.size(); // maximum number of roots to be found for P x rho
    nroot=0;
    drho=(rho2-rho1)/n;
    pcsaft_calculation(pcsaftCalcProps,T,rho=rho1,x);
    fp = pcsaftCalcProps.P/P - 1.;     
    for (i=0;i<n;i++) 
    {
      pcsaft_calculation(pcsaftCalcProps,T,rho+=drho,x);
	  fc = pcsaftCalcProps.P/P - 1.;        
      if (fc*fp <= 0.0)
      {
        rhob1[nroot]=rho-drho;
        rhob2[nroot++]=rho;
        if(nroot == nb) return;
      }
      fp=fc;
    }
  }

  auto newton_bissection(real &rho_c, real rho_a, real rho_b, real P, real T, 
    ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void
  {
    int nr, count = 1, max_ite = 100;
    auto a = rho_a, b = rho_b;
    real drho, fprime_c, f2prime_c , rho_new;
    pcsaft_calculation(pcsaftCalcProps,T,rho_a,x);
    auto fa = pcsaftCalcProps.P/P-1.;
    real fb;

    rho_c = (rho_a+rho_b)/2.;  
    pcsaft_calculation(pcsaftCalcProps,T,rho_c,x);
    auto fc = pcsaftCalcProps.P/P-1.;
    auto residue = fc;

    while ( (count < max_ite) && (abs(residue) > tol) )
    { 
      pcsaft_calculation(pcsaftCalcProps,T,rho_c,x);
      fprime_c = pcsaftCalcProps.PD/P;          
        
      if(!(abs(fprime_c) < 1.e-16))
      {
        drho = -fc/fprime_c; 
        rho_new = rho_c+drho;
      } 

      if (rho_new > rho_a && rho_new < rho_b)
      {
        rho_c = rho_new;
      } else 
      {
        if(fc*fa > 0.)
        {
          rho_a=rho_c;
          fa=fc; 
        } else 
        {
          rho_b=rho_c;
          fb=fc;
        }
        rho_c = (rho_a+rho_b)/2.;
      } 

      pcsaft_calculation(pcsaftCalcProps,T,rho_c,x); 
      fc = pcsaftCalcProps.P/P-1.;
      residue = fc;
      count++;
    }

    if (count == max_ite)
    {
      std::cout << "Calculation of density root values " << std::endl; 
      std::cout << "Density did not converge for P = " << P/1.e5 << " bar and T = " << T << " K" << std::endl;
      std::cout << "Residue = " << abs(residue) << std::endl;   
      std::cout << "Density between " << a << " mol/m³ and " << b << " mol/m³" << std::endl;
      exit(0);  
    }     

  }

  auto newton_bissectionD(real &rho_c, real rho_a, real rho_b, real T, 
    ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void
  {
    int nr, count = 1, max_ite = 100;
    auto a = rho_a, b = rho_b;
    real drho, fprime_c, rho_new;
    pcsaft_calculation(pcsaftCalcProps,T,rho_a,x);
    auto fa = pcsaftCalcProps.PD;
    real fb;

    rho_c = (rho_a+rho_b)/2.;  
    pcsaft_calculation(pcsaftCalcProps,T,rho_c,x);
    auto fc = pcsaftCalcProps.PD;
    auto residue = fc;

    while ( (count < max_ite) && (abs(residue) > tol) ) 
    { 
      pcsaft_calculation(pcsaftCalcProps,T,rho_c,x);  
      fprime_c = pcsaftCalcProps.PDD;      

      if(!(abs(fprime_c) < 1.e-16))
      {
        drho = -fc/fprime_c; 
        rho_new = rho_c+drho;
      } 

      if (rho_new > rho_a && rho_new < rho_b)
      {
        rho_c = rho_new;
      } else 
      {
        if(fc*fa > 0.)
        {
          rho_a=rho_c;
          fa=fc; 
        } else 
        {
          rho_b=rho_c;
          fb=fc;
        }
        rho_c = (rho_a+rho_b)/2.;
      } 
      
      pcsaft_calculation(pcsaftCalcProps,T,rho_c,x);
      fc = pcsaftCalcProps.PD;
      residue = fc;
      count++;
    }

    if (count == max_ite)
    {
      std::cout << "Calculation of extreme pressures values " << std::endl;
      std::cout << "Density did not converge for T = " << T << " K" << std::endl;
      std::cout << "Residue = " << abs(residue) << std::endl;   
      std::cout << "Density between " << a << " mol/m³ and " << b << " mol/m³" << std::endl;  
      exit(0);  
    }     

  }

  auto resizing_vectors(result_pcsaft_calculation &pcsaftCalcProps, int n_comp_New, int &n_comp_Old) -> void
  {
    if (n_comp_New > n_comp_Old)
    { 
      n_comp_Old = n_comp_New;
      auto order = n_comp_New*s;
      d.resize(n_comp_New);
      dDT.resize(n_comp_New);
      dDTT.resize(n_comp_New); // temperature-dependent diameter and its derivative wrt T
      p_dj.resize(n_comp_New);
      ghs.resize(n_comp_New*n_comp_New);
      ghsD.resize(n_comp_New*n_comp_New); // derivative of ghs wrt rho
      ghsDD.resize(n_comp_New*n_comp_New); // second derivative of ghs wrt rho
      ghsDDD.resize(n_comp_New*n_comp_New); // third derivative of ghs wrt rho
      ghsT.resize(n_comp_New*n_comp_New); // derivative of ghs wrt T  
      ghsTT.resize(n_comp_New*n_comp_New); // second derivative of ghs wrt T                                      
      Dghs.resize(n_comp_New*n_comp_New); //rho*Dghs_Drho
      DghsD.resize(n_comp_New*n_comp_New); // derivative of Dghs wrt rho
      DghsDD.resize(n_comp_New*n_comp_New); // second derivative of Dghs wrt rho
      DghsDDD.resize(n_comp_New*n_comp_New); // third derivative of Dghs wrt rho
      DghsT.resize(n_comp_New*n_comp_New); // derivative of Dghs wrt T
      DghsTT.resize(n_comp_New*n_comp_New); // second derivative of Dghs wrt T
      epsilon_ij.resize(n_comp_New*n_comp_New);
      sigma_ij.resize(n_comp_New*n_comp_New);
      mu_assoc.resize(n_comp_New);
      col.resize(order);  
      psi.resize(order);  
      dXA.resize(order);  // Auxiliary vector for calculating all XA derivatives
      Dghsii_Dx.resize(n_comp_New*n_comp_New);
      Dahs_Dx.resize(n_comp_New);
      Dadisp_Dx.resize(n_comp_New);
      Dahc_Dx.resize(n_comp_New);  
      mu_hc.resize(n_comp_New);
      mu_disp.resize(n_comp_New);
      mu.resize(n_comp_New);      
      S.resize(order);
      delta.resize(order*order);
      deltaD.resize(order*order);
      deltaDD.resize(order*order);
      deltaDDD.resize(order*order);
      deltaT.resize(order*order);
      deltaTT.resize(order*order); 
      delta_xm.resize(order*order*n_comp_New); 
      delta_xm_rho.resize(order*order*n_comp_New); 
      delta_xm_rho_rho.resize(order*order*n_comp_New);
      delta_xm_rho_rho_rho.resize(order*order*n_comp_New); 
      delta_xm_T.resize(order*order*n_comp_New); 
      delta_xm_T_T.resize(order*order*n_comp_New); 
      XA.resize(order);
      XA_old.resize(order);
      F.resize(order);        
      dif.resize(order);  
      indx.resize(order);    
      Fprime.resize(order*order);  
      lambda.resize(order*order);      
      inv_lambda.resize(order*order);
      XAT.resize(order);
      XATT.resize(order);
      XAD.resize(order);
      XADD.resize(order);
      XADDD.resize(order);
      psi2D.resize(order*n_comp_New); // Auxiliary vector for calculating XA derivatives wrt molar fraction
      XA_xm.resize(order*n_comp_New);
      XA_xm_T.resize(order*n_comp_New);
      XA_xm_T_T.resize(order*n_comp_New);
      XA_xm_rho.resize(order*n_comp_New);
      XA_xm_rho_rho.resize(order*n_comp_New);
      XA_xm_rho_rho_rho.resize(order*n_comp_New);  
      pcsaftCalcProps.ln_fugacity_coefficients.resize(n_comp_New);
    }
  }  

  auto calculate_XA(int n_comp, ArrayXrConstRef x,real num_den,Vec<int> &indx) -> void
  {
    // Creating and Calculating XA matrix by Newton's method
    auto count = 0;
    auto max_iter_XA = 150; // max nº of iterations for calculating XA, fraction of nonbonded sites 
    auto order = s*n_comp;
    std::fill(XA.begin(),XA.end(),0.5);
    int i,j,k,l,p,q;
    static real sum1,sum2,deltaljki,aux,norm0,norm,tolerance = tol*tol;

    // lambda function that represents sucessive substitution iterations
    auto suc_subs_iter = [&](int nssIter) -> void
    { 
      for(p = 0; p < nssIter; p++) 
      { 
        for(i = 0; i < order; i++) XA_old[i] = XA[i];

        for (i = 0; i < n_comp; i++) {
          for (j = 0; j < s; j++) {
            sum1 = 0.;
            for (k = 0; k < n_comp; k++){
              sum2 = 0.;
              for (l = 0; l < s; l++){
                if (j==l) continue;
                sum2 += S(l,k)*XA_old(l,k)*delta(l,k,j,i);
              }
              sum1 += x[k]*sum2;
            }
            XA(j,i) = 1./(1.+num_den*sum1);
          }
        }
      } 
    }; 

    // Initial Successive Substitution iterations
    suc_subs_iter(10);

    // Newton-Raphson Iterations
    do
    { 
      for(i = 0; i < order; i++) XA_old[i] = XA[i];

      // Creating F(X), the objetive function, and Fprime(X), the Jacobian Matrix  
      p = -1;
      q = -1;    
      for (i = 0; i < n_comp; i++) 
      {
        for (j = 0; j < s; j++) 
        {
          ++p;
          sum1 = 0.;
          for (k = 0; k < n_comp; k++)
          { 
            sum2 = 0.;
            for (l = 0; l < s; l++)
            { 
              ++q;
              if ((j==l) && (i==k)) deltaljki = 1.;
              else deltaljki = 0.;

              sum2 += S(l,k)*XA_old(l,k)*delta(l,k,j,i);
              Fprime(p,q) = deltaljki/XA_old(j,i) + XA_old(j,i)*num_den*x[k]*S(l,k)*
                delta(l,k,j,i);       
            }
            sum1 += x[k]*sum2;
          }
          F[p] = 1.-XA_old(j,i)*(1.+num_den*sum1);
          q = -1;
        }
      }     

      // LU decomposition and Back Substitution to solve the linear system:
      // Fprime * deltaXA = -F
      detail::ludcmp(Fprime,order,indx,aux);
      detail::lubksb(Fprime,order,indx,F);

      // Updating X
      p = -1;
      norm = 0.;
      for (i = 0; i < n_comp; i++) 
      {
        for (j = 0; j < s; j++)  
        { 
          XA(j,i) = XA_old(j,i) + F[++p];
          norm += F[p]*F[p];
        }
      }  

      // Successive Substitution iterations if Newton does not converge 
      // count > 0 in order to not compare initial Succe. Subs. to Newton
      if ((norm > norm0) && (count > 0)) suc_subs_iter(5);

      norm0 = norm;
      count++;

    } while ((count < max_iter_XA) && (norm > tol));

    if (count == max_iter_XA) 
      std::cout << "XA calculation did not converge" << std::endl;    

  }  
  
  auto pcsaft_calculation(result_pcsaft_calculation &pcsaftCalcProps, real T, real rho, ArrayXrConstRef x) -> void
  {
  // Input 
  //     T: Temperature (K) 
  //     rho: density (mol/m³)  
  //     x: mole fractions of species (vector) 
  //
  //
  // Output
  //     pcsaftCalcProps: Struct containing all thermodinamic properties
  //
  // --------------------------------------------------------------------------------  
  
  auto n_comp = nspecies; // number of components  
  
  //  Calculating d and its derivatives
  real aux;
  for (int i = 0; i < n_comp; i++){
    aux = exp(-3.*epsilon(i)/T);
    d[i] = sigma(i)*(1.-0.12*aux); // EQ. 3
    dDT[i] = -9.*epsilon(i)*sigma(i)/25./T/T*aux; // derivative of d[i] wrt T
    dDTT[i] = 9.*(2.*T-3.*epsilon(i))*epsilon(i)*sigma(i)/25./T/T/T/T*aux; // second derivative of d[i] wrt T
  }

  //  Calculating numerical density and its derivatives
  auto num_den = rho*n_av/1.e30; // numerical density, part/angs³
  auto num_denD = num_den/rho; // derivative of num_den wrt rho
  auto num_denDD = 0.; // second derivative of num_den wrt rho
  auto num_denDDD = 0.; // third derivative of num_den wrt rho

  //  Calculating zeta and its derivatives
  static real zeta[4],zetaD[4],zetaDD[4],zetaDDD[4],zetaT[4],zetaTT[4];
  real summation,summationT,summationTT;
  std::fill(p_dj.begin(),p_dj.end(),1.);
  for (int i = 0; i < 4; i++){
    summation = 0.;
    summationT = 0.;
    summationTT = 0.;
    for (int j = 0; j < n_comp; j++){
      aux = x[j]*m(j)*p_dj[j];
      summation += aux; 
      summationT += aux/d[j]*dDT[j];
      summationTT +=  i*x[j]*m(j)*dDT[j]*p_dj[j]/d[j]/d[j]*dDT[j]*(i-1.) +
        i*x[j]*m(j)*p_dj[j]/d[j]*dDTT[j];
      p_dj[j] *= d[j];   
    }
    zeta[i] = pi/6.*num_den*summation; // EQ. 9 , zetas are functions of T and rho
    zetaD[i] =  zeta[i]/rho; // derivative of zeta[i] wrt rho
    zetaDD[i] = 0.; // second derivative of zeta[i] wrt rho
    zetaDDD[i] = 0.; // third derivative of zeta[i] wrt rho
    zetaT[i] =  pi/6.*num_den*i*summationT;  // derivative of zeta[i] wrt T
    zetaTT[i] = pi/6.*num_den*summationTT; // second derivative of zeta[i] wrt T
  }

  //  Calculating eta and its derivatives
  auto eta = zeta[3]; // packing fraction, below EQ. 13
  auto etaD = zetaD[3]; // derivative of eta wrt rho
  auto etaDD = zetaDD[3]; // second derivative of eta wrt rho
  auto etaDDD = zetaDDD[3]; // third derivative of eta wrt rho
  auto etaT = zetaT[3]; // derivative of eta wrt T
  auto etaTT = zetaTT[3]; // second derivative of eta wrt T

  real m_avg = 0.;
  for (int i = 0; i < n_comp; i++){
    m_avg += x[i]*m(i); // EQ. 6
  }  


  // Auxiliary variables to calculate Zhs and its derivatives
  static real auxD,auxDD,auxDDD,aux21,aux21D,aux21DD,aux21DDD,aux22,aux22D,aux22DD,aux22DDD,
    aux2,aux2D,aux2DD,aux2DDD,aux31,aux31D,aux31DD,aux31DDD,aux32,aux32D,aux32DD,aux32DDD,
    aux3,aux3D,aux3DD,aux3DDD;
  aux = zeta[3]/(1.-zeta[3]);
  auxD = zetaD[3]/(1. - zeta[3]) + (zeta[3]*zetaD[3])/(1. - zeta[3])/(1. - zeta[3]);
  auxDD = (-2.*(zetaD[3]*zetaD[3]) + (-1. + 1.*zeta[3])*zetaDD[3])/(-1. + zeta[3])/(-1. + zeta[3])/(-1. + zeta[3]);//2*auxD*auxD/(1.+aux)+zetaDD[3]*aux/zetaD[3];
  auxDDD = (6.*((zetaD[3])*(zetaD[3])*(zetaD[3])) + (6. - 6.*zeta[3])*zetaD[3]*
    zetaDD[3] + ((1. - 1.*zeta[3])*(1. - 1.*zeta[3]))*zetaDDD[3])/
    ((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]));
  aux21 = 3.*zeta[1]/zeta[0];
  aux21D = 0.; // It's identically zero
  aux21DD = 0.;
  aux21DDD = 0.; 
  aux22 = zeta[2]/(1.-zeta[3])/(1.-zeta[3]);
  aux22D = zetaD[2]/(1. - zeta[3])/(1. - zeta[3]) + (2*zeta[2]*zetaD[3])/(1. - zeta[3])/(1. - zeta[3])/(1. - zeta[3]);
  aux22DD = ((4. - 4.*zeta[3])*zetaD[2]*zetaD[3] +(1. - 2.*zeta[3] + (zeta[3]*zeta[3]))*zetaDD[2] + 
    zeta[2]*(6.*(zetaD[3]*zetaD[3]) + (2. - 2.*zeta[3])*zetaDD[3]))/((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]));//aux22D*aux22D/aux22 + aux22*(-zetaD[2]*zetaD[2]/zeta[2]/zeta[2]+zetaDD[2]/zeta[2] + 2.*zetaD[3]*zetaD[3]/(1.-zeta[3])/(1.-zeta[3])+2.*zetaDD[3]/(1.-zeta[3]));
  aux22DDD = (-1.*(6.*((-1. + zeta[3])*(-1. + zeta[3]))*zetaD[3]*zetaDD[2] + 
    6.*(-1. + zeta[3])*zetaD[2]*(-3.*(zetaD[3]*zetaD[3]) + (-1. + zeta[3])*zetaDD[3]) - 
    1.*((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]))*zetaDDD[2] + 
    2*zeta[2]*(12*(zetaD[3]*zetaD[3]*zetaD[3]) - 9.*(-1. + zeta[3])*zetaD[3]*zetaDD[3] + 
    ((-1. + zeta[3])*(-1. + zeta[3]))*zetaDDD[3])))/
    ((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]));
  aux2 = aux21*aux22;
  aux2D = aux21D*aux22+aux21*aux22D;
  aux2DD = aux21DD*aux22+aux21D*aux22D+aux21D*aux22D+aux21*aux22DD;
  aux2DDD = aux21DDD*aux22+aux21DD*aux22D+aux21DD*aux22D+aux21D*aux22DD+aux21DD*aux22D+aux21D*aux22DD+aux21D*aux22DD+aux21*aux22DDD;   
  aux31 = zeta[2]*zeta[2]*zeta[2]/zeta[0];
  aux31D = (zeta[2]*zeta[2]*(-zeta[2]*zetaD[0] + 3*zeta[0]*zetaD[2]))/zeta[0]/zeta[0];
  aux31DD = (zeta[2]*(6*zeta[0]*zeta[0]*(zetaD[2]*zetaD[2]) + (zeta[2]*zeta[2])*(2*zetaD[0]*zetaD[0] - zeta[0]*zetaDD[0]) + 
    3.*zeta[0]*zeta[2]*(-2*zetaD[0]*zetaD[2] + zeta[0]*zetaDD[2])))/(zeta[0]*zeta[0]*zeta[0]);//aux31D*aux31D/aux31 + aux31*(zetaD[0]*zetaD[0]/zeta[0]/zeta[0] - zetaDD[0]/zeta[0] - 3.*zetaD[2]*zetaD[2]/(zeta[2])/(zeta[2])+3.*zetaDD[2]/(zeta[2]));; 
  aux31DDD = (6*(zeta[0]*zeta[0]*zeta[0])*((zetaD[2])*(zetaD[2])*(zetaD[2])) + 
    18*((zeta[0])*(zeta[0]))*zeta[2]*zetaD[2]*(-(zetaD[0]*zetaD[2]) + zeta[0]*zetaDD[2]) - 
    ((zeta[2])*(zeta[2])*(zeta[2]))*(6*((zetaD[0])*(zetaD[0])*(zetaD[0])) - 6*zeta[0]*zetaD[0]*zetaDD[0] + 
    ((zeta[0])*(zeta[0]))*zetaDDD[0]) + 3*zeta[0]*((zeta[2])*(zeta[2]))*(6*((zetaD[0])*(zetaD[0]))*zetaD[2] - 
    3*zeta[0]*zetaD[0]*zetaDD[2] + zeta[0]*(-3*zetaD[2]*zetaDD[0] + zeta[0]*zetaDDD[2])))/((zeta[0])*(zeta[0])*(zeta[0])*(zeta[0]));
  aux32 = (3.-zeta[3])/(1.-zeta[3])/(1.-zeta[3])/(1.-zeta[3]);
  aux32D = ((8.- 2.*zeta[3])*zetaD[3])/(-1. + zeta[3])/(-1. + zeta[3])/(-1. + zeta[3])/(-1. + zeta[3]);
  aux32DD = ((-30. + 6.*zeta[3])*(zetaD[3]*zetaD[3]) + (-8. + 10.*zeta[3] - 2.*(zeta[3]*zeta[3]))*zetaDD[3])/((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]));//aux32D*(aux32D/aux32+zetaD[3]/(1.-zeta[3])+zetaDD[3]/zetaD[3]);
  aux32DDD = ((144. - 24.*zeta[3])*((zetaD[3])*(zetaD[3])*(zetaD[3])) + 
    (90. - 108.*zeta[3] + 18.*((zeta[3])*(zeta[3])))*zetaD[3]*zetaDD[3] - 
    2.*((1. - 1.*zeta[3])*(1. - 1.*zeta[3]))*(-4. + 1.*zeta[3])*zetaDDD[3])/((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]));
  aux3 = aux31*aux32;
  aux3D = aux31D*aux32+aux31*aux32D;
  aux3DD = aux31DD*aux32+aux31D*aux32D+aux31D*aux32D+aux31*aux32DD;
  aux3DDD = aux31DDD*aux32+aux31DD*aux32D+aux31DD*aux32D+aux31D*aux32DD+aux31DD*aux32D+aux31D*aux32DD+aux31D*aux32DD+aux31*aux32DDD;
 
  static real auxT,auxTT,aux21T,aux21TT,aux22T,aux22TT,
    aux2T,aux2TT,aux31T,aux31TT,aux32T,aux32TT,aux3T,aux3TT;
  auxT = zetaT[3]/(1. - zeta[3]) + (zeta[3]*zetaT[3])/(1. - zeta[3])/(1. - zeta[3]);
  auxTT = (-2.*(zetaT[3]*zetaT[3]) + (-1. + 1.*zeta[3])*zetaTT[3])/(-1. + zeta[3])/(-1. + zeta[3])/(-1. + zeta[3]);//2*auxT*auxT/(1.+aux)+zetaTT[3]*aux/zetaT[3];
  aux21T = (3.*(-zeta[1]*zetaT[0] + zeta[0]*zetaT[1]))/zeta[0]/zeta[0];
  aux21TT = (6*zeta[1]*zetaT[0]*zetaT[0] - 6*zeta[0]*zetaT[0]*zetaT[1] - 
    3*zeta[0]*zeta[1]*zetaTT[0] + 3*zeta[0]*zeta[0]*zetaTT[1])/(zeta[0]*zeta[0]*zeta[0]);//(-2.*aux21T*zetaT[0]-aux21*zetaTT[0]+3.*zetaTT[1])/zeta[0];
  aux22T = zetaT[2]/(1. - zeta[3])/(1. - zeta[3]) + (2*zeta[2]*zetaT[3])/(1. - zeta[3])/(1. - zeta[3])/(1. - zeta[3]);
  aux22TT = ((4. - 4.*zeta[3])*zetaT[2]*zetaT[3] +(1. - 2.*zeta[3] + (zeta[3]*zeta[3]))*zetaTT[2] + 
    zeta[2]*(6.*(zetaT[3]*zetaT[3]) + (2. - 2.*zeta[3])*zetaTT[3]))/((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]));//aux22T*aux22T/aux22 + aux22*(-zetaT[2]*zetaT[2]/zeta[2]/zeta[2]+zetaTT[2]/zeta[2] + 2.*zetaT[3]*zetaT[3]/(1.-zeta[3])/(1.-zeta[3])+2.*zetaTT[3]/(1.-zeta[3]));
  aux2T = aux21T*aux22+aux21*aux22T;
  aux2TT = aux21TT*aux22+aux21T*aux22T+aux21T*aux22T+aux21*aux22TT;
  aux31T = (zeta[2]*zeta[2]*(-zeta[2]*zetaT[0] + 3*zeta[0]*zetaT[2]))/zeta[0]/zeta[0];
  aux31TT = (zeta[2]*(6*zeta[0]*zeta[0]*(zetaT[2]*zetaT[2]) + (zeta[2]*zeta[2])*(2*zetaT[0]*zetaT[0] - zeta[0]*zetaTT[0]) + 
    3*zeta[0]*zeta[2]*(-2*zetaT[0]*zetaT[2] + zeta[0]*zetaTT[2])))/(zeta[0]*zeta[0]*zeta[0]);//aux31T*aux31T/aux31 + aux31*(zetaT[0]*zetaT[0]/zeta[0]/zeta[0] - zetaTT[0]/zeta[0] - 3.*zetaT[2]*zetaT[2]/(zeta[2])/(zeta[2])+3.*zetaTT[2]/(zeta[2]));; 
  aux32T = ((8.- 2.*zeta[3])*zetaT[3])/(-1. + zeta[3])/(-1. + zeta[3])/(-1. + zeta[3])/(-1. + zeta[3]);
  aux32TT = ((-30. + 6.*zeta[3])*(zetaT[3]*zetaT[3]) + (-8. + 10.*zeta[3] - 2.*(zeta[3]*zeta[3]))*zetaTT[3])/((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]));//aux32T*(aux32T/aux32+zetaT[3]/(1.-zeta[3])+zetaTT[3]/zetaT[3]);
  aux3T = aux31T*aux32+aux31*aux32T;
  aux3TT = aux31TT*aux32+aux31T*aux32T+aux31T*aux32T+aux31*aux32TT;  
  
  auto Zhs = aux+aux2+aux3;
  auto ZhsD = auxD+aux2D+aux3D;
  auto ZhsDD = auxDD+aux2DD+aux3DD;
  auto ZhsDDD = auxDDD+aux2DDD+aux3DDD;
  auto ZhsT = auxT+aux2T+aux3T; 
  auto ZhsTT = auxTT+aux2TT+aux3TT;
  
  real m2epsilonsigma3 = 0.;
  real m2epsilonsigma3T = 0.; // derivative of m2epsilonsigma3 wrt T
  real m2epsilonsigma3TT = 0.; // second derivative of m2epsilonsigma3 wrt T
  real m2epsilon2sigma3 = 0.; 
  real m2epsilon2sigma3T = 0.; // derivative of m2epsilon2sigma3 wrt T
  real m2epsilon2sigma3TT = 0.; // derivative of m2epsilon2sigma3 wrt T

  // Calculation of ghs and its derivatives,
  // Calculation of Dghs and its derivatives, 
  // Calculation of m2epsilonsigma3 and m2epsilon2sigma3 (and its derivatives),
  // Calculation of epsilon_ij and sigma_ij
  int idx = 0;
  for (int i = 0; i < n_comp; i++)
  {
    for (int j = 0; j < n_comp; j++)
    {
      // Calculation of auxiliary variables to obtain ghs and its derivatives
      aux = 1./(1.-zeta[3]);
      aux21 = (d[i]*d[j]/(d[i]+d[j]));
      aux22 = 3.*zeta[2]/(1.-zeta[3])/(1.-zeta[3]);
      aux2 = aux21*aux22;
      aux31 = (d[i]*d[j]/(d[i]+d[j]))*(d[i]*d[j]/(d[i]+d[j]));
      aux32 = 2.*zeta[2]*zeta[2]/(1.-zeta[3])/(1.-zeta[3])/(1.-zeta[3]);
      aux3 = aux31*aux32;
      auxD = zetaD[3]/((-1 + zeta[3])*(-1 + zeta[3]));
      auxDD = (-2*(zetaD[3]*zetaD[3]) + (-1 + zeta[3])*zetaDD[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      auxDDD = (6*((zetaD[3])*(zetaD[3])*(zetaD[3])) - 6*(-1 + zeta[3])*zetaD[3]*
        zetaDD[3] + ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDDD[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux21D = 0.;
      aux21DD = 0.;
      aux21DDD = 0;
      aux22D = (3*(-1 + zeta[3])*zetaD[2] - 6*zeta[2]*zetaD[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux22DD = (3*(-4*(-1 + zeta[3])*zetaD[2]*zetaD[3] + 6*zeta[2]*(zetaD[3]*zetaD[3]) + ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDD[2] - 
        2*zeta[2]*(-1 + zeta[3])*zetaDD[3]))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux22DDD = (-3*(-6*(-1 + zeta[3])*zetaD[2]*(3*((zetaD[3])*(zetaD[3])) - (-1 + zeta[3])*zetaDD[3]) + 
        ((-1 + zeta[3])*(-1 + zeta[3]))*(6*zetaD[3]*zetaDD[2] - (-1 + zeta[3])*zetaDDD[2]) + 
        2*zeta[2]*(12*((zetaD[3])*(zetaD[3])*(zetaD[3])) - 9*(-1 + zeta[3])*zetaD[3]*zetaDD[3] + 
        ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDDD[3])))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux31D = 0.;
      aux31DD = 0.;
      aux31DDD = 0;
      aux32D = (zeta[2]*((4. - 4.*zeta[3])*zetaD[2] + 6.*zeta[2]*zetaD[3]))/((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]));
      aux32DD = (-4.*((1. - zeta[3])*(1. - zeta[3]))*(zetaD[2]*zetaD[2]) + zeta[2]*(-24. + 24.*zeta[3])*zetaD[2]*zetaD[3] + 
        zeta[2]*(-4.*((1. - zeta[3])*(1. - zeta[3]))*zetaDD[2] + zeta[2]*(-24.*(zetaD[3]*zetaD[3]) + (-6. + 6.*zeta[3])*zetaDD[3])))
        /((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));        
      aux32DDD = (2*(18.*((-1. + zeta[3])*(-1. + zeta[3]))*zetaD[3]*((zetaD[2]*zetaD[2]) + zeta[2]*zetaDD[2]) + 
        18.*zeta[2]*(-1. + zeta[3])*zetaD[2]*(-4.*(zetaD[3]*zetaD[3]) + (-1. + zeta[3])*zetaDD[3]) - 
        2.*((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]))*(3.*zetaD[2]*zetaDD[2] + 
        zeta[2]*zetaDDD[2]) + 3*(zeta[2]*zeta[2])*(20*(zetaD[3]*zetaD[3]*zetaD[3]) - 
        12.*(-1. + zeta[3])*zetaD[3]*zetaDD[3] + ((-1. + zeta[3])*(-1. + zeta[3]))*zetaDDD[3])))/((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]));;
      aux2D = aux21D*aux22+aux21*aux22D;
      aux2DD = aux21DD*aux22+aux21D*aux22D+aux21D*aux22D+aux21*aux22DD; 
      aux2DDD = aux21DDD*aux22+aux21DD*aux22D+aux21DD*aux22D+aux21D*aux22DD+aux21DD*aux22D+aux21D*aux22DD+aux21D*aux22DD+aux21*aux22DDD;
      aux3D = aux31D*aux32+aux31*aux32D;
      aux3DD = aux31DD*aux32+aux31D*aux32D+aux31D*aux32D+aux31*aux32DD;
      aux3DDD = aux31DDD*aux32+aux31DD*aux32D+aux31DD*aux32D+aux31D*aux32DD+aux31DD*aux32D+aux31D*aux32DD+aux31D*aux32DD+aux31*aux32DDD;
        

      auxT = zetaT[3]/((-1 + zeta[3])*(-1 + zeta[3]));
      auxTT = (-2*(zetaT[3]*zetaT[3]) + (-1 + zeta[3])*zetaTT[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux21T = ((d[j]*d[j])*dDT[i] + (d[i]*d[i])*dDT[j])/((d[i] + d[j])*(d[i] + d[j]));
      aux21TT = ((d[j]*d[j]*d[j])*dDTT[i] + (d[j]*d[j])*(-2*(dDT[i]*dDT[i]) + d[i]*dDTT[i]) + 
        d[i]*d[j]*(4*dDT[i]*dDT[j] + d[i]*dDTT[j]) + (d[i]*d[i])*(-2*(dDT[j]*dDT[j]) + d[i]*dDTT[j]))/
        ((d[i] + d[j])*(d[i] + d[j])*(d[i] + d[j]));
      aux22T = (3*(-1 + zeta[3])*zetaT[2] - 6*zeta[2]*zetaT[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])); 
      aux22TT = (3*(-4*(-1 + zeta[3])*zetaT[2]*zetaT[3] + 6*zeta[2]*(zetaT[3]*zetaT[3]) + ((-1 + zeta[3])*(-1 + zeta[3]))*zetaTT[2] - 
        2*zeta[2]*(-1 + zeta[3])*zetaTT[3]))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux31T = (2*(d[i]*(d[j]*d[j]*d[j])*dDT[i] + (d[i]*d[i]*d[i])*d[j]*dDT[j]))/((d[i] + d[j])*(d[i] + d[j])*(d[i] + d[j])); 
      aux31TT = (2*((d[j]*d[j]*d[j]*d[j])*(dDT[i]*dDT[i]) + d[i]*(d[j]*d[j]*d[j])*(-2*(dDT[i]*dDT[i]) + d[j]*dDTT[i]) + 
        (d[i]*d[i])*(d[j]*d[j])*(6*dDT[i]*dDT[j] + d[j]*dDTT[i]) + (d[i]*d[i]*d[i])*d[j]*
        (-2*(dDT[j]*dDT[j]) + d[j]*dDTT[j]) + (d[i]*d[i]*d[i]*d[i])*((dDT[j]*dDT[j]) + d[j]*dDTT[j])))/
        ((d[i] + d[j])*(d[i] + d[j])*(d[i] + d[j])*(d[i] + d[j]));
      aux32T = (zeta[2]*((4. - 4.*zeta[3])*zetaT[2] + 6.*zeta[2]*zetaT[3]))/((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])); 
      aux32TT = (-4.*((1. - zeta[3])*(1. - zeta[3]))*(zetaT[2]*zetaT[2]) + zeta[2]*(-24. + 24.*zeta[3])*zetaT[2]*zetaT[3] + 
        zeta[2]*(-4.*((1. - zeta[3])*(1. - zeta[3]))*zetaTT[2] + zeta[2]*(-24.*(zetaT[3]*zetaT[3]) + (-6. + 6.*zeta[3])*zetaTT[3])))
        /((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));        
      aux2T = aux21T*aux22+aux21*aux22T;
      aux2TT = aux21TT*aux22+aux21T*aux22T+aux21T*aux22T+aux21*aux22TT;
      aux3T = aux31T*aux32+aux31*aux32T;
      aux3TT = aux31TT*aux32+aux31T*aux32T+aux31T*aux32T+aux31*aux32TT;  

      ghs[idx] = aux+aux2+aux3;      
      ghsD[idx] = auxD+aux2D+aux3D; 
      ghsDD[idx] = auxDD+aux2DD+aux3DD; 
      ghsDDD[idx] = auxDDD+aux2DDD+aux3DDD; 
      ghsT[idx] = auxT+aux2T+aux3T; 
      ghsTT[idx] = auxTT+aux2TT+aux3TT; 
      

      // Calculation of auxiliary variables to obtain Dghs and its derivatives
      aux = zeta[3]/(1.-zeta[3])/(1.-zeta[3]);
      aux21 = (d[i]*d[j]/(d[i]+d[j]));
      aux22 = 3.*zeta[2]/(1.-zeta[3])/(1.-zeta[3]) + 6.*zeta[2]*zeta[3]/(1.-zeta[3])/(1.-zeta[3])/(1.-zeta[3]);
      aux2 = aux21*aux22;
      aux31 = (d[i]*d[j]/(d[i]+d[j]))*(d[i]*d[j]/(d[i]+d[j]));
      aux32 = 4.*zeta[2]*zeta[2]/(1.-zeta[3])/(1.-zeta[3])/(1.-zeta[3])+6.*zeta[2]*zeta[2]*zeta[3]/(1.-zeta[3])/(1.-zeta[3])/(1.-zeta[3])/(1.-zeta[3]);
      aux3 = aux31*aux32;
      auxD = -(((1 + zeta[3])*zetaD[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])));
      auxDD = (2*(2 + zeta[3])*(zetaD[3]*zetaD[3]) - (-1 + (zeta[3]*zeta[3]))*zetaDD[3])/
        ((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      auxDDD = (-6*(3 + zeta[3])*((zetaD[3])*(zetaD[3])*(zetaD[3])) + 6*(-2 + zeta[3] + ((zeta[3])*(zeta[3])))*zetaD[3]*zetaDD[3] - 
        ((-1 + zeta[3])*(-1 + zeta[3]))*(1 + zeta[3])*zetaDDD[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux21D = 0.;
      aux21DD = 0.;
      aux21DDD = 0;
      aux22D = (-3*(-1 + (zeta[3]*zeta[3]))*zetaD[2] + 6*zeta[2]*(2 + zeta[3])*zetaD[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux22DD = (12*(-2 + zeta[3] + (zeta[3]*zeta[3]))*zetaD[2]*zetaD[3] - 3*((-1 + zeta[3])*(-1 + zeta[3]))*(1 + zeta[3])*zetaDD[2] + 
        6*zeta[2]*(-3*(3 + zeta[3])*(zetaD[3]*zetaD[3]) + (-2 + zeta[3] + (zeta[3]*zeta[3]))*zetaDD[3]))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux22DDD = (3*(-6*(-1 + zeta[3])*zetaD[2]*(3*(3 + zeta[3])*((zetaD[3])*(zetaD[3])) - 
        (-2 + zeta[3] + ((zeta[3])*(zeta[3])))*zetaDD[3]) + ((-1 + zeta[3])*(-1 + zeta[3]))*(6*(2 + zeta[3])*zetaD[3]*zetaDD[2] - 
        (-1 + ((zeta[3])*(zeta[3])))*zetaDDD[2]) + 2*zeta[2]*(12*(4 + zeta[3])*((zetaD[3])*(zetaD[3])*(zetaD[3])) - 
        9*(-3 + 2*zeta[3] + ((zeta[3])*(zeta[3])))*zetaD[3]*zetaDD[3] + (2 - 3*zeta[3] + ((zeta[3])*(zeta[3])*(zeta[3])))*zetaDDD[3])))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux31D = 0.;
      aux31DD = 0.;
      aux31DDD = 0.;
      aux32D = (4*zeta[2]*(-2 + zeta[3] + (zeta[3]*zeta[3]))*zetaD[2] - 6*(zeta[2]*zeta[2])*(3 + zeta[3])*zetaD[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux32DD = (4*(2 - 3*zeta[3] + (zeta[3]*zeta[3]*zeta[3]))*(zetaD[2]*zetaD[2]) - 
        24*zeta[2]*(-3 + 2*zeta[3] + (zeta[3]*zeta[3]))*zetaD[2]*zetaD[3] + 4*zeta[2]*(2 - 3*zeta[3] + (zeta[3]*zeta[3]*zeta[3]))*zetaDD[2] + 6*(zeta[2]*zeta[2])*
        (4*(4 + zeta[3])*(zetaD[3]*zetaD[3]) - (-3 + 2*zeta[3] + (zeta[3]*zeta[3]))*zetaDD[3]))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));        
      aux32DDD = (-2*(18*((-1 + zeta[3])*(-1 + zeta[3]))*(3 + zeta[3])*((zetaD[2])*(zetaD[2]))*
        zetaD[3] - 6*(-1 + zeta[3])*zetaD[2]*((2 - 3*zeta[3] + ((zeta[3])*(zeta[3])*(zeta[3])))*zetaDD[2] + 
        3*zeta[2]*(4*(4 + zeta[3])*((zetaD[3])*(zetaD[3])) - (-3 + 2*zeta[3] + ((zeta[3])*(zeta[3])))*zetaDD[3])) + 
        zeta[2]*(-2*((-1 + zeta[3])*(-1 + zeta[3]))*(-9*(3 + zeta[3])*zetaD[3]*zetaDD[2] + 
        (-2 + zeta[3] + ((zeta[3])*(zeta[3])))*zetaDDD[2]) + 3*zeta[2]*(20*(5 + zeta[3])*((zetaD[3])*(zetaD[3])*(zetaD[3])) - 
        12*(-4 + 3*zeta[3] + ((zeta[3])*(zeta[3])))*zetaD[3]*zetaDD[3] + ((-1 + zeta[3])*(-1 + zeta[3]))*(3 + zeta[3])*zetaDDD[3]))))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));   
      aux2D = aux21D*aux22+aux21*aux22D;
      aux2DD = aux21DD*aux22+aux21D*aux22D+aux21D*aux22D+aux21*aux22DD; 
      aux2DDD = aux21DDD*aux22+aux21DD*aux22D+aux21DD*aux22D+aux21D*aux22DD+aux21DD*aux22D+aux21D*aux22DD+aux21D*aux22DD+aux21*aux22DDD;
      aux3D = aux31D*aux32+aux31*aux32D;
      aux3DD = aux31DD*aux32+aux31D*aux32D+aux31D*aux32D+aux31*aux32DD; 
      aux3DDD = aux31DDD*aux32+aux31DD*aux32D+aux31DD*aux32D+aux31D*aux32DD+aux31DD*aux32D+aux31D*aux32DD+aux31D*aux32DD+aux31*aux32DDD;

      auxT = -(((1 + zeta[3])*zetaT[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])));
      auxTT = (2*(2 + zeta[3])*(zetaT[3]*zetaT[3]) - (-1 + (zeta[3]*zeta[3]))*zetaTT[3])/
        ((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      //aux21T,  the same as in ghs calculations
      //aux21TT, the same as in ghs calculations
      aux22T = (-3*(-1 + (zeta[3]*zeta[3]))*zetaT[2] + 6*zeta[2]*(2 + zeta[3])*zetaT[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux22TT = (12*(-2 + zeta[3] + (zeta[3]*zeta[3]))*zetaT[2]*zetaT[3] - 3*((-1 + zeta[3])*(-1 + zeta[3]))*(1 + zeta[3])*zetaTT[2] + 
        6*zeta[2]*(-3*(3 + zeta[3])*(zetaT[3]*zetaT[3]) + (-2 + zeta[3] + (zeta[3]*zeta[3]))*zetaTT[3]))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      //aux31T, the same as in ghs calculations
      //aux31TT, the same as in ghs calculations
      aux32T = (4*zeta[2]*(-2 + zeta[3] + (zeta[3]*zeta[3]))*zetaT[2] - 6*(zeta[2]*zeta[2])*(3 + zeta[3])*zetaT[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
      aux32TT = (4*(2 - 3*zeta[3] + (zeta[3]*zeta[3]*zeta[3]))*(zetaT[2]*zetaT[2]) - 
        24*zeta[2]*(-3 + 2*zeta[3] + (zeta[3]*zeta[3]))*zetaT[2]*zetaT[3] + 4*zeta[2]*(2 - 3*zeta[3] + (zeta[3]*zeta[3]*zeta[3]))*
        zetaTT[2] + 6*(zeta[2]*zeta[2])*(4*(4 + zeta[3])*(zetaT[3]*zetaT[3]) - 
        (-3 + 2*zeta[3] + (zeta[3]*zeta[3]))*zetaTT[3]))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));          
      aux2T = aux21T*aux22+aux21*aux22T;
      aux2TT = aux21TT*aux22+aux21T*aux22T+aux21T*aux22T+aux21*aux22TT;
      aux3T = aux31T*aux32+aux31*aux32T;
      aux3TT = aux31TT*aux32+aux31T*aux32T+aux31T*aux32T+aux31*aux32TT;

      Dghs[idx] = aux+aux2+aux3; 
      DghsD[idx] = auxD+aux2D+aux3D; 
      DghsDD[idx] = auxDD+aux2DD+aux3DD;
      DghsDDD[idx] = auxDDD+aux2DDD+aux3DDD; 
      DghsT[idx] = auxT+aux2T+aux3T; 
      DghsTT[idx] = auxTT+aux2TT+aux3TT; 

     // Calculation of m2epsilonsigma3 and m2epsilon2sigma3 (and its derivatives),
     // Calculation of epsilon_ij and sigma_ij      
      sigma_ij[idx] = (sigma(i) + sigma(j))/2.; // EQ. 23 
      epsilon_ij[idx] = sqrt(epsilon(i) * epsilon(j)) * (1. - kij(idx)); // EQ. 24      
      aux = x[i]*x[j]*m(i)*m(j)*epsilon_ij[idx]*(sigma_ij[idx]*sigma_ij[idx]*sigma_ij[idx])/T;
      m2epsilonsigma3 +=  aux;// EQ. A12 
      m2epsilonsigma3T -= aux/T; 
      m2epsilonsigma3TT += 2.*aux/T/T;  
      aux = x[i]*x[j]*m(i)*m(j)*((epsilon_ij[idx]/T)*(epsilon_ij[idx]/T))*(sigma_ij[idx]*sigma_ij[idx]*sigma_ij[idx]);       
      m2epsilon2sigma3 +=  aux;// EQ. A13
      m2epsilon2sigma3T -= 2.*aux/T; 
      m2epsilon2sigma3TT += 6.*aux/T/T;         
      
      idx += 1;  
    }
  }

  // Calculation of a and b vectors
  static real a[7]; 
  static real b[7]; 
  for (int i = 0; i < 7; i++){
    a[i] = a0[i] + (m_avg-1.)/m_avg*a1[i] + (m_avg-1.)/m_avg*(m_avg-2.)/m_avg*a2[i]; // EQ. A18
    b[i] = b0[i] + (m_avg-1.)/m_avg*b1[i] + (m_avg-1.)/m_avg*(m_avg-2.)/m_avg*b2[i]; // EQ. A19
  }

  // Calculation of I1, I2 (and its derivatives)
  static real I1,I1D,I1DD,I1DDD,I1T,I1TT,I2,I2D,I2DD,I2DDD,I2T,I2TT,
    DI1,DI1D,DI1DD,DI1DDD,DI1T,DI1TT,DI2,DI2D,DI2DD,DI2DDD,DI2T,DI2TT;
  I1 = 0.;  // I1
  I1D = 0.; // derivative of I1 wrt rho
  I1DD = 0.; // second derivative of I1 wrt rho
  I1DDD = 0.; // third derivative of I1 wrt rho  
  I1T = 0.; // derivative of I1 wrt T
  I1TT = 0.; // second derivative of I1 wrt T
  I2 = 0.;  // I2
  I2D = 0.; // derivative of I2 wrt rho
  I2DD = 0.; // second derivative of I2 wrt rho
  I2DDD = 0.; // third derivative of I2 wrt rho   
  I2T = 0.; // derivative of I2 wrt T
  I2TT = 0.; // second derivative of I2 wrt T  
  DI1 = 0.; // DetaI1_Deta
  DI1D = 0.; // derivative of DI1 wrt rho
  DI1DD = 0.; // second derivative of DI1D wrt rho
  DI1DDD = 0.; // third derivative of DI1D wrt rho  
  DI1T = 0.; // derivative of DI1 wrt T
  DI1TT = 0.; // second derivative of DI1 wrt T
  DI2 = 0.; // DetaI2_Deta
  DI2D = 0.; // derivative of DI2 wrt rho
  DI2DD = 0.; // second derivative of DI2D wrt rho
  DI2DDD = 0.; // third derivative of DI2D wrt rho   
  DI2T = 0.; // derivative of DI2 wrt T 
  DI2TT = 0.; // second derivative of DI2 wrt T
  real p_eta = 1.;
  auto p_eta1 = 1./eta;
  auto p_eta2 = 1./eta/eta;
  auto p_eta3 = 1./eta/eta/eta;
  for (int i = 0; i < 7; i++) {
    aux = a[i]*p_eta;
    I1 += aux; // EQ. A16
    I1D += i*aux/eta*etaD;
    I1DD += i*a[i]*etaD*(-(p_eta2*etaD) + i*p_eta2*etaD) + i*a[i]*p_eta1*etaDD;
    I1DDD += 2*i*a[i]*(-(p_eta2*etaD) + i*p_eta2*etaD)*etaDD + i*a[i]*etaD*(-(etaD*(-2*p_eta3*etaD + i*p_eta3*etaD)) + 
      i*etaD*(-2*p_eta3*etaD + i*p_eta3*etaD) - p_eta2*etaDD + i*p_eta2*etaDD) + i*a[i]*p_eta1*etaDDD;
    I1T += i*aux/eta*etaT;
    I1TT += i*a[i]*etaT*(-(p_eta2*etaT) + i*p_eta2*etaT) + i*a[i]*p_eta1*etaTT;
    aux = b[i]*p_eta;
    I2 += aux; // EQ. A17
    I2D += i*aux/eta*etaD;
    I2DD += i*b[i]*etaD*(-(p_eta2*etaD) + i*p_eta2*etaD) + i*b[i]*p_eta1*etaDD;
    I2DDD += 2*i*b[i]*(-(p_eta2*etaD) + i*p_eta2*etaD)*etaDD + i*b[i]*etaD*(-(etaD*(-2*p_eta3*etaD + i*p_eta3*etaD)) + 
      i*etaD*(-2*p_eta3*etaD + i*p_eta3*etaD) - p_eta2*etaDD + i*p_eta2*etaDD) + i*b[i]*p_eta1*etaDDD;    
    I2T += i*aux/eta*etaT;
    I2TT += i*b[i]*etaT*(-(p_eta2*etaT) + i*p_eta2*etaT) + i*b[i]*p_eta1*etaTT;
    aux =  a[i]*(i+1.)*p_eta;   
    DI1 += aux; // EQ. A29
    DI1D += i*aux/eta*etaD;
    DI1DD += i*(1 + i)*a[i]*etaD*(-(p_eta2*etaD) + i*p_eta2*etaD) + i*(1 + i)*a[i]*p_eta1*etaDD; 
    DI1DDD += 2*i*(1 + i)*a[i]*(-(p_eta2*etaD) + i*p_eta2*etaD)*etaDD + i*(1 + i)*a[i]*etaD*(-(etaD*(-2*p_eta3*etaD + 
      i*p_eta3*etaD)) + i*etaD*(-2*p_eta3*etaD + i*p_eta3*etaD) - p_eta2*etaDD + i*p_eta2*etaDD) + i*(1 + i)*a[i]*p_eta1*etaDDD;
    DI1T += i*aux/eta*etaT;
    DI1TT += i*(1 + i)*a[i]*etaT*(-(p_eta2*etaT) + i*p_eta2*etaT) + i*(1 + i)*a[i]*p_eta1*etaTT; 
    aux = b[i]*(i+1.)*p_eta;
    DI2 += aux; // EQ. A30
    DI2D += i*aux/eta*etaD;
    DI2DD += i*(1 + i)*b[i]*etaD*(-(p_eta2*etaD) + i*p_eta2*etaD) + i*(1 + i)*b[i]*p_eta1*etaDD;
    DI2DDD += 2*i*(1 + i)*b[i]*(-(p_eta2*etaD) + i*p_eta2*etaD)*etaDD + i*(1 + i)*b[i]*etaD*(-(etaD*(-2*p_eta3*etaD + 
      i*p_eta3*etaD)) + i*etaD*(-2*p_eta3*etaD + i*p_eta3*etaD) - p_eta2*etaDD + i*p_eta2*etaDD) + i*(1 + i)*b[i]*p_eta1*etaDDD;
    DI2T += i*aux/eta*etaT;
    DI2TT += i*(1 + i)*b[i]*etaT*(-(p_eta2*etaT) + i*p_eta2*etaT) + i*(1 + i)*b[i]*p_eta1*etaTT;
    p_eta*=eta;
    p_eta1*=eta;
    p_eta2*=eta;
    p_eta3*=eta;
  }

  // Calculation of C1, C2 (and its derivatives)
  static real C1,C1D,C1DD,C1DDD,C1T,C1TT,C2,C2D,C2DD,C2DDD,C2T,C2TT;
  C1 = 1./(1. + m_avg*(8.*eta-2.*eta*eta)/((1.-eta)*(1.-eta)*(1.-eta)*(1.-eta)) + 
      (1.-m_avg)*(20.*eta-27.*eta*eta+12.*(eta*eta*eta)-2.*(eta*eta*eta*eta))/((1.-eta)*(2.-eta)*(1.-eta)*(2.-eta))); // EQ. A11
  
  C1D =  (-2*(-2 + eta)*((eta-1.)*(eta-1.)*(eta-1.))*(((eta-1.)*(eta-1.))*(20 - 24*eta + 6*(eta*eta) + 
      (eta*eta*eta)) + m_avg*(12 + 96*eta - 186*(eta*eta) + 115*(eta*eta*eta) - 26*(eta*eta*eta*eta) + 
      (eta*eta*eta*eta*eta)))*etaD)/((((eta-1.)*(eta-1.))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + 
      (eta*eta*eta*eta)) - m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 
      2*(eta*eta*eta*eta*eta)))*(((eta-1.)*(eta-1.))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + 
      (eta*eta*eta*eta)) - m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 
      2*(eta*eta*eta*eta*eta))));  // derivative of C1 wrt rho

  C1DD =  (-2*((eta-1.)*(eta-1.))*(-((-(((1.-eta)*(1.-eta)*(1.-eta)*(1.-eta))*(-1072 + 3936*eta - 6456*(eta*eta) + 6200*(eta*eta*eta) - 3744*(eta*eta*eta*eta) + 
      1368*(eta*eta*eta*eta*eta) - 250*(eta*eta*eta*eta*eta*eta) + 6*(eta*eta*eta*eta*eta*eta*eta) + 3*(eta*eta*eta*eta*eta*eta*eta*eta))) + m_avg*((eta-1.)*(eta-1.))*(528 + 8208*eta - 39252*(eta*eta) + 
      77200*(eta*eta*eta) - 88460*(eta*eta*eta*eta) + 65596*(eta*eta*eta*eta*eta) - 32137*(eta*eta*eta*eta*eta*eta) + 10004*(eta*eta*eta*eta*eta*eta*eta) - 1756*(eta*eta*eta*eta*eta*eta*eta*eta) + 120*(eta*eta*eta*eta*eta*eta*eta*eta*eta) + 3*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)) + (m_avg*m_avg)*(576 + 5040*eta + 8172*(eta*eta) - 95280*(eta*eta*eta) + 233112*(eta*eta*eta*eta) - 311268*(eta*eta*eta*eta*eta) + 269557*(eta*eta*eta*eta*eta*eta) - 160048*(eta*eta*eta*eta*eta*eta*eta) + 65442*(eta*eta*eta*eta*eta*eta*eta*eta) - 17812*(eta*eta*eta*eta*eta*eta*eta*eta*eta) + 2971*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 252*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 6*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)))*(etaD*etaD)) + 
      (2 - 3*eta + (eta*eta))*(-(((1.-eta)*(1.-eta)*(1.-eta)*(1.-eta))*(-80 - 64*eta + 448*(eta*eta) - 508*(eta*eta*eta) + 240*(eta*eta*eta*eta) - 
      46*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta*eta))) + m_avg*((eta-1.)*(eta-1.))*(48 + 720*eta + 108*(eta*eta) - 4276*(eta*eta*eta) + 6858*(eta*eta*eta*eta) - 4963*(eta*eta*eta*eta*eta) + 1908*(eta*eta*eta*eta*eta*eta) - 378*(eta*eta*eta*eta*eta*eta*eta) + 28*(eta*eta*eta*eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta*eta*eta*eta)) + (m_avg*m_avg)*eta*(144 + 1476*eta - 480*(eta*eta) - 9750*(eta*eta*eta) + 
      20517*(eta*eta*eta*eta) - 19738*(eta*eta*eta*eta*eta) + 10880*(eta*eta*eta*eta*eta*eta) - 3608*(eta*eta*eta*eta*eta*eta*eta) + 697*(eta*eta*eta*eta*eta*eta*eta*eta) - 68*(eta*eta*eta*eta*eta*eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)))*etaDD))/
      ((-(((eta-1.)*(eta-1.))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta))) + m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta)))*(-(((eta-1.)*(eta-1.))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta))) + m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta)))*(-(((eta-1.)*(eta-1.))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta))) + m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta))));  // second derivative of C1 wrt rho

  C1DDD = (-2*(-1 + eta)*(12*(-2 + eta)*(((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta))*
    (1728 - 5936*eta + 9376*(eta*eta) - 8936*(eta*eta*eta) + 5752*(eta*eta*eta*eta) - 2600*(eta*eta*eta*eta*eta) + 784*(eta*eta*eta*eta*eta*eta) - 128*(eta*eta*eta*eta*eta*eta*eta) + 4*(eta*eta*eta*eta*eta*eta*eta*eta) + 
    (eta*eta*eta*eta*eta*eta*eta*eta*eta)) - m_avg*((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta))*(-592 - 20848*eta + 86904*(eta*eta) - 161556*(eta*eta*eta) + 181188*(eta*eta*eta*eta) - 
    137532*(eta*eta*eta*eta*eta) + 74378*(eta*eta*eta*eta*eta*eta) - 28687*(eta*eta*eta*eta*eta*eta*eta) + 7428*(eta*eta*eta*eta*eta*eta*eta*eta) - 1102*(eta*eta*eta*eta*eta*eta*eta*eta*eta) + 56*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 3*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)) + 
    (m_avg*m_avg)*((-1 + eta)*(-1 + eta))*(144 + 6864*eta + 50652*(eta*eta) - 527612*(eta*eta*eta) + 1780992*(eta*eta*eta*eta) - 3404740*(eta*eta*eta*eta*eta) + 4322670*(eta*eta*eta*eta*eta*eta) - 3928155*(eta*eta*eta*eta*eta*eta*eta) + 
    2658426*(eta*eta*eta*eta*eta*eta*eta*eta) - 1363355*(eta*eta*eta*eta*eta*eta*eta*eta*eta) + 528174*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 150445*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 29674*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 3597*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 200*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)) + 
    2*(m_avg*m_avg*m_avg)*(432 + 4536*eta + 9630*(eta*eta) - 28836*(eta*eta*eta) - 279732*(eta*eta*eta*eta) + 1417776*(eta*eta*eta*eta*eta) - 3227854*(eta*eta*eta*eta*eta*eta) + 4647053*(eta*eta*eta*eta*eta*eta*eta) - 4730052*(eta*eta*eta*eta*eta*eta*eta*eta) + 3588948*(eta*eta*eta*eta*eta*eta*eta*eta*eta) - 
    2084234*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 934340*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 320758*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 82110*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 14884*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 1739*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 108*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)))*(etaD*etaD*etaD) - 
    3*(-1 + eta)*(((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta))*(4288 - 7168*eta - 20672*(eta*eta) + 88384*(eta*eta*eta) - 149696*(eta*eta*eta*eta) + 
    153952*(eta*eta*eta*eta*eta) - 106016*(eta*eta*eta*eta*eta*eta) + 49792*(eta*eta*eta*eta*eta*eta*eta) - 15512*(eta*eta*eta*eta*eta*eta*eta*eta) + 2928*(eta*eta*eta*eta*eta*eta*eta*eta*eta) - 244*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 12*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 3*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)) - 
    3*m_avg*((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta))*(-704 - 16640*eta + 39008*(eta*eta) + 73600*(eta*eta*eta) - 430704*(eta*eta*eta*eta) + 854032*(eta*eta*eta*eta*eta) - 
    1016720*(eta*eta*eta*eta*eta*eta) + 822256*(eta*eta*eta*eta*eta*eta*eta) - 470152*(eta*eta*eta*eta*eta*eta*eta*eta) + 190604*(eta*eta*eta*eta*eta*eta*eta*eta*eta) - 53394*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 9676*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 959*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 
    22*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 3*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)) + 3*(m_avg*m_avg)*((-1 + eta)*(-1 + eta))*(768 + 10368*eta + 59232*(eta*eta) - 223072*(eta*eta*eta) - 198524*(eta*eta*eta*eta) + 
    2059472*(eta*eta*eta*eta*eta) - 4797576*(eta*eta*eta*eta*eta*eta) + 6509088*(eta*eta*eta*eta*eta*eta*eta) -6019223*(eta*eta*eta*eta*eta*eta*eta*eta) + 4011784*(eta*eta*eta*eta*eta*eta*eta*eta*eta) - 1962780*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 
    702204*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 179371*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 31076*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 3282*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 160*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)) + (m_avg*m_avg*m_avg)*eta*(6912 + 76032*eta + 193824*(eta*eta) - 1246140*(eta*eta*eta) - 
    99432*(eta*eta*eta*eta) + 9565692*(eta*eta*eta*eta*eta) - 26467344*(eta*eta*eta*eta*eta*eta) + 40575759*(eta*eta*eta*eta*eta*eta*eta) - 42200002*(eta*eta*eta*eta*eta*eta*eta*eta) + 31950469*(eta*eta*eta*eta*eta*eta*eta*eta*eta) - 
    18124108*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 7761457*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 2490282*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 585199*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 96432*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 10280*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 600*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 12*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)))*etaD*etaDD + 
    (-2 + eta)*((-1 + eta)*(-1 + eta))*(((-1 + eta)*(-1 + eta))*(20 - 24*eta + 6*(eta*eta) + (eta*eta*eta)) +m_avg*(12 + 96*eta - 186*(eta*eta) + 115*(eta*eta*eta) - 26*(eta*eta*eta*eta) + 
    (eta*eta*eta*eta*eta)))*(((-1 + eta)*(-1 + eta))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta)) - m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 
    2*(eta*eta*eta*eta*eta))*((-1 + eta)*(-1 + eta))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta)) - m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 
    2*(eta*eta*eta*eta*eta)))*etaDDD))/((((-1 + eta)*(-1 + eta))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta)) -m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 
    2*(eta*eta*eta*eta*eta)))*(((-1 + eta)*(-1 + eta))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta)) - m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 
    2*(eta*eta*eta*eta*eta)))*(((-1 + eta)*(-1 + eta))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta)) - m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 
    2*(eta*eta*eta*eta*eta)))*(((-1 + eta)*(-1 + eta))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta)) - m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta)))); // third derivative of C1 wrt rho

  C1T =  (-2*(-2 + eta)*((eta-1.)*(eta-1.)*(eta-1.))*(((eta-1.)*(eta-1.))*(20 - 24*eta + 6*(eta*eta) + (eta*eta*eta)) + 
      m_avg*(12 + 96*eta - 186*(eta*eta) + 115*(eta*eta*eta) - 26*(eta*eta*eta*eta) + (eta*eta*eta*eta*eta)))*
      etaT)/((((eta-1.)*(eta-1.))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta)) - 
      m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta)))*(((eta-1.)*(eta-1.))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta)) - 
      m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta))));  // derivative of C1 wrt T
  
  C1TT = (-2*((eta-1.)*(eta-1.))*(-((-(((1.-eta)*(1.-eta)*(1.-eta)*(1.-eta))*(-1072 + 3936*eta - 6456*(eta*eta) + 6200*(eta*eta*eta) - 3744*(eta*eta*eta*eta) + 
      1368*(eta*eta*eta*eta*eta) - 250*(eta*eta*eta*eta*eta*eta) + 6*(eta*eta*eta*eta*eta*eta*eta) + 3*(eta*eta*eta*eta*eta*eta*eta*eta))) + m_avg*((eta-1.)*(eta-1.))*(528 + 8208*eta - 39252*(eta*eta) + 
      77200*(eta*eta*eta) - 88460*(eta*eta*eta*eta) + 65596*(eta*eta*eta*eta*eta) - 32137*(eta*eta*eta*eta*eta*eta) + 10004*(eta*eta*eta*eta*eta*eta*eta) - 1756*(eta*eta*eta*eta*eta*eta*eta*eta) + 120*(eta*eta*eta*eta*eta*eta*eta*eta*eta) + 3*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)) + (m_avg*m_avg)*(576 + 5040*eta + 8172*(eta*eta) - 95280*(eta*eta*eta) + 233112*(eta*eta*eta*eta) - 311268*(eta*eta*eta*eta*eta) + 269557*(eta*eta*eta*eta*eta*eta) - 160048*(eta*eta*eta*eta*eta*eta*eta) + 65442*(eta*eta*eta*eta*eta*eta*eta*eta) - 17812*(eta*eta*eta*eta*eta*eta*eta*eta*eta) + 2971*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) - 252*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta) + 6*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)))*(etaT*etaT)) + 
      (2 - 3*eta + (eta*eta))*(-(((1.-eta)*(1.-eta)*(1.-eta)*(1.-eta))*(-80 - 64*eta + 448*(eta*eta) - 508*(eta*eta*eta) + 240*(eta*eta*eta*eta) - 
      46*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta*eta))) + m_avg*((eta-1.)*(eta-1.))*(48 + 720*eta + 108*(eta*eta) - 4276*(eta*eta*eta) + 6858*(eta*eta*eta*eta) - 4963*(eta*eta*eta*eta*eta) + 1908*(eta*eta*eta*eta*eta*eta) - 378*(eta*eta*eta*eta*eta*eta*eta) + 28*(eta*eta*eta*eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta*eta*eta*eta)) + (m_avg*m_avg)*eta*(144 + 1476*eta - 480*(eta*eta) - 9750*(eta*eta*eta) + 
      20517*(eta*eta*eta*eta) - 19738*(eta*eta*eta*eta*eta) + 10880*(eta*eta*eta*eta*eta*eta) - 3608*(eta*eta*eta*eta*eta*eta*eta) + 697*(eta*eta*eta*eta*eta*eta*eta*eta) - 68*(eta*eta*eta*eta*eta*eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta*eta*eta*eta*eta*eta)))*etaTT))/
      ((-(((eta-1.)*(eta-1.))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta))) + m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta)))*(-(((eta-1.)*(eta-1.))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta))) + m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta)))*(-(((eta-1.)*(eta-1.))*(-4 - 8*eta + 14*(eta*eta) - 6*(eta*eta*eta) + (eta*eta*eta*eta))) + m_avg*eta*(12 + 27*eta - 70*(eta*eta) + 51*(eta*eta*eta) - 16*(eta*eta*eta*eta) + 2*(eta*eta*eta*eta*eta)))); // second derivative of C1 wrt T

  C2 = -1.*C1*C1*(m_avg*(-4.*eta*eta+20.*eta+8.)/((1.-eta)*(1.-eta)*(1.-eta)*(1.-eta)*(1.-eta)) + 
      (1.-m_avg)*(2.*(eta*eta*eta)+12.*eta*eta-48.*eta+40.)/((1.-eta)*(2.-eta)*(1.-eta)*(2.-eta)*(1.-eta)*(2.-eta))); // EQ. A31

  C2D = (2*C1*(-2*(2 - 3*eta + (eta*eta))*(((eta-1.)*(eta-1.))*(20 - 24*eta + 6*(eta*eta) + (eta*eta*eta)) + 
      m_avg*(12 + 96*eta - 186*(eta*eta) + 115*(eta*eta*eta) - 26*(eta*eta*eta*eta) + (eta*eta*eta*eta*eta)))*C1D + 
      3*C1*(((eta-1.)*(eta-1.))*(-44 + 80*eta - 48*(eta*eta) + 8*(eta*eta*eta) + (eta*eta*eta*eta)) + m_avg*(-116 - 40*eta + 428*(eta*eta) - 456*(eta*eta*eta) + 
      197*(eta*eta*eta*eta) - 34*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta)))*etaD))/(((-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta))*((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta))); // derivative of C2 wrt rho

  C2DD = (-4*((2 - 3*eta + (eta*eta))*(2 - 3*eta + (eta*eta)))*(((eta-1.)*(eta-1.))*(20 - 24*eta + 6*(eta*eta) + (eta*eta*eta)) + 
      m_avg*(12 + 96*eta - 186*(eta*eta) + 115*(eta*eta*eta) - 26*(eta*eta*eta*eta) + (eta*eta*eta*eta*eta)))*(C1D*C1D) + 
      24*C1*(2 - 3*eta + (eta*eta))*(((eta-1.)*(eta-1.))*(-44 + 80*eta - 48*(eta*eta) + 8*(eta*eta*eta) + (eta*eta*eta*eta)) + m_avg*(-116 - 40*eta + 428*(eta*eta) - 
      456*(eta*eta*eta) + 197*(eta*eta*eta*eta) - 34*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta)))*C1D*etaD - 2*C1*(2*((2 - 3*eta + (eta*eta))*(2 - 3*eta + (eta*eta)))*(((eta-1.)*(eta-1.))*(20 - 24*eta + 6*(eta*eta) + (eta*eta*eta)) + 
      m_avg*(12 + 96*eta - 186*(eta*eta) + 115*(eta*eta*eta) - 26*(eta*eta*eta*eta) + (eta*eta*eta*eta*eta)))*C1DD + 3*C1*(4*(((eta-1.)*(eta-1.))*(92 - 220*eta + 200*(eta*eta) - 
      80*(eta*eta*eta) + 10*(eta*eta*eta*eta) + (eta*eta*eta*eta*eta)) + m_avg*(484 - 588*eta - 476*(eta*eta) + 1260*(eta*eta*eta) - 910*(eta*eta*eta*eta) + 301*(eta*eta*eta*eta*eta) - 
      42*(eta*eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta*eta)))*(etaD*etaD) - (2 - 3*eta + (eta*eta))*(((eta-1.)*(eta-1.))*(-44 + 80*eta - 48*(eta*eta) + 8*(eta*eta*eta) + (eta*eta*eta*eta)) + 
      m_avg*(-116 - 40*eta + 428*(eta*eta) - 456*(eta*eta*eta) + 197*(eta*eta*eta*eta) - 34*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta)))*etaDD)))/(((-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta))*((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta))); // second derivative of C2 wrt rho

  C2DDD = (36*(((-1 + eta)*(-1 + eta))*(-44 + 80*eta - 48*(eta*eta) + 8*(eta*eta*eta) + (eta*eta*eta*eta)) + m_avg*(-116 - 40*eta + 428*(eta*eta) - 456*(eta*eta*eta) + 197*(eta*eta*eta*eta) - 
    34*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta)))*etaD*((C1D*C1D) + C1*C1DD))/(((-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta))*((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta))) + (36*C1*C1D*(-4*(((-1 + eta)*(-1 + eta))*
    (92 - 220*eta + 200*(eta*eta) - 80*(eta*eta*eta) + 10*(eta*eta*eta*eta) + (eta*eta*eta*eta*eta))+ m_avg*(484 - 588*eta - 476*(eta*eta) + 1260*(eta*eta*eta) - 910*(eta*eta*eta*eta) + 301*(eta*eta*eta*eta*eta) - 42*(eta*eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta*eta)))*(etaD*etaD) + 
    (2 - 3*eta + (eta*eta))*(((-1 + eta)*(-1 + eta))*(-44 + 80*eta - 48*(eta*eta) + 8*(eta*eta*eta) + (eta*eta*eta*eta)) + m_avg*(-116 - 40*eta + 428*(eta*eta) - 456*(eta*eta*eta) + 197*(eta*eta*eta*eta) - 
    34*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta)))*etaDD))/(((-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta))*((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta))) - 2*((4*m_avg*(-2 - 5*eta + (eta*eta)))/((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)) - 
    (2*(-1 + m_avg)*(20 - 24*eta + 6*(eta*eta) + (eta*eta*eta)))/(((-2 + eta)*(-2 + eta)*(-2 + eta))*((-1 + eta)*(-1 + eta)*(-1 + eta))))*(3*C1D*C1DD + C1*C1DDD) + 
    (6*(C1*C1)*(20*(((-1 + eta)*(-1 + eta))*(-188 + 552*eta - 660*(eta*eta) + 400*(eta*eta*eta) - 120*(eta*eta*eta*eta) + 12*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta)) + 
    m_avg*(-1604 + 3424*eta - 1568*(eta*eta) - 2016*(eta*eta*eta) + 2940*(eta*eta*eta*eta) - 1596*(eta*eta*eta*eta*eta) + 427*(eta*eta*eta*eta*eta*eta) - 50*(eta*eta*eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta*eta*eta)))*(etaD*etaD*etaD) - 
    12*(2 - 3*eta + (eta*eta))*(((-1 + eta)*(-1 + eta))*(92 - 220*eta + 200*(eta*eta) - 80*(eta*eta*eta) + 10*(eta*eta*eta*eta) + (eta*eta*eta*eta*eta))+ m_avg*(484 - 588*eta - 476*(eta*eta) + 1260*(eta*eta*eta) - 910*(eta*eta*eta*eta) + 
    301*(eta*eta*eta*eta*eta) - 42*(eta*eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta*eta)))*etaD*etaDD + ((2 - 3*eta + (eta*eta))*(2 - 3*eta + (eta*eta)))*
    (((-1 + eta)*(-1 + eta))*(-44 + 80*eta - 48*(eta*eta) + 8*(eta*eta*eta) + (eta*eta*eta*eta)) + m_avg*(-116 - 40*eta + 428*(eta*eta) - 456*(eta*eta*eta) + 197*(eta*eta*eta*eta) - 
    34*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta)))*etaDDD))/(((-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta))*((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta))); // third derivative of C2 wrt rho
  
  C2T = (2*C1*(-2*(2 - 3*eta + (eta*eta))*(((eta-1.)*(eta-1.))*
      (20 - 24*eta + 6*(eta*eta) + (eta*eta*eta)) + m_avg*(12 + 96*eta - 186*(eta*eta) + 115*(eta*eta*eta) - 26*(eta*eta*eta*eta) + (eta*eta*eta*eta*eta)))*
      C1T + 3*C1*(((eta-1.)*(eta-1.))*(-44 + 80*eta - 48*(eta*eta) + 8*(eta*eta*eta) + (eta*eta*eta*eta)) + 
      m_avg*(-116 - 40*eta + 428*(eta*eta) - 456*(eta*eta*eta) + 197*(eta*eta*eta*eta) - 
      34*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta)))*etaT))/(((-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta))*((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta))); // derivative of C2 wrt T   

  C2TT = (-4*((2 - 3*eta + (eta*eta))*(2 - 3*eta + (eta*eta)))*(((eta-1.)*(eta-1.))*(20 - 24*eta + 6*(eta*eta) + (eta*eta*eta)) + 
      m_avg*(12 + 96*eta - 186*(eta*eta) + 115*(eta*eta*eta) - 26*(eta*eta*eta*eta) + (eta*eta*eta*eta*eta)))*(C1T*C1T) + 
      24*C1*(2 - 3*eta + (eta*eta))*(((eta-1.)*(eta-1.))*(-44 + 80*eta - 48*(eta*eta) + 8*(eta*eta*eta) + (eta*eta*eta*eta)) + m_avg*(-116 - 40*eta + 428*(eta*eta) - 
      456*(eta*eta*eta) + 197*(eta*eta*eta*eta) - 34*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta)))*C1T*etaT - 2*C1*(2*((2 - 3*eta + (eta*eta))*(2 - 3*eta + (eta*eta)))*(((eta-1.)*(eta-1.))*(20 - 24*eta + 6*(eta*eta) + (eta*eta*eta)) + 
      m_avg*(12 + 96*eta - 186*(eta*eta) + 115*(eta*eta*eta) - 26*(eta*eta*eta*eta) + (eta*eta*eta*eta*eta)))*C1TT + 3*C1*(4*(((eta-1.)*(eta-1.))*(92 - 220*eta + 200*(eta*eta) - 
      80*(eta*eta*eta) + 10*(eta*eta*eta*eta) + (eta*eta*eta*eta*eta)) + m_avg*(484 - 588*eta - 476*(eta*eta) + 1260*(eta*eta*eta) - 910*(eta*eta*eta*eta) + 301*(eta*eta*eta*eta*eta) - 
      42*(eta*eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta*eta)))*(etaT*etaT) - (2 - 3*eta + (eta*eta))*(((eta-1.)*(eta-1.))*(-44 + 80*eta - 48*(eta*eta) + 8*(eta*eta*eta) + (eta*eta*eta*eta)) + 
      m_avg*(-116 - 40*eta + 428*(eta*eta) - 456*(eta*eta*eta) + 197*(eta*eta*eta*eta) - 34*(eta*eta*eta*eta*eta) + (eta*eta*eta*eta*eta*eta)))*etaTT)))/(((-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta)*(-2 + eta))*((-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta)*(-1 + eta))); // second derivative of C2 wrt T

  summation = 0.;
  real summationD = 0.,summationDD = 0.,summationDDD = 0.;
  summationT = 0.;
  summationTT = 0.;
  for (int i = 0; i < n_comp; i++){
    aux = (ghs[i*n_comp+i]*ghs[i*n_comp+i]);
    summation += x[i]*(m(i)-1.)/ghs[i*n_comp+i]*Dghs[i*n_comp+i];
    summationD += ((-1 + m(i))*x[i]*DghsD[i*n_comp+i])/ghs[i*n_comp+i] - ((-1 + m(i))*x[i]*Dghs[i*n_comp+i]*ghsD[i*n_comp+i])/aux; 
    summationDD += (-2*(-1 + m(i))*x[i]*DghsD[i*n_comp+i]*ghsD[i*n_comp+i])/(aux) + (2*(-1 + m(i))*x[i]*Dghs[i*n_comp+i]*
    (ghsD[i*n_comp+i]*ghsD[i*n_comp+i]))/(ghs[i*n_comp+i]*ghs[i*n_comp+i]*ghs[i*n_comp+i]) + ((-1 + m(i))*x[i]*DghsDD[i*n_comp+i])/ghs[i*n_comp+i] - 
    ((-1 + m(i))*x[i]*Dghs[i*n_comp+i]*ghsDD[i*n_comp+i])/(aux);
    summationDDD += ((-1 + m(i))*x[i]*(((ghs[i*n_comp+i])*(ghs[i*n_comp+i])*(ghs[i*n_comp+i]))*DghsDDD[i*n_comp+i] - 
       6*Dghs[i*n_comp+i]*((ghsD[i*n_comp+i])*(ghsD[i*n_comp+i]))*ghsD[i*n_comp+i] + 
       2*ghs[i*n_comp+i]*(((ghsD[i*n_comp+i])*(ghsD[i*n_comp+i]))*DghsD[i*n_comp+i] + Dghs[i*n_comp+i]*ghsDD[i*n_comp+i]*ghsD[i*n_comp+i] + 
       2*ghsD[i*n_comp+i]*(DghsD[i*n_comp+i]*ghsD[i*n_comp+i] + Dghs[i*n_comp+i]*ghsDD[i*n_comp+i])) - 
       ((ghs[i*n_comp+i])*(ghs[i*n_comp+i]))*(ghsDD[i*n_comp+i]*DghsD[i*n_comp+i] + 2*ghsD[i*n_comp+i]*DghsDD[i*n_comp+i] + 
       DghsDD[i*n_comp+i]*ghsD[i*n_comp+i] + 2*DghsD[i*n_comp+i]*ghsDD[i*n_comp+i] + 
       Dghs[i*n_comp+i]*ghsDDD[i*n_comp+i])))/((ghs[i*n_comp+i])*(ghs[i*n_comp+i])*(ghs[i*n_comp+i])*(ghs[i*n_comp+i]));
    summationT += ((-1 + m(i))*x[i]*DghsT[i*n_comp+i])/ghs[i*n_comp+i] - ((-1 + m(i))*x[i]*Dghs[i*n_comp+i]*ghsT[i*n_comp+i])/aux;
    summationTT += (-2*(-1 + m(i))*x[i]*DghsT[i*n_comp+i]*ghsT[i*n_comp+i])/(aux) + (2*(-1 + m(i))*x[i]*Dghs[i*n_comp+i]*
    ((ghsT[i*n_comp+i])*(ghsT[i*n_comp+i])))/(ghs[i*n_comp+i]*ghs[i*n_comp+i]*ghs[i*n_comp+i]) + ((-1 + m(i))*x[i]*DghsTT[i*n_comp+i])/ghs[i*n_comp+i] - 
    ((-1 + m(i))*x[i]*Dghs[i*n_comp+i]*ghsTT[i*n_comp+i])/(aux);        
  }

  static real Zid,ZidD,ZidDD,ZidDDD,ZidT,ZidTT;
  Zid = 1.; 
  ZidD = 0.;  // derivative of Zid wrt rho
  ZidDD = 0.;  // second derivative of Zid wrt rho
  ZidDDD = 0.;  // third derivative of Zid wrt rho
  ZidT = 0.; //derivative of Zid wrt T
  ZidTT = 0.; // second derivative of Zid wrt T

  static real Zhc,ZhcD,ZhcDD,ZhcDDD,ZhcT,ZhcTT;
  Zhc = m_avg*Zhs - summation; // EQ. A25
  ZhcD = m_avg*ZhsD - summationD; // derivative of Zhc wrt rho
  ZhcDD = m_avg*ZhsDD - summationDD; // second derivative of Zhc wrt rho
  ZhcDDD = m_avg*ZhsDDD - summationDDD; // third derivative of Zhc wrt rho
  ZhcT = m_avg*ZhsT - summationT; // derivative of Zhc wrt T
  ZhcTT = m_avg*ZhsTT - summationTT; // second derivative of Zhc wrt T  

  static real Zdisp,ZdispD,ZdispDD,ZdispDDD,ZdispT,ZdispTT;
  Zdisp = -2.*pi*num_den*DI1*m2epsilonsigma3 - 
    pi*num_den*m_avg*(C1*DI2 + C2*eta*I2)*m2epsilon2sigma3; // EQ. A28

  ZdispD = pi*(-2*m2epsilonsigma3*num_den*DI1D - m2epsilon2sigma3*m_avg*num_den*(DI2*C1D + C1*DI2D + C2*eta*I2D + 
    I2*(eta*C2D + C2*etaD)) - 2*m2epsilonsigma3*DI1*num_denD - m2epsilon2sigma3*m_avg*(C1*DI2 + C2*I2*eta)*num_denD);  // derivative of Zdisp wrt rho
 
  ZdispDD = -2*m2epsilonsigma3*pi*(2*DI1D*num_denD + num_den*DI1DD + DI1*num_denDD) - 
    m2epsilon2sigma3*pi*m_avg*(2*(DI2*C1D + C1*DI2D + C2*eta*I2D + I2*(eta*C2D + C2*etaD))*num_denD + 
    num_den*(2*C1D*DI2D + 2*(I2*C2D + C2*I2D)*etaD + DI2*C1DD + C1*DI2DD + eta*(2*C2D*I2D + I2*C2DD + C2*I2DD) + 
    C2*I2*etaDD) + (C1*DI2 + C2*I2*eta)*num_denDD); // second derivative of Zdisp wrt rho

  ZdispDDD = -2*m2epsilonsigma3*pi*(3*num_denD*DI1DD + 3*DI1D*num_denDD + num_den*DI1DDD + 
    DI1*num_denDDD) - m2epsilon2sigma3*pi*m_avg*(3*num_denD*(2*C1D*DI2D + 2*(I2*C2D + C2*I2D)*etaD + 
    DI2*C1DD + C1*DI2DD + eta*(2*C2D*I2D + I2*C2DD + C2*I2DD) + C2*I2*etaDD) + 3*(DI2*C1D + C1*DI2D + 
    C2*eta*I2D + I2*(eta*C2D + C2*etaD))*num_denDD + num_den*(3*DI2D*C1DD + 3*C1D*DI2DD + 3*etaD*(2*C2D*I2D + 
    I2*C2DD + C2*I2DD) + 3*(I2*C2D + C2*I2D)*etaDD + DI2*C1DDD + C1*DI2DDD + eta*(3*I2D*C2DD + 
    3*C2D*I2DD + I2*C2DDD + C2*I2DDD) + C2*I2*etaDDD) + (C1*DI2 + C2*I2*eta)*num_denDDD);

  ZdispT = pi*num_den*(-2*m2epsilonsigma3*DI1T - m_avg*(C1*DI2 + C2*I2*eta)*m2epsilon2sigma3T - 2*DI1*m2epsilonsigma3T - 
    m2epsilon2sigma3*m_avg*(DI2*C1T + C1*DI2T + C2*eta*I2T + I2*(eta*C2T + C2*etaT)));  // derivative of Zdisp wrt T

  ZdispTT = pi*num_den*(-2*(2*DI1T*m2epsilonsigma3T + m2epsilonsigma3*DI1TT + DI1*m2epsilonsigma3TT) - m_avg*
    (2*m2epsilon2sigma3T*(DI2*C1T + C1*DI2T + C2*eta*I2T + I2*(eta*C2T + C2*etaT)) + (C1*DI2 + C2*I2*eta)*m2epsilon2sigma3TT + 
    m2epsilon2sigma3*(2*C1T*DI2T + 2*(I2*C2T + C2*I2T)*etaT + DI2*C1TT + C1*DI2TT + eta*(2*C2T*I2T + I2*C2TT + 
    C2*I2TT) + C2*I2*etaTT)));  //second derivative of Zdisp wrt T  


  // Association term
  // Derivatives will be calculated wrt numerical density, they will be transformed to derivatives wrt rho (mol/m³) at the final step of the derivative calculation
  static real Zassoc,ZassocD,ZassocDD,ZassocDDD,ZassocT,ZassocTT,a_assoc,a_assocT,a_assocTT;
  Zassoc = 0;
  ZassocD = 0;
  ZassocDD = 0;
  ZassocDDD = 0;
  ZassocT = 0;
  ZassocTT = 0;
  a_assoc = 0.;
  a_assocT = 0.;
  a_assocTT = 0.;

  // Testing if there is any associative species
  bool test_assoc = false;
  for (auto i=0; i<n_comp; i++) {if (types[pcsaft_type[i]]==A_PCSAFT) test_assoc = true; }

  real sum1,sum1T,sum1TT,sum2,sum2T,sum2TT,sum3;
  if (test_assoc) 
  { 
    int order = n_comp*s;

    // Creating S matrix 
    for (int i = 0; i < s; i++)
    {
      for (int j = 0; j < n_comp; j++){
        S(i,j) = n_donor_sites[j];
        S(i,j) = n_acceptor_sites[j];
      } 
    }

    // Creating delta, deltaD, deltaDD, deltaDDD, deltaDT and deltaTT matrices
    static real eassoc_ki, kappa_ki, sigma_ki;                     
    for (int l = 0; l < s; l++) {
      for (int k = 0; k < n_comp; k++) {
        for (int j = 0; j < s; j++){
          for (int i = 0; i < n_comp; i++){
            if (l==j) continue;
            eassoc_ki = (eassoc(k)+eassoc(i))/2.; // EQ. 2 Association term paper Gross, Sadowski 
            aux = sqrt(sigma(k)*sigma(i))/(0.5*(sigma(k)+sigma(i)));
            kappa_ki = sqrt(kappa(k)*kappa(i))*aux*aux*aux; // EQ. 3 Association term paper Gross, Sadowski 
            sigma_ki = (sigma(k)+sigma(i))/2.;

            // Calculation of auxiliary variables to obtain delta and its derivatives 
            aux = 1./(1.-zeta[3]);
            aux21 = (d[k]*d[i]/(d[k]+d[i]));
            aux22 = 3.*zeta[2]/(1.-zeta[3])/(1.-zeta[3]);
            aux2 = aux21*aux22;
            aux31 = (d[k]*d[i]/(d[k]+d[i]))*(d[k]*d[i]/(d[k]+d[i]));
            aux32 = 2.*zeta[2]*zeta[2]/(1.-zeta[3])/(1.-zeta[3])/(1.-zeta[3]);
            aux3 = aux31*aux32;
            auxD = zetaD[3]/((-1 + zeta[3])*(-1 + zeta[3]));
            auxDD = (-2*(zetaD[3]*zetaD[3]) + (-1 + zeta[3])*zetaDD[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
            auxDDD = (6*((zetaD[3])*(zetaD[3])*(zetaD[3])) - 6*(-1 + zeta[3])*zetaD[3]*
              zetaDD[3] + ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDDD[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
            aux21D = 0.;
            aux21DD = 0.;
            aux21DDD = 0;
            aux22D = (3*(-1 + zeta[3])*zetaD[2] - 6*zeta[2]*zetaD[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
            aux22DD = (3*(-4*(-1 + zeta[3])*zetaD[2]*zetaD[3] + 6*zeta[2]*(zetaD[3]*zetaD[3]) + ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDD[2] - 
              2*zeta[2]*(-1 + zeta[3])*zetaDD[3]))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
            aux22DDD = (-3*(-6*(-1 + zeta[3])*zetaD[2]*(3*((zetaD[3])*(zetaD[3])) - (-1 + zeta[3])*zetaDD[3]) + 
              ((-1 + zeta[3])*(-1 + zeta[3]))*(6*zetaD[3]*zetaDD[2] - (-1 + zeta[3])*zetaDDD[2]) + 
              2*zeta[2]*(12*((zetaD[3])*(zetaD[3])*(zetaD[3])) - 9*(-1 + zeta[3])*zetaD[3]*zetaDD[3] + 
              ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDDD[3])))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
            aux31D = 0.;
            aux31DD = 0.;
            aux31DDD = 0;
            aux32D = (zeta[2]*((4. - 4.*zeta[3])*zetaD[2] + 6.*zeta[2]*zetaD[3]))/((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]));
            aux32DD = (-4.*((1. - zeta[3])*(1. - zeta[3]))*(zetaD[2]*zetaD[2]) + zeta[2]*(-24. + 24.*zeta[3])*zetaD[2]*zetaD[3] + 
              zeta[2]*(-4.*((1. - zeta[3])*(1. - zeta[3]))*zetaDD[2] + zeta[2]*(-24.*(zetaD[3]*zetaD[3]) + (-6. + 6.*zeta[3])*zetaDD[3])))
              /((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));        
            aux32DDD = (2*(18.*((-1. + zeta[3])*(-1. + zeta[3]))*zetaD[3]*((zetaD[2]*zetaD[2]) + zeta[2]*zetaDD[2]) + 
              18.*zeta[2]*(-1. + zeta[3])*zetaD[2]*(-4.*(zetaD[3]*zetaD[3]) + (-1. + zeta[3])*zetaDD[3]) - 
              2.*((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]))*(3.*zetaD[2]*zetaDD[2] + 
              zeta[2]*zetaDDD[2]) + 3*(zeta[2]*zeta[2])*(20*(zetaD[3]*zetaD[3]*zetaD[3]) - 
              12.*(-1. + zeta[3])*zetaD[3]*zetaDD[3] + ((-1. + zeta[3])*(-1. + zeta[3]))*zetaDDD[3])))/((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3]));;
            aux2D = aux21D*aux22+aux21*aux22D;
            aux2DD = aux21DD*aux22+aux21D*aux22D+aux21D*aux22D+aux21*aux22DD; 
            aux2DDD = aux21DDD*aux22+aux21DD*aux22D+aux21DD*aux22D+aux21D*aux22DD+aux21DD*aux22D+aux21D*aux22DD+aux21D*aux22DD+aux21*aux22DDD;
            aux3D = aux31D*aux32+aux31*aux32D;
            aux3DD = aux31DD*aux32+aux31D*aux32D+aux31D*aux32D+aux31*aux32DD;
            aux3DDD = aux31DDD*aux32+aux31DD*aux32D+aux31DD*aux32D+aux31D*aux32DD+aux31DD*aux32D+aux31D*aux32DD+aux31D*aux32DD+aux31*aux32DDD;
      

            auxT = zetaT[3]/((-1 + zeta[3])*(-1 + zeta[3]));
            auxTT = (-2*(zetaT[3]*zetaT[3]) + (-1 + zeta[3])*zetaTT[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
            aux21T = ((d[k]*d[k])*dDT[i] + (d[i]*d[i])*dDT[k])/((d[i] + d[k])*(d[i] + d[k]));
            aux21TT = ((d[k]*d[k]*d[k])*dDTT[i] + (d[k]*d[k])*(-2*(dDT[i]*dDT[i]) + d[i]*dDTT[i]) + 
            d[i]*d[k]*(4*dDT[i]*dDT[k] + d[i]*dDTT[k]) + (d[i]*d[i])*(-2*(dDT[k]*dDT[k]) + d[i]*dDTT[k]))/
            ((d[i] + d[k])*(d[i] + d[k])*(d[i] + d[k]));
            aux22T = (3*(-1 + zeta[3])*zetaT[2] - 6*zeta[2]*zetaT[3])/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])); 
            aux22TT = (3*(-4*(-1 + zeta[3])*zetaT[2]*zetaT[3] + 6*zeta[2]*(zetaT[3]*zetaT[3]) + ((-1 + zeta[3])*(-1 + zeta[3]))*zetaTT[2] - 
            2*zeta[2]*(-1 + zeta[3])*zetaTT[3]))/((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));
            aux31T = (2*(d[i]*(d[k]*d[k]*d[k])*dDT[i] + (d[i]*d[i]*d[i])*d[k]*dDT[k]))/((d[i] + d[k])*(d[i] + d[k])*(d[i] + d[k])); 
            aux31TT = (2*((d[k]*d[k]*d[k]*d[k])*(dDT[i]*dDT[i]) + d[i]*(d[k]*d[k]*d[k])*(-2*(dDT[i]*dDT[i]) + d[k]*dDTT[i]) + 
            (d[i]*d[i])*(d[k]*d[k])*(6*dDT[i]*dDT[k] + d[k]*dDTT[i]) + (d[i]*d[i]*d[i])*d[k]*
            (-2*(dDT[k]*dDT[k]) + d[k]*dDTT[k]) + (d[i]*d[i]*d[i]*d[i])*((dDT[k]*dDT[k]) + d[k]*dDTT[k])))/
            ((d[i] + d[k])*(d[i] + d[k])*(d[i] + d[k])*(d[i] + d[k]));
            aux32T = (zeta[2]*((4. - 4.*zeta[3])*zetaT[2] + 6.*zeta[2]*zetaT[3]))/((-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])*(-1. + zeta[3])); 
            aux32TT = (-4.*((1. - zeta[3])*(1. - zeta[3]))*(zetaT[2]*zetaT[2]) + zeta[2]*(-24. + 24.*zeta[3])*zetaT[2]*zetaT[3] + 
            zeta[2]*(-4.*((1. - zeta[3])*(1. - zeta[3]))*zetaTT[2] + zeta[2]*(-24.*(zetaT[3]*zetaT[3]) + (-6. + 6.*zeta[3])*zetaTT[3])))
            /((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]));        
            aux2T = aux21T*aux22+aux21*aux22T;
            aux2TT = aux21TT*aux22+aux21T*aux22T+aux21T*aux22T+aux21*aux22TT;
            aux3T = aux31T*aux32+aux31*aux32T;
            aux3TT = aux31TT*aux32+aux31T*aux32T+aux31T*aux32T+aux31*aux32TT;

            auto ghs_ki = aux+aux2+aux3;
            auto dghs_ki = (auxD+aux2D+aux3D)/num_denD;
            auto d2ghs_ki = (auxDD+aux2DD+aux3DD)/num_denD/num_denD;
            auto d3ghs_ki = (auxDDD+aux2DDD+aux3DDD)/num_denD/num_denD/num_denD;
            auto Tghs_ki = auxT+aux2T+aux3T;
            auto T2ghs_ki = auxTT+aux2TT+aux3TT;
            
            auto exponential = exp(eassoc_ki/T);
            delta(l,k,j,i) = kappa_ki*(exponential - 1.)*ghs_ki*sigma_ki*sigma_ki*sigma_ki; // EQ. A13, Huang and Radosz (1991)  
            deltaD(l,k,j,i) = delta(l,k,j,i)/ghs_ki*dghs_ki; 
            deltaDD(l,k,j,i) = delta(l,k,j,i)/ghs_ki*d2ghs_ki;
            deltaDDD(l,k,j,i) = delta(l,k,j,i)/ghs_ki*d3ghs_ki;
            deltaT(l,k,j,i) = kappa_ki*sigma_ki*sigma_ki*sigma_ki*
              ((exponential - 1.)*Tghs_ki-exponential*eassoc_ki/T/T*ghs_ki);
            deltaTT(l,k,j,i) = kappa_ki*sigma_ki*sigma_ki*sigma_ki*((exponential*eassoc_ki*(eassoc_ki + 2*T)*ghs_ki)/T/T/T/T - (2*exponential*eassoc_ki*Tghs_ki)/T/T + 
              (-1 + exponential)*T2ghs_ki);              
          }
        }
      }
    } 
             
    auto C = 6.*zeta[2]/pi/num_den;
    auto D = 6.*zeta[3]/pi/num_den;    
    for (int l = 0; l < s; l++) {
      for (int k = 0; k < n_comp; k++) {
        for (int j = 0; j < s; j++){
          for (int i = 0; i < n_comp; i++){
            for (int mm = 0; mm < n_comp; mm++){
              if (l==j) continue;
              eassoc_ki = (eassoc(k)+eassoc(i))/2.; // EQ. 2 Association term paper Gross, Sadowski 
              aux = sqrt(sigma(k)*sigma(i))/(0.5*(sigma(k)+sigma(i)));
              kappa_ki = sqrt(kappa(k)*kappa(i))*aux*aux*aux; // EQ. 3 Association term paper Gross, Sadowski 
              sigma_ki = (sigma(k)+sigma(i))/2.;
  
              // Calculation of auxiliary variables to obtain delta_xm and its derivatives
              auto dghs_xm = num_den*pi/6.*m(mm)*((d[mm]*d[mm]*d[mm])/(1-zeta[3])/(1-zeta[3]) + 3*d[k]*d[i]/
                    (d[k]+d[i])*(d[mm]*d[mm]/(1-zeta[3])/(1-zeta[3])+2*(d[mm]*d[mm]*d[mm])*
                    zeta[2]/((1-zeta[3])*(1-zeta[3])*(1-zeta[3]))) + 2*((d[k]*d[i]/(d[k]+d[i]))*(d[k]*d[i]/(d[k]+d[i])))*
                    (2*d[mm]*d[mm]*zeta[2]/((1-zeta[3])*(1-zeta[3])*(1-zeta[3]))+3*(d[mm]*d[mm]*d[mm])*zeta[2]*zeta[2]
                    /((1-zeta[3])*(1-zeta[3])*(1-zeta[3])*(1-zeta[3])))); // EQ. A88, Huang and Radosz (1991)
                    
              auto d2ghs_xm_rho = (-6*pi*(d[mm]*d[mm])*(((-6 + D*pi*num_den)*(-6 + D*pi*num_den))*(6 + D*pi*num_den)*
                  (d[k]*d[k])*d[mm] + (-6 + D*pi*num_den)*d[i]*d[k]*(2*(-36 + (D*D)*(pi*pi)*(num_den*num_den))*d[mm] + 
                  3*d[k]*(-36 + (D*D)*(pi*pi)*(num_den*num_den) - 2*C*pi*num_den*(12 + D*pi*num_den)*d[mm])) + 
                  (d[i]*d[i])*(((-6 + D*pi*num_den)*(-6 + D*pi*num_den))*(6 + D*pi*num_den)*d[mm] + 3*(-6 + D*pi*num_den)*d[k]*(-36 + (D*D)*(pi*pi)*(num_den*num_den) - 
                  2*C*pi*num_den*(12 + D*pi*num_den)*d[mm]) + 2*C*pi*num_den*(d[k]*d[k])*(-2*(-72 + 6*D*pi*num_den + (D*D)*(pi*pi)*(num_den*num_den)) + 
                  3*C*pi*num_den*(18 + D*pi*num_den)*d[mm])))*m(mm))/(((-6 + D*pi*num_den)*(-6 + D*pi*num_den)*(-6 + D*pi*num_den)*(-6 + D*pi*num_den)*(-6 + D*pi*num_den))*((d[i] + d[k])*(d[i] + d[k]))); // derivative of dghs_xm wrt numerical density  
                  
                  
              auto d3ghs_xm_rho_rho = (12*(pi*pi)*(d[mm]*d[mm])*(D*((-6 + D*pi*num_den)*(-6 + D*pi*num_den))*(12 + D*pi*num_den)*(d[k]*d[k])*d[mm] + 
                  (-6 + D*pi*num_den)*d[i]*d[k]*(2*D*(-72 + 6*D*pi*num_den + (D*D)*(pi*pi)*(num_den*num_den))*d[mm] + 
                  3*d[k]*(D*(-72 + 6*D*pi*num_den + (D*D)*(pi*pi)*(num_den*num_den)) - 
                  2*C*(36 + 24*D*pi*num_den + (D*D)*(pi*pi)*(num_den*num_den))*d[mm])) + 
                  (d[i]*d[i])*(D*((-6 + D*pi*num_den)*(-6 + D*pi*num_den))*(12 + D*pi*num_den)*d[mm] + 
                  3*(-6 + D*pi*num_den)*d[k]*(D*(-72 + 6*D*pi*num_den + (D*D)*(pi*pi)*(num_den*num_den)) - 
                  2*C*(36 + 24*D*pi*num_den + (D*D)*(pi*pi)*(num_den*num_den))*d[mm]) + 
                  2*C*(d[k]*d[k])*(432 + 216*D*pi*num_den - 36*(D*D)*(pi*pi)*(num_den*num_den) - 
                  2*(D*D*D)*(pi*pi*pi)*(num_den*num_den*num_den) + 3*C*pi*num_den*(108 + 36*D*pi*num_den + (D*D)*(pi*pi)*(num_den*num_den))*d[mm])))*m(mm))/
                  (((-6 + D*pi*num_den)*(-6 + D*pi*num_den)*(-6 + D*pi*num_den)*(-6 + D*pi*num_den)*(-6 + D*pi*num_den)*(-6 + D*pi*num_den))*((d[i] + d[k])*(d[i] + d[k]))); // derivative of dghs_xm_rho wrt numerical density 

              auto d4ghs_xm_rho_rho_rho = 1./num_denD/num_denD/num_denD*((d[mm]*d[mm])*m(mm)*pi*(6*((-1 + zeta[3])*(-1 + zeta[3]))*
                  (((d[i] + d[k])*(d[i] + d[k]))*d[mm]*((-1 + zeta[3])*(-1 + zeta[3]))*zetaD[3] + 2*(d[i]*d[i])*(d[k]*d[k])*(-((1 + 3*d[mm]*zeta[2] - zeta[3])*(-1 + zeta[3])*
                  zetaD[2]) + 3*zeta[2]*(1 + 2*d[mm]*zeta[2] - zeta[3])*zetaD[3]) + 3*d[i]*d[k]*(d[i] + d[k])*(1 - zeta[3])*
                  ((d[mm] - d[mm]*zeta[3])*zetaD[2] + (1 + 3*d[mm]*zeta[2] - zeta[3])*zetaD[3]))*num_denDD + 
                  6*(1 - zeta[3])*num_denD*(((d[i] + d[k])*(d[i] + d[k]))*d[mm]*((-1 + zeta[3])*(-1 + zeta[3]))*
                  (3*(zetaD[3]*zetaD[3]) - (-1 + zeta[3])*zetaDD[3]) + 3*d[i]*d[k]*(d[i] + d[k])*(1 - zeta[3])*
                  (-6*d[mm]*(-1 + zeta[3])*zetaD[2]*zetaD[3] + 3*(1 + 4*d[mm]*zeta[2] - zeta[3])*(zetaD[3]*zetaD[3]) + 
                  (-1 + zeta[3])*(d[mm]*(-1 + zeta[3])*zetaDD[2] + (-1 - 3*d[mm]*zeta[2] + zeta[3])*zetaDD[3])) + 
                  (d[i]*d[i])*(d[k]*d[k])*(2*(1 - zeta[3])*(-6*(-1 + zeta[3])*zetaD[2]*zetaD[3] + 
                  ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDD[2] + 3*zeta[2]*(4*(zetaD[3]*zetaD[3]) - 
                  (-1 + zeta[3])*zetaDD[3])) + 3*d[mm]*(-16*zeta[2]*(-1 + zeta[3])*zetaD[2]*zetaD[3] + 
                  2*((-1 + zeta[3])*(-1 + zeta[3]))*((zetaD[2]*zetaD[2]) + zeta[2]*zetaDD[2]) + 
                  4*(zeta[2]*zeta[2])*(5*(zetaD[3]*zetaD[3]) - (-1 + zeta[3])*zetaDD[3])))) + 
                  (2*(d[i]*d[i])*(d[k]*d[k])*zeta[2]*(2 + 3*d[mm]*zeta[2] - 2*zeta[3]) - 3*d[i]*d[k]*(d[i] + d[k])*(1 + 2*d[mm]*zeta[2] - zeta[3])*(-1 + zeta[3]) + 
                  ((d[i] + d[k])*(d[i] + d[k]))*d[mm]*((-1 + zeta[3])*(-1 + zeta[3])))*((1 - zeta[3])*(1 - zeta[3])*(1 - zeta[3]))*num_denDDD + 
                  2*num_den*(((d[i] + d[k])*(d[i] + d[k]))*d[mm]*((-1 + zeta[3])*(-1 + zeta[3]))*
                  (12*(zetaD[3]*zetaD[3]*zetaD[3]) - 9*(-1 + zeta[3])*zetaD[3]*zetaDD[3] + ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDDD[3]) + 
                  3*d[i]*d[k]*(d[i] + d[k])*(1 - zeta[3])*((1 - zeta[3])*(12*(zetaD[3]*zetaD[3]*zetaD[3]) - 
                  9*(-1 + zeta[3])*zetaD[3]*zetaDD[3] + ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDDD[3]) + 
                  d[mm]*(-9*(-1 + zeta[3])*zetaD[2]*(4*(zetaD[3]*zetaD[3]) - (-1 + zeta[3])*zetaDD[3]) + 
                  ((-1 + zeta[3])*(-1 + zeta[3]))*(9*zetaD[3]*zetaDD[2] - (-1 + zeta[3])*zetaDDD[2]) + 
                  3*zeta[2]*(20*(zetaD[3]*zetaD[3]*zetaD[3]) - 12*(-1 + zeta[3])*zetaD[3]*zetaDD[3] + 
                  ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDDD[3]))) + (d[i]*d[i])*(d[k]*d[k])*(3*d[mm]*(24*((-1 + zeta[3])*(-1 + zeta[3]))*zetaD[3]*
                  ((zetaD[2]*zetaD[2]) + zeta[2]*zetaDD[2]) + 24*zeta[2]*(1 - zeta[3])*zetaD[2]*
                  (5*(zetaD[3]*zetaD[3]) - (-1 + zeta[3])*zetaDD[3]) - 2*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*(3*zetaD[2]*zetaDD[2] + 
                  zeta[2]*zetaDDD[2]) + 4*(zeta[2]*zeta[2])*(30*(zetaD[3]*zetaD[3]*zetaD[3]) - 15*(-1 + zeta[3])*zetaD[3]*zetaDD[3] + 
                  ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDDD[3])) + 2*(1 - zeta[3])*(-9*(-1 + zeta[3])*zetaD[2]*
                  (4*(zetaD[3]*zetaD[3]) - (-1 + zeta[3])*zetaDD[3]) + ((-1 + zeta[3])*(-1 + zeta[3]))*(9*zetaD[3]*zetaDD[2] - 
                  (-1 + zeta[3])*zetaDDD[2]) + 3*zeta[2]*(20*(zetaD[3]*zetaD[3]*zetaD[3]) - 12*(-1 + zeta[3])*zetaD[3]*zetaDD[3] + 
                  ((-1 + zeta[3])*(-1 + zeta[3]))*zetaDDD[3]))))))/(6.*((d[i] + d[k])*(d[i] + d[k]))*((1 - zeta[3])*(1 - zeta[3])*(1 - zeta[3])*(1 - zeta[3])*(1 - zeta[3])*(1 - zeta[3])*(1 - zeta[3])));

              auto d2ghs_xm_T = (m(mm)*pi*num_den*d[mm]*(-4*d[i]*(d[k]*d[k])*(d[i] + d[k])*d[mm]*zeta[2]*
                  (2 + 3*d[mm]*zeta[2] - 2*zeta[3])*(-1 + zeta[3])*dDT[i] + 3*d[k]*((d[i] + d[k])*(d[i] + d[k]))*d[mm]*(1 + 2*d[mm]*zeta[2] - zeta[3])*((-1 + zeta[3])*(-1 + zeta[3]))*
                  dDT[i] - 4*(d[i]*d[i])*d[k]*(d[i] + d[k])*d[mm]*zeta[2]*(2 + 3*d[mm]*zeta[2] - 2*zeta[3])*(-1 + zeta[3])*dDT[k] + 
                  3*d[i]*((d[i] + d[k])*(d[i] + d[k]))*d[mm]*(1 + 2*d[mm]*zeta[2] - zeta[3])*((-1 + zeta[3])*(-1 + zeta[3]))*
                  dDT[k] + 4*(d[i]*d[i])*(d[k]*d[k])*d[mm]*zeta[2]*(2 + 3*d[mm]*zeta[2] - 2*zeta[3])*(-1 + zeta[3])*(dDT[i] + dDT[k]) - 
                  3*d[i]*d[k]*(d[i] + d[k])*d[mm]*(1 + 2*d[mm]*zeta[2] - zeta[3])*((-1 + zeta[3])*(-1 + zeta[3]))*
                  (dDT[i] + dDT[k]) - 3*((d[i] + d[k])*(d[i] + d[k])*(d[i] + d[k]))*d[mm]*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*dDT[mm] + 
                  2*((d[i] + d[k])*(d[i] + d[k])*(d[i] + d[k]))*(d[mm]*d[mm])*((-1 + zeta[3])*(-1 + zeta[3]))*zetaT[3] + 2*(d[i]*d[i])*(d[k]*d[k])*(d[i] + d[k])*
                  (2*d[mm]*((-1 + zeta[3])*(-1 + zeta[3]))*zetaT[2] + 3*d[mm]*(zeta[2]*zeta[2])*(-3*(-1 + zeta[3])*dDT[mm] + 4*d[mm]*zetaT[3]) + 
                  2*zeta[2]*(-1 + zeta[3])*(2*(-1 + zeta[3])*dDT[mm] - 3*d[mm]*(d[mm]*zetaT[2] + zetaT[3]))) + 
                  6*d[i]*d[k]*((d[i] + d[k])*(d[i] + d[k]))*(1 - zeta[3])*(-((1 + 3*d[mm]*zeta[2] - zeta[3])*(-1 + zeta[3])*dDT[mm]) + 
                  d[mm]*(-((-1 + zeta[3])*zetaT[3]) + d[mm]*(-((-1 + zeta[3])*zetaT[2]) + 3*zeta[2]*zetaT[3])))))/
                  (6.*((d[i] + d[k])*(d[i] + d[k])*(d[i] + d[k]))*((1 - zeta[3])*(1 - zeta[3])*(1 - zeta[3])*(1 - zeta[3])*(1 - zeta[3]))); // derivative of dghs_xm wrt Temperature 
                  
              auto d3ghs_xm_T_T = (m(mm)*pi*num_den*(-12*(d[mm]*d[mm])*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*dDT[mm]*zetaT[3] - 
                  (8*d[i]*d[k]*d[mm]*(1 - zeta[3])*(d[k]*dDT[i] + d[i]*dDT[k])*(-2*d[mm]*zeta[2]*(2 + 3*d[mm]*zeta[2] - 2*zeta[3])*(-1 + zeta[3])*
                  (dDT[i] + dDT[k]) - (d[i] + d[k])*(2*d[mm]*((-1 + zeta[3])*(-1 + zeta[3]))*zetaT[2] + 3*d[mm]*(zeta[2]*zeta[2])*(-3*(-1 + zeta[3])*dDT[mm] + 
                  4*d[mm]*zetaT[3]) + 2*zeta[2]*(-1 + zeta[3])*(2*(-1 + zeta[3])*dDT[mm] - 3*d[mm]*(d[mm]*zetaT[2] + zetaT[3])))))/((d[i] + d[k])*(d[i] + d[k])*(d[i] + d[k]))
                  - (6*d[mm]*((-1 + zeta[3])*(-1 + zeta[3]))*(d[k]*dDT[i] + d[i]*dDT[k])*(-(d[mm]*(1 + 2*d[mm]*zeta[2] - zeta[3])*(-1 + zeta[3])*
                  (dDT[i] + dDT[k])) - 2*(d[i] + d[k])*(-((1 + 3*d[mm]*zeta[2] - zeta[3])*(-1 + zeta[3])*dDT[mm]) + 
                  d[mm]*(-((-1 + zeta[3])*zetaT[3]) + d[mm]*(-((-1 + zeta[3])*zetaT[2]) + 3*zeta[2]*zetaT[3])))))/
                  ((d[i] + d[k])*(d[i] + d[k])) + (3*(d[mm]*d[mm])*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*(-1 - 2*d[mm]*zeta[2] + zeta[3])*(2*dDT[i]*dDT[k] + d[k]*dDTT[i] + 
                  d[i]*dDTT[k]))/(d[i] + d[k]) + (4*(d[mm]*d[mm])*zeta[2]*(2 + 3*d[mm]*zeta[2] - 2*zeta[3])*((-1 + zeta[3])*(-1 + zeta[3]))*
                  ((d[i]*d[i])*(dDT[k]*dDT[k]) + (d[k]*d[k])*((dDT[i]*dDT[i]) + d[i]*dDTT[i]) + d[i]*d[k]*(4*dDT[i]*dDT[k] + d[i]*dDTT[k])))/
                  ((d[i] + d[k])*(d[i] + d[k])) + 3*d[mm]*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*(2*
                  (dDT[mm]*dDT[mm]) + d[mm]*dDTT[mm]) + 2*(d[mm]*d[mm]*d[mm])*((-1 + zeta[3])*(-1 + zeta[3]))*
                  (3*(zetaT[3]*zetaT[3]) - (-1 + zeta[3])*zetaTT[3]) - (3*d[i]*d[k]*(1 - zeta[3])*(4*(d[i] + d[k])*d[mm]*(1 - zeta[3])*
                  (dDT[i] + dDT[k])*(-((1 + 3*d[mm]*zeta[2] - zeta[3])*(-1 + zeta[3])*dDT[mm]) + d[mm]*(-((-1 + zeta[3])*zetaT[3]) + 
                  d[mm]*(-((-1 + zeta[3])*zetaT[2]) + 3*zeta[2]*zetaT[3]))) - (d[mm]*d[mm])*(1 + 2*d[mm]*zeta[2] - zeta[3])*((-1 + zeta[3])*(-1 + zeta[3]))*
                  (2*((dDT[i] + dDT[k])*(dDT[i] + dDT[k])) - (d[i] + d[k])*(dDTT[i] + dDTT[k])) - 2*((d[i] + d[k])*(d[i] + d[k]))*(4*d[mm]*((-1 + zeta[3])*(-1 + zeta[3]))*dDT[mm]*
                  zetaT[3] + 6*(d[mm]*d[mm])*(1 - zeta[3])*dDT[mm]*(-((-1 + zeta[3])*zetaT[2]) + 3*zeta[2]*zetaT[3]) - ((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*((dDT[mm]*dDT[mm]) + d[mm]*dDTT[mm]) + 
                  3*d[mm]*zeta[2]*((-1 + zeta[3])*(-1 + zeta[3]))*(2*(dDT[mm]*dDT[mm]) + d[mm]*dDTT[mm]) + (d[mm]*d[mm])*(1 - zeta[3])*
                  (3*(zetaT[3]*zetaT[3]) - (-1 + zeta[3])*zetaTT[3]) + (d[mm]*d[mm]*d[mm])*(-6*(-1 + zeta[3])*zetaT[2]*zetaT[3] + 
                  ((-1 + zeta[3])*(-1 + zeta[3]))*zetaTT[2] + 3*zeta[2]*(4*(zetaT[3]*zetaT[3]) - (-1 + zeta[3])*zetaTT[3])))))/
                  ((d[i] + d[k])*(d[i] + d[k])*(d[i] + d[k])) - (2*(d[i]*d[i])*(d[k]*d[k])*(4*(d[i] + d[k])*d[mm]*(1 - zeta[3])*(dDT[i] + dDT[k])*
                  (2*d[mm]*((-1 + zeta[3])*(-1 + zeta[3]))*zetaT[2] + 3*d[mm]*(zeta[2]*zeta[2])*(-3*(-1 + zeta[3])*dDT[mm] + 
                  4*d[mm]*zetaT[3]) + 2*zeta[2]*(-1 + zeta[3])*(2*(-1 + zeta[3])*dDT[mm] - 3*d[mm]*(d[mm]*zetaT[2] + zetaT[3]))) - 
                  2*(d[mm]*d[mm])*zeta[2]*(2 + 3*d[mm]*zeta[2] - 2*zeta[3])*((-1 + zeta[3])*(-1 + zeta[3]))*(3*((dDT[i] + dDT[k])*(dDT[i] + dDT[k])) - 
                  (d[i] + d[k])*(dDTT[i] + dDTT[k])) - ((d[i] + d[k])*(d[i] + d[k]))*(-36*(d[mm]*d[mm])*zeta[2]*(-1 + zeta[3])*dDT[mm]*
                  (-((-1 + zeta[3])*zetaT[2]) + 2*zeta[2]*zetaT[3]) + 8*d[mm]*((-1 + zeta[3])*(-1 + zeta[3]))*dDT[mm]*(-((-1 + zeta[3])*zetaT[2]) + 3*zeta[2]*zetaT[3]) - 
                  4*zeta[2]*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*((dDT[mm]*dDT[mm]) + d[mm]*dDTT[mm]) + 9*d[mm]*(zeta[2]*zeta[2])*((-1 + zeta[3])*(-1 + zeta[3]))*
                  (2*(dDT[mm]*dDT[mm]) + d[mm]*dDTT[mm]) + 2*(d[mm]*d[mm])*(1 - zeta[3])*(-6*(-1 + zeta[3])*zetaT[2]*zetaT[3] + 
                  ((-1 + zeta[3])*(-1 + zeta[3]))*zetaTT[2] + 3*zeta[2]*(4*(zetaT[3]*zetaT[3]) - (-1 + zeta[3])*zetaTT[3])) + 
                  3*(d[mm]*d[mm]*d[mm])*(-16*zeta[2]*(-1 + zeta[3])*zetaT[2]*zetaT[3] + 2*((-1 + zeta[3])*(-1 + zeta[3]))*((zetaT[2]*zetaT[2]) + 
                  zeta[2]*zetaTT[2]) + 4*(zeta[2]*zeta[2])*(5*(zetaT[3]*zetaT[3]) - (-1 + zeta[3])*zetaTT[3])))))/((d[i] + d[k])*(d[i] + d[k])*(d[i] + d[k])*(d[i] + d[k]))))/
                  (6.*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))); // derivative of dghs_xm_T wrt Temperature 

              auto exponential = exp(eassoc_ki/T); 
              delta_xm(l,k,j,i,mm) = kappa_ki*(exponential - 1.)*dghs_xm*sigma_ki*sigma_ki*sigma_ki; // EQ. A87, Huang and Radosz (1991)
              delta_xm_rho(l,k,j,i,mm) = delta_xm(l,k,j,i,mm)/dghs_xm*d2ghs_xm_rho;
              delta_xm_rho_rho(l,k,j,i,mm) = delta_xm(l,k,j,i,mm)/dghs_xm*d3ghs_xm_rho_rho;
              delta_xm_rho_rho_rho(l,k,j,i,mm) = delta_xm(l,k,j,i,mm)/dghs_xm*d4ghs_xm_rho_rho_rho;
              delta_xm_T(l,k,j,i,mm) = kappa_ki*sigma_ki*sigma_ki*sigma_ki*
                ((exponential - 1.)*d2ghs_xm_T-exponential*eassoc_ki/T/T*dghs_xm);
              delta_xm_T_T(l,k,j,i,mm) = kappa_ki*sigma_ki*sigma_ki*sigma_ki*
                ((exponential*eassoc_ki*(eassoc_ki + 2*T)*dghs_xm)/T/T/T/T - 
                (2*exponential*eassoc_ki*d2ghs_xm_T)/T/T + (-1 + exponential)*d3ghs_xm_T_T);                
            } 
          }
        }
      }
    } 
         
    calculate_XA(n_comp,x,num_den,indx);

    //Calculating the elements of the matrix lambda(p,q) (Tan et al. (2004))
    int p = -1;
    int q = -1;
    for (int j = 0; j < s; j++) {
      for (int i = 0; i < n_comp; i++) {
        ++p;
        for (int l = 0; l < s; l++){
          for (int k = 0; k < n_comp; k++){
            ++q;
            if(l==j) {
              if (k==i) {lambda(p,q) = 1.;} else {lambda(p,q) = 0.;} continue; 
            } 
            lambda(p,q) = x[k]*num_den*S(l,k)*delta(l,k,j,i)*XA(j,i)*XA(j,i);
          } 
        }
        q = -1;
      }
    }

    //Inverting lambda(p,q) using LU decomposition (Tan et al. (2004))  
    detail::ludcmp(lambda,order,indx,aux);
    for(int j=0;j<order;j++) {
      for(int i=0;i<order;i++) col[i]=0.;
      col[j]=1.;
      detail::lubksb(lambda,order,indx,col);
      for(int i=0;i<order;i++) inv_lambda(i,j)=col[i];
    }   

    // Calculating psi_T (lambda * XAT = psi_T) (Tan et al. (2004))
    p = -1;
    for (int j = 0; j < s; j++) {
      for (int i = 0; i < n_comp; i++) {
        sum1 = 0.;
        for (int k = 0; k < n_comp; k++){
          sum2 = 0.;
          for (int l = 0; l < s; l++){
            if(l==j) continue;
            sum2 += S(l,k)*XA(l,k)*deltaT(l,k,j,i);
          } 
          sum1 += x[k]*sum2;
        }
        psi[++p] = -XA(j,i)*XA(j,i)*num_den*sum1;
      }
    } 

    // Creating and Calculating XAT matrix (derivative of XA wrt Temperature) (Tan et al. (2004))     
    std::fill(dXA.begin(),dXA.end(),0.);
    for (p = 0; p < order; p++) {
      for (q = 0; q < order; q++) {
        dXA[p] += inv_lambda(p,q)*psi[q];
      }
    }
    p = -1;
    for (int j = 0; j < s; j++) {
      for (int i = 0; i < n_comp; i++) {
        XAT(j,i) = dXA[++p];
      }
    }

    // Calculating psi_T_T (lambda * XATT = psi_T_T)  (Tan et al. (2004)
    p = -1;
    for (int j = 0; j < s; j++) {
      for (int i = 0; i < n_comp; i++) {
        sum1 = 0.;
        for (int k = 0; k < n_comp; k++){
          sum2 = 0.;
          for (int l = 0; l < s; l++){
            if(l==j) continue;
            sum2 += S(l,k)*(XA(l,k)*deltaTT(l,k,j,i)+2.*XAT(l,k)*deltaT(l,k,j,i)); 
          } 
          sum1 += x[k]*sum2;
        }
        psi[++p] = 2./XA(j,i)*XAT(j,i)*XAT(j,i)-XA(j,i)*XA(j,i)*num_den*sum1; 
      }
    } 

    // Creating and Calculating XATT matrix (second derivative of XA wrt Temperature) (Tan et al. (2004))
    std::fill(dXA.begin(),dXA.end(),0.);
    for (p = 0; p < order; p++) {
      for (q = 0; q < order; q++) {
        dXA[p] += inv_lambda(p,q)*psi[q];
      }
    }
    p = -1;
    for (int j = 0; j < s; j++) {
      for (int i = 0; i < n_comp; i++) {
        XATT(j,i) = dXA[++p];
      }
    }    

    // Calculating psi_rho (lambda * XAD = psi_rho) (Tan et al. (2004))
    p = -1;
    for (int j = 0; j < s; j++) {
      for (int i = 0; i < n_comp; i++) {
        sum1 = 0.;
        for (int k = 0; k < n_comp; k++){
          sum2 = 0.;
          for (int l = 0; l < s; l++){
            if(l==j) continue;
            sum2 += S(l,k)*XA(l,k)*deltaD(l,k,j,i);
          } 
          sum1 += x[k]*sum2;
        }
        psi[++p] = -XA(j,i)*XA(j,i)*(num_den*sum1+(1./XA(j,i)-1.)/num_den);
      }
    }           

    // Creating and Calculating XAD matrix (derivative of XA wrt rho (not rhoi)) (Tan et al. (2004)) 
    std::fill(dXA.begin(),dXA.end(),0.);
    for (p = 0; p < order; p++) {
      for (q = 0; q < order; q++) {
        dXA[p] += inv_lambda(p,q)*psi[q];
      }
    }
    p = -1;
    for (int j = 0; j < s; j++) {
      for (int i = 0; i < n_comp; i++) {
        XAD(j,i) = dXA[++p];
      }
    }  

    // Calculating psi_rho_rho (lambda * XADD = psi_rho_rho) (Tan et al. (2004))
    p = -1;
    for (int j = 0; j < s; j++) {
      for (int i = 0; i < n_comp; i++) {
        sum1 = 0.;
        for (int k = 0; k < n_comp; k++){
          sum2 = 0.;
          for (int l = 0; l < s; l++){
            if(l==j) continue;
            sum2 += S(l,k)*(XA(l,k)*deltaDD(l,k,j,i)+2.*XAD(l,k)*deltaD(l,k,j,i));
          } 
          sum1 += x[k]*sum2;
        }
        psi[++p] = 2./XA(j,i)*XAD(j,i)*XAD(j,i)+2./num_den*XAD(j,i)+
          XA(j,i)*XA(j,i)*(2.*(1./XA(j,i)-1.)/num_den/num_den-num_den*sum1);
      }
    }  

    // Creating and Calculating XADD matrix (derivative of XA wrt rho wrt rho ) (Tan et al. (2004)) 
    std::fill(dXA.begin(),dXA.end(),0.);
    for (p = 0; p < order; p++) {
      for (q = 0; q < order; q++) {
        dXA[p] += inv_lambda(p,q)*psi[q];
      }
    }
    p = -1;
    for (int j = 0; j < s; j++) {
      for (int i = 0; i < n_comp; i++) {
        XADD(j,i) = dXA[++p];
      }
    } 


    // Calculating psi_rho_rho_rho (lambda * XADDD = psi_rho_rho_rho)
    p = -1;
    for (int j = 0; j < s; j++) {
      for (int i = 0; i < n_comp; i++) {
        sum1 = 0.;
        for (int k = 0; k < n_comp; k++){
          sum2 = 0.;
          for (int l = 0; l < s; l++){
            if(l==j) continue;
            sum2 += S(l,k)*(XA(l,k)*deltaDDD(l,k,j,i)+3.*XAD(l,k)*deltaDD(l,k,j,i)+
              3.*XADD(l,k)*deltaD(l,k,j,i));
          } 
          sum1 += x[k]*sum2;
        }
        psi[++p] = 3./num_den*XADD(j,i)-6./num_den/num_den/num_den*XA(j,i)*(1.-XA(j,i))
          +6./XA(j,i)*(XAD(j,i)*XADD(j,i)-1./XA(j,i)*XAD(j,i)*XAD(j,i)*XAD(j,i))
          -6./num_den*XAD(j,i)*(1./num_den+1./XA(j,i)*XAD(j,i))-num_den*XA(j,i)*XA(j,i)*sum1;
      }
    }  

    // Creating and Calculating XADDD matrix (derivative of XA wrt rho wrt rho wrt rho ) 
    std::fill(dXA.begin(),dXA.end(),0.);
    for (p = 0; p < order; p++) {
      for (q = 0; q < order; q++) {
        dXA[p] += inv_lambda(p,q)*psi[q];
      }
    }
    p = -1;
    for (int j = 0; j < s; j++) {
      for (int i = 0; i < n_comp; i++) {
        XADDD(j,i) = dXA[++p];
      }
    } 
  
    // Calculating psi_xm (lambda * XA_xm = psi__xm) (Tan et al. (2004))
    for (int mm = 0; mm < n_comp; mm++) {
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++){
          sum1 = 0;
          for (int l = 0; l < s; l++){
            if (l==j) continue;
            sum2 = 0;
            for (int k = 0; k < n_comp; k++){              
              sum2 += x[k]*S(l,k)*XA(l,k)*delta_xm(l,k,j,i,mm);
            } 
            sum1 += S(l,mm)*XA(l,mm)*delta(l,mm,j,i) + sum2;
          }
          psi2D(++p,mm) = -XA(j,i)*XA(j,i)*num_den*sum1;
        }
      }
    }         


    // Creating and Calculating XA_xm matrix (derivative of XA wrt xm) (Tan et al. (2004))
    for (int mm = 0; mm < n_comp; mm++) {
      std::fill(dXA.begin(),dXA.end(),0.);
      for (p = 0; p < order; p++) {
        for (q = 0; q < order; q++) {
          dXA[p] += inv_lambda(p,q)*psi2D(q,mm);
        }
      }      
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++) {
          XA_xm(j,i,mm) = dXA[++p];
        }
      }         
    }    

    // Calculating psi_xm_T (lambda * XA_xm_T = psi_xm_T)
    for (int mm = 0; mm < n_comp; mm++) {
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++){
          sum1 = 0;
          for (int l = 0; l < s; l++){
            if (l==j) continue;
            sum2 = 0;
            for (int k = 0; k < n_comp; k++){              
              sum2 += x[k]*S(l,k)*(XAT(l,k)*delta_xm(l,k,j,i,mm)+XA(l,k)*delta_xm_T(l,k,j,i,mm)+
                XA_xm(l,k,mm)*deltaT(l,k,j,i));
            } 
            sum1 += S(l,mm)*(XAT(l,mm)*delta(l,mm,j,i)+XA(l,mm)*deltaT(l,mm,j,i)) + sum2;
          }
          psi2D(++p,mm) = 2.*XAT(j,i)*XA_xm(j,i,mm)/XA(j,i)-XA(j,i)*XA(j,i)*num_den*sum1;
        }
      }
    }     

    // Creating and Calculating XA_xm_T matrix (derivative of XA wrt xm wrt T)
    for (int mm = 0; mm < n_comp; mm++) {
      std::fill(dXA.begin(),dXA.end(),0.);
      for (p = 0; p < order; p++) {
        for (q = 0; q < order; q++) {
          dXA[p] += inv_lambda(p,q)*psi2D(q,mm);
        }
      }      
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++) {
          XA_xm_T(j,i,mm) = dXA[++p];
        }
      }         
    } 

    // Calculating psi_xm_T_T (lambda * XA_xm_T_T = psi_xm_T_T)
    for (int mm = 0; mm < n_comp; mm++) {
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++){
          sum1 = 0;
          for (int l = 0; l < s; l++){
            if (l==j) continue;
            sum2 = 0;
            for (int k = 0; k < n_comp; k++){              
              sum2 += x[k]*S(l,k)*(XATT(l,k)*delta_xm(l,k,j,i,mm)+2.*XAT(l,k)*delta_xm_T(l,k,j,i,mm)
                + XA(l,k)*delta_xm_T_T(l,k,j,i,mm)+2.*XA_xm_T(l,k,mm)*deltaT(l,k,j,i)
                + XA_xm(l,k,mm)*deltaTT(l,k,j,i)); 
            } 
            sum1 += S(l,mm)*(XATT(l,mm)*delta(l,mm,j,i)+2.*XAT(l,mm)*deltaT(l,mm,j,i) 
              + XA(l,mm)*deltaTT(l,mm,j,i)) + sum2;
          }
          psi2D(++p,mm) = 4.*XAT(j,i)*XA_xm_T(j,i,mm)/XA(j,i)+2.*XATT(j,i)*XA_xm(j,i,mm)/XA(j,i)
           - 6.*XAT(j,i)*XAT(j,i)*XA_xm(j,i,mm)/XA(j,i)/XA(j,i)-XA(j,i)*XA(j,i)*num_den*sum1;
        }
      }
    }     

    // Creating and Calculating XA_xm_T_T matrix (derivative of XA wrt xm wrt T wrt T)
    for (int mm = 0; mm < n_comp; mm++) {
      std::fill(dXA.begin(),dXA.end(),0.);
      for (p = 0; p < order; p++) {
        for (q = 0; q < order; q++) {
          dXA[p] += inv_lambda(p,q)*psi2D(q,mm);
        }
      }      
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++) {
          XA_xm_T_T(j,i,mm) = dXA[++p];
        }
      }         
    } 

    // Calculating psi_xm_rho (lambda * XA_xm_rho = psi__xm_rho)
    for (int mm = 0; mm < n_comp; mm++) {
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++){
          sum1 = 0;
          for (int l = 0; l < s; l++){
            if (l==j) continue;
            sum2 = 0;
            for (int k = 0; k < n_comp; k++){              
              sum2 += x[k]*S(l,k)*(XA(l,k)*delta_xm_rho(l,k,j,i,mm)+XA_xm(l,k,mm)*deltaD(l,k,j,i)
              + XAD(l,k)*delta_xm(l,k,j,i,mm));
            } 
            sum1 += S(l,mm)*(XA(l,mm)*deltaD(l,mm,j,i) + XAD(l,mm)*delta(l,mm,j,i)) + sum2;
          }
          psi2D(++p,mm) = 2./XA(j,i)*XA_xm(j,i,mm)*XAD(j,i) + XA_xm(j,i,mm)/num_den
            -XA(j,i)*XA(j,i)*num_den*sum1;
        }
      }
    }

    // Creating and Calculating XA_xm_rho matrix (derivative of XA wrt xm wrt rho)
    for (int mm = 0; mm < n_comp; mm++) {
      std::fill(dXA.begin(),dXA.end(),0.);
      for (p = 0; p < order; p++) {
        for (q = 0; q < order; q++) {
          dXA[p] += inv_lambda(p,q)*psi2D(q,mm);
        }
      }      
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++) {
          XA_xm_rho(j,i,mm) = dXA[++p];
        }
      }         
    }    

    // Calculating psi_xm_rho_rho (lambda * XA_xm_rho_rho = psi__xm_rho_rho)
    for (int mm = 0; mm < n_comp; mm++) {
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++){
          sum1 = 0;
          for (int l = 0; l < s; l++){
            if (l==j) continue;
            sum2 = 0;
            for (int k = 0; k < n_comp; k++){              
              sum2 += x[k]*S(l,k)*(XA(l,k)*delta_xm_rho_rho(l,k,j,i,mm)+XADD(l,k)*delta_xm(l,k,j,i,mm)+
                2.*XAD(l,k)*delta_xm_rho(l,k,j,i,mm)+2.*XA_xm_rho(l,k,mm)*deltaD(l,k,j,i)+
                XA_xm(l,k,mm)*deltaDD(l,k,j,i));
            } 
            sum1 += S(l,mm)*(XADD(l,mm)*delta(l,mm,j,i) + 2.*XAD(l,mm)*deltaD(l,mm,j,i)+
              XA(l,mm)*deltaDD(l,mm,j,i)) + sum2;
          }
          psi2D(++p,mm) = -6./XA(j,i)/XA(j,i)*XAD(j,i)*XAD(j,i)*XA_xm(j,i,mm)+
            2./XA(j,i)*XADD(j,i)*XA_xm(j,i,mm)+4./XA(j,i)*XAD(j,i)*XA_xm_rho(j,i,mm)-
            2./num_den/num_den*XA_xm(j,i,mm)+2./num_den*XA_xm_rho(j,i,mm)-
            4./num_den/XA(j,i)*XAD(j,i)*XA_xm(j,i,mm)-XA(j,i)*XA(j,i)*num_den*sum1;
        }
      }
    }      

    // Creating and Calculating XA_xm_rho_rho matrix (derivative of XA wrt xm wrt rho wrt rho)
    for (int mm = 0; mm < n_comp; mm++) {
      std::fill(dXA.begin(),dXA.end(),0.);
      for (p = 0; p < order; p++) {
        for (q = 0; q < order; q++) {
          dXA[p] += inv_lambda(p,q)*psi2D(q,mm);
        }
      }      
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++) {
          XA_xm_rho_rho(j,i,mm) = dXA[++p];
        }
      }         
    }     


    // Calculating psi_xm_rho_rho_rho (lambda * XA_xm_rho_rho_rho = psi__xm_rho_rho_rho)
    for (int mm = 0; mm < n_comp; mm++) {
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++){
          sum1 = 0;
          for (int l = 0; l < s; l++){
            if (l==j) continue;
            sum2 = 0;
            for (int k = 0; k < n_comp; k++){              
              sum2 += x[k]*S(l,k)*(XA(l,k)*delta_xm_rho_rho_rho(l,k,j,i,mm)+3.*XAD(l,k)*delta_xm_rho_rho(l,k,j,i,mm)+
                XADDD(l,k)*delta_xm(l,k,j,i,mm)+3.*XADD(l,k)*delta_xm_rho(l,k,j,i,mm)+
                3.*XA_xm_rho_rho(l,k,mm)*deltaD(l,k,j,i)+3.*XA_xm_rho(l,k,mm)*deltaDD(l,k,j,i)+
                XA_xm(l,k,mm)*deltaDDD(l,k,j,i));
            } 
            sum1 += S(l,mm)*(XADDD(l,mm)*delta(l,mm,j,i) + 3.*XADD(l,mm)*deltaD(l,mm,j,i)+
              3.*XAD(l,mm)*deltaDD(l,mm,j,i)+XA(l,mm)*deltaDDD(l,mm,j,i)) + sum2;
          }
          psi2D(++p,mm) = 2./XA(j,i)*(XADDD(j,i)*XA_xm(j,i,mm)+3.*XADD(j,i)*XA_xm_rho(j,i,mm))
            +24./XA(j,i)/XA(j,i)/XA(j,i)*XAD(j,i)*XAD(j,i)*XAD(j,i)*XA_xm(j,i,mm)
            -18./XA(j,i)/XA(j,i)*XAD(j,i)*(XADD(j,i)*XA_xm(j,i,mm)+XAD(j,i)*XA_xm_rho(j,i,mm))
            +6./XA(j,i)*XAD(j,i)*XA_xm_rho_rho(j,i,mm) + 6./num_den/num_den/num_den*XA_xm(j,i,mm)
            -6./num_den/num_den*XA_xm_rho(j,i,mm)+3./num_den*XA_xm_rho_rho(j,i,mm) 
            +12./num_den/num_den/XA(j,i)*XAD(j,i)*XA_xm(j,i,mm)
            +18./num_den/XA(j,i)/XA(j,i)*XAD(j,i)*XAD(j,i)*XA_xm(j,i,mm)
            -6./num_den/XA(j,i)*XADD(j,i)*XA_xm(j,i,mm)
            -12./num_den/XA(j,i)*XAD(j,i)*XA_xm_rho(j,i,mm)-XA(j,i)*XA(j,i)*num_den*sum1;
        }
      }
    }      

    // Creating and Calculating XA_xm_rho_rho matrix (derivative of XA wrt xm wrt rho wrt rho wrt rho)
    for (int mm = 0; mm < n_comp; mm++) {
      std::fill(dXA.begin(),dXA.end(),0.);
      for (p = 0; p < order; p++) {
        for (q = 0; q < order; q++) {
          dXA[p] += inv_lambda(p,q)*psi2D(q,mm);
        }
      }      
      p = -1;
      for (int j = 0; j < s; j++) {
        for (int i = 0; i < n_comp; i++) {
          XA_xm_rho_rho_rho(j,i,mm) = dXA[++p];
        }
      }         
    }     

    // Calculating mu_assoc - EQ. A2 (Chapman et al. (1990))
    for (int i = 0; i < n_comp; i++) {
      sum1 = 0.;
      for (int j = 0; j < s; j++) {
        sum1 += S(j,i)*(log(XA(j,i)) - XA(j,i)/2. + 0.5);
      }
      
      sum2 = 0.;
      for (int k = 0; k < n_comp; k++) {
        sum3 = 0.;
        for (int j = 0; j < s; j++) {
          sum3 += S(j,k)*XA_xm(j,k,i)*(1./XA(j,k)-0.5);
        }
        sum2 += x[k]*sum3;
      }

      mu_assoc[i] = sum1+sum2;
    }

    // Calculating a_assoc - EQ. 1 (Tan et al. (2004)) and a_assocT
    sum1 = 0.;
    sum1T = 0.;
    sum1TT = 0.;
    for (int i = 0; i < n_comp; i++) {
      sum2 = 0.;
      sum2T = 0.;
      sum2TT = 0.;
      for (int j = 0; j < s; j++) {
        sum2 += S(j,i)*(log(XA(j,i)) - XA(j,i)/2. + 0.5);
        sum2T += S(j,i)*(1./XA(j,i)*XAT(j,i) - XAT(j,i)/2.);
        sum2TT += S(j,i)*(XATT(j,i)*(1./XA(j,i)-0.5)-XAT(j,i)*XAT(j,i)/XA(j,i)/XA(j,i));
      }
      sum1 += x[i]*sum2;
      sum1T += x[i]*sum2T;
      sum1TT += x[i]*sum2TT;
    } 
    a_assoc = sum1;    
    a_assocT = sum1T; // Derivative of a_assoc wrt temperature 
    a_assocTT = sum1TT; // Second derivative of a_assoc wrt temperature 

    // Calculating Zassoc - EQ. A10 (Chapman et al. (1990))
    for (int i = 0; i < n_comp; i++) {
      Zassoc += x[i]*mu_assoc[i];
    }  
    Zassoc -= a_assoc;

    // Calculating ZassocD
    for (int i = 0; i < n_comp; i++){
      sum1 = 0.;
      for (int k = 0; k < n_comp; k++){
        sum2 = 0.;
        for (int j = 0; j < s; j++){
          sum2 += S(j,k)*(XA_xm_rho(j,k,i)*(1./XA(j,k)-0.5)-XA_xm(j,k,i)*XAD(j,k)/XA(j,k)/XA(j,k));
        }
        sum1 += x[k]*sum2;
      }      
      ZassocD += x[i]*sum1;
    }
    ZassocD *= num_denD; // Transform derivative wrt numerical density to derivative wrt rho (mol/m³)

    // Calculating ZassocDD
    for (int i = 0; i < n_comp; i++){
      sum1 = 0.;
      for (int k = 0; k < n_comp; k++){
        sum2 = 0.;
        for (int j = 0; j < s; j++){
          sum2 += S(j,k)*(XA_xm_rho_rho(j,k,i)*(1./XA(j,k)-0.5)-2./XA(j,k)/XA(j,k)*XAD(j,k)*XA_xm_rho(j,k,i)+
            XA_xm(j,k,i)*(2.*XAD(j,k)*XAD(j,k)/XA(j,k)/XA(j,k)/XA(j,k)-
            XADD(j,k)/XA(j,k)/XA(j,k)));
        }
        sum1 += x[k]*sum2;
      }      
      ZassocDD += x[i]*sum1;
    }
    ZassocDD *= num_denD*num_denD; // Transform derivative wrt numerical density to derivative wrt rho (mol/m³)       

    // Calculating ZassocDDD
    for (int i = 0; i < n_comp; i++){
      sum1 = 0.;
      for (int k = 0; k < n_comp; k++){
        sum2 = 0.;
        for (int j = 0; j < s; j++)
        {
          sum2 += S(j,k)*(XA_xm_rho_rho_rho(j,k,i)*(1./XA(j,k)-0.5)-3./XA(j,k)/XA(j,k)*XAD(j,k)*XA_xm_rho_rho(j,k,i)+
            3.*XA_xm_rho(j,k,i)/XA(j,k)/XA(j,k)*(2.*XAD(j,k)*XAD(j,k)/XA(j,k)-XADD(j,k))+
            6.*XA_xm(j,k,i)*XAD(j,k)/XA(j,k)/XA(j,k)/XA(j,k)*(XADD(j,k)-XAD(j,k)*XAD(j,k)/XA(j,k))-
            XA_xm(j,k,i)/XA(j,k)/XA(j,k)*XADDD(j,k));
        }
        sum1 += x[k]*sum2;
      }      
      ZassocDDD += x[i]*sum1;
    }
    ZassocDDD *= num_denD*num_denD*num_denD; // Transform derivative wrt numerical density to derivative wrt rho (mol/m³) 

    // Calculating ZassocT
    for (int i = 0; i < n_comp; i++){
      sum1 = 0.;
      for (int k = 0; k < n_comp; k++){
        sum2 = 0.;
        for (int j = 0; j < s; j++){
          sum2 += S(j,k)*(XA_xm_T(j,k,i)*(1./XA(j,k)-0.5)-XA_xm(j,k,i)*XAT(j,k)/XA(j,k)/XA(j,k));
        }
        sum1 += x[k]*sum2;
      }      
      ZassocT += x[i]*sum1;
    }

    // Calculating ZassocTT
    for (int i = 0; i < n_comp; i++){
      sum1 = 0.;
      for (int k = 0; k < n_comp; k++){
        sum2 = 0.;
        for (int j = 0; j < s; j++){
          sum2 += S(j,k)*(XA_xm_T_T(j,k,i)*(1./XA(j,k)-0.5)-2./XA(j,k)/XA(j,k)*XAT(j,k)*XA_xm_T(j,k,i)+
            XA_xm(j,k,i)*(2.*XAT(j,k)*XAT(j,k)/XA(j,k)/XA(j,k)/XA(j,k)-
            XATT(j,k)/XA(j,k)/XA(j,k)));
        }
        sum1 += x[k]*sum2;
      }      
      ZassocTT += x[i]*sum1;
    }    

  } // end association term

  // Calculate Z and its derivatives
  pcsaftCalcProps.Z = Zid + Zhc + Zdisp + Zassoc; // // EQ. A24
  pcsaftCalcProps.ZD = ZidD + ZhcD + ZdispD + ZassocD; // derivative of Z wrt rho
  pcsaftCalcProps.ZDD = ZidDD + ZhcDD + ZdispDD + ZassocDD; // second derivative of Z wrt rho
  pcsaftCalcProps.ZDDD = ZidDDD + ZhcDDD + ZdispDDD + ZassocDDD; //just ZhcDDD has been implemented // third derivative of Z wrt rho
  pcsaftCalcProps.ZT = ZidT + ZhcT + ZdispT + ZassocT; // derivative of Z wrt T
  pcsaftCalcProps.ZTT = ZidTT + ZhcTT + ZdispTT + ZassocTT; // second derivative of Z wrt T

  // Calculate pressure and its derivatives
  pcsaftCalcProps.P = pcsaftCalcProps.Z*kb*T*num_den*1.e30; // Pa, EQ. A23 
  pcsaftCalcProps.PD = pcsaftCalcProps.P*(pcsaftCalcProps.ZD/pcsaftCalcProps.Z + num_denD/num_den); // derivative of P wrt rho (Pa * m³ / mol)
  pcsaftCalcProps.PDD = pcsaftCalcProps.P*(pcsaftCalcProps.ZDD/pcsaftCalcProps.Z+2.*pcsaftCalcProps.ZD/pcsaftCalcProps.Z*num_denD/num_den+num_denDD/num_den); // second derivative of P wrt rho (Pa * m⁶ / mol²)
  pcsaftCalcProps.PDDD = pcsaftCalcProps.P*(pcsaftCalcProps.ZDDD/pcsaftCalcProps.Z+3.*pcsaftCalcProps.ZDD/pcsaftCalcProps.Z*num_denD/num_den+3.*pcsaftCalcProps.ZD/pcsaftCalcProps.Z*num_denDD/num_den+num_denDDD/num_den); // third derivative of P wrt rho (Pa * m⁶ / mol²)
  pcsaftCalcProps.PT = pcsaftCalcProps.P*(pcsaftCalcProps.ZT/pcsaftCalcProps.Z + 1./T);  // derivative of P wrt T  (Pa / K)
  pcsaftCalcProps.PTT = pcsaftCalcProps.P*(pcsaftCalcProps.ZTT/pcsaftCalcProps.Z+2.*pcsaftCalcProps.ZT/pcsaftCalcProps.Z/T);  // second derivative of P wrt T (Pa / K / K)  
  
  // Begin calculating ln of fugacity coefficients and Chemical Potential  
  static real a_hs,a_hsT,a_hsTT;  
  auto logvar = log(1.-zeta[3]);
  a_hs = 1./zeta[0]*(3.*zeta[1]*zeta[2]/(1.-zeta[3]) + (zeta[2]*zeta[2]*zeta[2])/(zeta[3]*((1.-zeta[3])*(1.-zeta[3])))
      + ((zeta[2]*zeta[2]*zeta[2])/(zeta[3]*zeta[3]) - zeta[0])*logvar); // EQ. A6,  Helmholtz free energy of the hard-sphere fluid on a per-segment basis
  
  a_hsT = 1/zeta[0]*(3*(zetaT[1]*zeta[2] + zeta[1]*zetaT[2])/(1-zeta[3])
      + 3*zeta[1]*zeta[2]*zetaT[3]/((1.-zeta[3])*(1.-zeta[3])) + 3*(zeta[2]*zeta[2])*zetaT[2]/zeta[3]/((1.-zeta[3])*(1.-zeta[3]))
      + (zeta[2]*zeta[2]*zeta[2])*zetaT[3]*(3*zeta[3]-1)/(zeta[3]*zeta[3])/((1.-zeta[3])*(1.-zeta[3])*(1.-zeta[3]))
      + (3*(zeta[2]*zeta[2])*zetaT[2]*zeta[3] - 2*(zeta[2]*zeta[2]*zeta[2])*zetaT[3])/(zeta[3]*zeta[3]*zeta[3])
      * logvar + (zeta[0]-(zeta[2]*zeta[2]*zeta[2])/(zeta[3]*zeta[3]))*zetaT[3]/(1-zeta[3])); // Derivative of a_hs wrt Temperature

  a_hsTT = (-2*zeta[0]*(1 - zeta[3])*zeta[3]*zetaT[0]*(3*zeta[2]*((-1 + zeta[3])*(-1 + zeta[3]))*(zeta[3]*zeta[3]*zeta[3])*
      zetaT[1] - 3*(zeta[2]*zeta[2])*(-1 + zeta[3])*(zeta[3]*zeta[3])*zetaT[2] + 3*zeta[1]*((-1 + zeta[3])*(-1 + zeta[3]))*(zeta[3]*zeta[3]*zeta[3])*
      zetaT[2] + (zeta[2]*zeta[2]*zeta[2])*(-1 + zeta[3])*zeta[3]*zetaT[3] + 2*(zeta[2]*zeta[2]*zeta[2])*(zeta[3]*zeta[3])*zetaT[3] - 3*zeta[1]*zeta[2]*(-1 + zeta[3])*(zeta[3]*zeta[3]*zeta[3])*
      zetaT[3] + ((-1 + zeta[3])*(-1 + zeta[3]))*zeta[3]*(-(zeta[2]*zeta[2]*zeta[2]) + zeta[0]*(zeta[3]*zeta[3]))*zetaT[3] + logvar*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*
      ((zeta[3]*zeta[3]*zeta[3])*zetaT[0] - 3*(zeta[2]*zeta[2])*zeta[3]*zetaT[2] + 2*(zeta[2]*zeta[2]*zeta[2])*zetaT[3])) + ((-1 + zeta[3])*(-1 + zeta[3]))*(zeta[3]*zeta[3])*
      ((zeta[2]*zeta[2]*zeta[2])*zeta[3] - 3*zeta[1]*zeta[2]*(-1 + zeta[3])*(zeta[3]*zeta[3]) + logvar*((-1 + zeta[3])*(-1 + zeta[3]))*
      ((zeta[2]*zeta[2]*zeta[2]) - zeta[0]*(zeta[3]*zeta[3])))*(2*(zetaT[0]*zetaT[0]) - zeta[0]*zetaTT[0]) + (zeta[0]*zeta[0])*(6*((-1 + zeta[3])*(-1 + zeta[3]))*(zeta[3]*zeta[3]*zeta[3]*zeta[3])*
      zetaT[1]* (-((-1 + zeta[3])*zetaT[2]) + zeta[2]*zetaT[3]) + 2*(zeta[2]*zeta[2])*(-1 + zeta[3])*(zeta[3]*zeta[3])*zetaT[3]*
      (-3*(-1 + zeta[3])*zetaT[2] + 2*zeta[2]*zetaT[3]) - 2*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*zeta[3]*zetaT[3]*((zeta[3]*zeta[3]*zeta[3])*zetaT[0] - 
      3*(zeta[2]*zeta[2])*zeta[3]*zetaT[2] + 2*(zeta[2]*zeta[2]*zeta[2])*zetaT[3]) - 3*zeta[2]*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*(zeta[3]*zeta[3]*zeta[3]*zeta[3])*
      zetaTT[1] - ((-1 + zeta[3])*(-1 + zeta[3]))*(zeta[3]*zeta[3])*((zeta[2]*zeta[2]*zeta[2]) - zeta[0]*(zeta[3]*zeta[3]))*((zetaT[3]*zetaT[3]) - 
      (-1 + zeta[3])*zetaTT[3]) - (zeta[2]*zeta[2]*zeta[2])*((-1 + zeta[3])*(-1 + zeta[3]))*zeta[3]*(-2*(zetaT[3]*zetaT[3]) + zeta[3]*zetaTT[3]) + 
      3*zeta[1]*(1 - zeta[3])*(zeta[3]*zeta[3]*zeta[3]*zeta[3])*(-2*(-1 + zeta[3])*zetaT[2]*zetaT[3] + ((-1 + zeta[3])*(-1 + zeta[3]))*zetaTT[2] + 
      zeta[2]*(2*(zetaT[3]*zetaT[3]) - (-1 + zeta[3])*zetaTT[3])) + zeta[2]*(zeta[3]*zeta[3]*zeta[3])*(-12*zeta[2]*(-1 + zeta[3])*zetaT[2]*
      zetaT[3] + 3*((-1 + zeta[3])*(-1 + zeta[3]))*(2*(zetaT[2]*zetaT[2]) + zeta[2]*zetaTT[2]) + 2*(zeta[2]*zeta[2])*(3*(zetaT[3]*zetaT[3]) - 
      (-1 + zeta[3])*zetaTT[3])) - logvar*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*(-6*zeta[2]*(zeta[3]*zeta[3])*(zetaT[2]*zetaT[2]) + 
      (zeta[3]*zeta[3]*zeta[3]*zeta[3])*zetaTT[0] - 3*(zeta[2]*zeta[2])*zeta[3]*(-4*zetaT[2]*zetaT[3] + zeta[3]*zetaTT[2]) + (zeta[2]*zeta[2]*zeta[2])*(-6*(zetaT[3]*zetaT[3]) + 
      2*zeta[3]*zetaTT[3]))))/((zeta[0]*zeta[0]*zeta[0])*((-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3])*(-1 + zeta[3]))*(zeta[3]*zeta[3]*zeta[3]*zeta[3])); // second derivative of a_hs wrt Temperature    

  summation = 0.;
  summationT = 0.;
  summationTT = 0.;
  for (int i = 0; i < n_comp; i++) {
    aux = x[i]*(m(i)-1.);
    summation += aux*log(ghs[i*n_comp+i]);
    summationT += aux/ghs[i*n_comp+i]*ghsT[i*n_comp+i];
    summationTT += aux*(ghsTT[i*n_comp+i]/ghs[i*n_comp+i]-ghsT[i*n_comp+i]*ghsT[i*n_comp+i]/ghs[i*n_comp+i]/ghs[i*n_comp+i]);
  }

  static real a_hc,a_hcT,a_hcTT,a_disp,a_dispT,a_dispTT;
  a_hc = m_avg*a_hs - summation;  // EQ. A4, hard-chain reference contribution to Helmholtz Free Energy
  a_hcT = m_avg*a_hsT - summationT; // Derivative of a_hc wrt Temperature
  a_hcTT = m_avg*a_hsTT - summationTT; // second derivative of a_hc wrt Temperature
  a_disp = -2.*pi*num_den*I1*m2epsilonsigma3 - pi*num_den*m_avg*C1*I2*m2epsilon2sigma3; // EQ. A10, dispersion contribution to Helmholtz Free Energy
  a_dispT = -2.*pi*num_den*(I1T*m2epsilonsigma3 + I1*m2epsilonsigma3T) 
    - pi*num_den*m_avg*(C1T*I2*m2epsilon2sigma3 + C1*(I2T*m2epsilon2sigma3+I2*m2epsilon2sigma3T));
  a_dispTT = num_den*pi*(-(m_avg*(2*(I2*C1T + C1*I2T)*m2epsilon2sigma3T + m2epsilon2sigma3*(2*C1T*I2T + 
    I2*C1TT + C1*I2TT) + C1*I2*m2epsilon2sigma3TT)) - 2*(2*I1T*m2epsilonsigma3T + m2epsilonsigma3*I1TT + I1*m2epsilonsigma3TT));    

  static ArrayXr Dzeta_Dx(4);
  idx = 0;
  static real p_di;
  for (int i = 0; i < n_comp; i++) {
    p_di = 1.; 
    for (int j = 0; j < 4; j++) {
      Dzeta_Dx[j] = pi/6.*num_den*m(i)*p_di; // EQ A.34
      p_di *= d[i];
    }

    Dahs_Dx[i] = -Dzeta_Dx[0]/zeta[0]*a_hs + 1./zeta[0]*(3.*(Dzeta_Dx[1]*zeta[2]
      + zeta[1]*Dzeta_Dx[2])/(1.-zeta[3]) + 3.*zeta[1]*zeta[2]*Dzeta_Dx[3]
      /(1.-zeta[3])/(1.-zeta[3]) + 3.*zeta[2]*zeta[2]*Dzeta_Dx[2]/zeta[3]/(1.-zeta[3])/(1.-zeta[3])
      + (zeta[2]*zeta[2]*zeta[2])*Dzeta_Dx[3]*(3*zeta[3]-1.)/zeta[3]/zeta[3]/((1.-zeta[3])*(1.-zeta[3])*(1.-zeta[3]))
      + logvar*((3.*zeta[2]*zeta[2]*Dzeta_Dx[2]*zeta[3] -
      2.*(zeta[2]*zeta[2]*zeta[2])*Dzeta_Dx[3])/(zeta[3]*zeta[3]*zeta[3]) - Dzeta_Dx[0]) +
      (zeta[0]-(zeta[2]*zeta[2]*zeta[2])/zeta[3]/zeta[3])*Dzeta_Dx[3]/(1.-zeta[3])); // EQ. 36
    

    for (int j = 0; j < n_comp; j++) {
      Dghsii_Dx[idx] = Dzeta_Dx[3]/(1.-zeta[3])/(1.-zeta[3]) + (d[j]*d[j]/(d[j]+d[j]))*
        (3.*Dzeta_Dx[2]/(1.-zeta[3])/(1.-zeta[3]) + 6.*zeta[2]*Dzeta_Dx[3]/((1.-zeta[3])*(1.-zeta[3])*(1.-zeta[3])))
        + ((d[j]*d[j]/(d[j]+d[j]))*(d[j]*d[j]/(d[j]+d[j])))*(4.*zeta[2]*Dzeta_Dx[2]/((1.-zeta[3])*(1.-zeta[3])*(1.-zeta[3]))
        + 6.*zeta[2]*zeta[2]*Dzeta_Dx[3]/((1.-zeta[3])*(1.-zeta[3])*(1.-zeta[3])*(1.-zeta[3]))); // EQ. A37

      idx += 1;  
    }
    
  }

  std::fill(Dahc_Dx.begin(),Dahc_Dx.end(),0.);
  static real Dzeta3_Dx, Da_Dx, Db_Dx, DI1_Dx, DI2_Dx, Dm2es3_Dx, Dm2e2s3_Dx, DC1_Dx;
  for (int i = 0; i < n_comp; i++) {
    Dzeta3_Dx = pi/6.*num_den*m(i)*(d[i]*d[i]*d[i]);
    DI1_Dx = 0.;
    DI2_Dx = 0.;
    Dm2es3_Dx = 0.;
    Dm2e2s3_Dx = 0.;
    p_eta = 1.;
    p_eta1 = 1./eta;    
    for (int j = 0; j < 7; j++) {
      Da_Dx = m(i)/m_avg/m_avg*a1[j] + m(i)/m_avg/m_avg*(3.-4./m_avg)*a2[j]; // EQ. 44
      Db_Dx = m(i)/m_avg/m_avg*b1[j] + m(i)/m_avg/m_avg*(3.-4./m_avg)*b2[j]; // EQ. 45
      DI1_Dx += a[j]*j*Dzeta3_Dx*p_eta1 + Da_Dx*p_eta; // EQ. A42
      DI2_Dx += b[j]*j*Dzeta3_Dx*p_eta1 + Db_Dx*p_eta; // EQ. A43 
      p_eta *= eta;
      p_eta1 *= eta;
    }
    for (int j = 0; j < n_comp; j++) {
      Dm2es3_Dx += x[j]*m(j)*(epsilon_ij[i*n_comp+j]/T)*((sigma_ij[i*n_comp+j])*(sigma_ij[i*n_comp+j])*(sigma_ij[i*n_comp+j]));
      Dm2e2s3_Dx += x[j]*m(j)*((epsilon_ij[i*n_comp+j]/T)*(epsilon_ij[i*n_comp+j]/T))*((sigma_ij[i*n_comp+j])*(sigma_ij[i*n_comp+j])*(sigma_ij[i*n_comp+j]));
      Dahc_Dx[i] += x[j]*(m(j)-1.)/ghs[j*n_comp+j]*Dghsii_Dx[i*n_comp+j]; 
    }
  
    Dm2es3_Dx = 2.*m(i)*Dm2es3_Dx; // EQ. A39
    Dm2e2s3_Dx = 2.*m(i)*Dm2e2s3_Dx; // EQ. A40 
    Dahc_Dx[i] = m(i)*a_hs + m_avg*Dahs_Dx[i] - Dahc_Dx[i] - (m(i)-1.)*log(ghs[i*n_comp+i]); // EQ. A35 
    DC1_Dx = C2*Dzeta3_Dx - C1*C1*(m(i)*(8.*eta-2.*eta*eta)/((1.-eta)*(1.-eta)*(1.-eta)*(1.-eta)) -
      m(i)*(20.*eta-27.*eta*eta+12.*(eta*eta*eta)-2.*(eta*eta*eta*eta))/(((1.-eta)*(2.-eta))*((1.-eta)*(2.-eta)))); // EQ. 41

    Dadisp_Dx[i] = -2.*pi*num_den*(DI1_Dx*m2epsilonsigma3 + I1*Dm2es3_Dx) - pi*num_den
      *((m(i)*C1*I2 + m_avg*DC1_Dx*I2 + m_avg*C1*DI2_Dx)*m2epsilon2sigma3
      + m_avg*C1*I2*Dm2e2s3_Dx); // EQ. A38
  }

  std::fill(mu_hc.begin(),mu_hc.end(),0.);
  std::fill(mu_disp.begin(),mu_disp.end(),0.);
  for (int i = 0; i < n_comp; i++) {
    for (int j = 0; j < n_comp; j++) {
      mu_hc[i] += x[j]*Dahc_Dx[j];
      mu_disp[i] += x[j]*Dadisp_Dx[j];
    }
    mu_hc[i] = a_hc + Zhc + Dahc_Dx[i] - mu_hc[i]; // EQ. A33 
    mu_disp[i] = a_disp + Zdisp + Dadisp_Dx[i] - mu_disp[i]; // EQ. A33 
    mu[i] = mu_hc[i] + mu_disp[i] + mu_assoc[i]; // Chemical Potential of component i in the phase (dimensionless)
    pcsaftCalcProps.ln_fugacity_coefficients[i] = mu[i] - log(pcsaftCalcProps.Z); // EQ. A32 , ln(phi) of component i in the phase (dimensionless)
  }

  // End calculating ln of fugacity coefficients and Chemical Potential   

  // Calculating Residual Gibbs Free energy and Residual Helmholtz Free energy
  pcsaftCalcProps.Ares = a_hc + a_disp + a_assoc; // EQ. A3 , Residual Helmholtz Free energy (dimensionless)
  pcsaftCalcProps.AresT = a_hcT + a_dispT + a_assocT; // Derivative of Ares wrt Temperature (1/K)
  pcsaftCalcProps.AresTT = a_hcTT + a_dispTT + a_assocTT; // second derivative of Ares wrt Temperature (1/K/K)
  pcsaftCalcProps.Hres = -T*pcsaftCalcProps.AresT + pcsaftCalcProps.Z - 1; // Dimensionless
  pcsaftCalcProps.Hres  = pcsaftCalcProps.Hres*kb*n_av*T; // EQ. A46 , Residual Enthalpy (J/mol)
  pcsaftCalcProps.Cvres = -2*pcsaftCalcProps.AresT-T*pcsaftCalcProps.AresTT; // 1/K
  pcsaftCalcProps.Cvres = pcsaftCalcProps.Cvres*kb*n_av*T; // Isochoric heat capacity (J/mol/K)
  pcsaftCalcProps.kappaT = 1./rho/pcsaftCalcProps.PD; // Isothermal compressibility (m²/N)
  pcsaftCalcProps.alpha = 1./rho*(pcsaftCalcProps.PT/pcsaftCalcProps.PD); // Isobaric expansivity (1/K)
  pcsaftCalcProps.Cpres = pcsaftCalcProps.Cvres - kb*n_av + T*pcsaftCalcProps.alpha*pcsaftCalcProps.alpha/pcsaftCalcProps.kappaT/rho; // Isobaric heat capacity (J/mol/K) 
  pcsaftCalcProps.Gres = pcsaftCalcProps.Ares + (pcsaftCalcProps.Z - 1.) - log(pcsaftCalcProps.Z); // Dimensionless
  pcsaftCalcProps.Gres = pcsaftCalcProps.Gres*kb*n_av*T; // EQ. A50, Residual Gibbs Free energy (J/mol)  

} // end pcsaft_calculation() calculation 


};

PcsaftEOS::PcsaftEOS(const Args& args)
: pimpl(new Impl(args))
{}

PcsaftEOS::PcsaftEOS(const PcsaftEOS& other)
: pimpl(new Impl(*other.pimpl))
{}

PcsaftEOS::~PcsaftEOS()
{}

auto PcsaftEOS::operator=(PcsaftEOS other) -> PcsaftEOS&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto PcsaftEOS::compute(PcsaftEOSProps& props, real T, real P, ArrayXrConstRef x) -> void
{
    return pimpl->compute(props,T,P,x);
}

auto PcsaftEOS::pcsaft_calculation(result_pcsaft_calculation &pcsaftCalcProps, real T, real rho, ArrayXrConstRef x) -> void
{
    return pimpl->pcsaft_calculation(pcsaftCalcProps,T,rho,x);
}

auto PcsaftEOS::zbrakPD(const real rho1, const real rho2, const int n,
    ArrayXr &rhob1, ArrayXr &rhob2, int &nroot,
    real T, ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void
{
    return pimpl->zbrakPD(rho1, rho2, n, rhob1, rhob2, nroot, T, x, pcsaftCalcProps);
}    


auto PcsaftEOS::zbrakP(const real rho1, const real rho2, const int n,
    ArrayXr &rhob1, ArrayXr &rhob2, int &nroot,
    real P, real T, ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void
{
    return pimpl->zbrakP(rho1, rho2, n, rhob1, rhob2, nroot, P, T, x, pcsaftCalcProps);
}    

auto PcsaftEOS::newton_bissection(real &rho_c, real rho_a, real rho_b, real P, real T, 
    ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void
{
    return pimpl->newton_bissection(rho_c,rho_a, rho_b, P, T, x, pcsaftCalcProps);
}      

auto PcsaftEOS::newton_bissectionD(real &rho_c, real rho_a, real rho_b, real T, 
    ArrayXrConstRef x, result_pcsaft_calculation &pcsaftCalcProps) -> void   
{
    return pimpl->newton_bissectionD(rho_c,rho_a, rho_b, T, x, pcsaftCalcProps);
}  

auto PcsaftEOS::calculate_XA(int n_comp, ArrayXrConstRef x,real num_den,Vec<int> &indx) -> void
{
    return pimpl->calculate_XA(n_comp,x,num_den,indx);
} 
    
auto PcsaftEOS::resizing_vectors(result_pcsaft_calculation &pcsaftCalcProps, 
      int n_comp_New, int &n_comp_Old) -> void
{
    return pimpl->resizing_vectors(pcsaftCalcProps,n_comp_New,n_comp_Old);
} 

} // namespace Reaktoro
