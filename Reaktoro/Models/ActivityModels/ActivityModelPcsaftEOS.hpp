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

#pragma once

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Models/ActivityModels/Support/PcsaftEOS.hpp>
#include "Reaktoro/Common/NamingUtils.hpp"

namespace Reaktoro {

/// Return the activity model for fluid phases based on the Pcsaft equation of state.
auto ActivityModelPcsaft(String file, MatrixXr kij) -> ActivityModelGenerator;

auto correctName(String name) -> String;

auto split_string (const String &s) -> std::vector<String>;

auto get_param(std::ifstream &file) -> double;

auto read_pcsaft_database(const SpeciesList& species, String file,
      Map<String,PcsaftEOSModel> types,Vec<String> &pcsaft_type, ArrayXr &m, ArrayXr &sigma, 
      ArrayXr &epsilon, ArrayXr &kappa, ArrayXr &eassoc, Vec<unsigned> &n_donor_sites, 
      Vec<unsigned> &n_acceptor_sites) -> void;

} // namespace Reaktoro
