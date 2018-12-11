// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "NMFSU5_info.hpp"

#include <iostream>

namespace flexiblesusy {

namespace NMFSU5_info {
   const double normalization_g5 = 1;
   const double normalization_gX = 1;

   const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities = {6, 3, 3,
      3};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_names = {"Fv",
      "Fd", "Fu", "Fe"};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names = {
      "\\nu", "d", "u", "e"};

   const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names = { "g5",
      "gX", "lam1", "lam2", "lam3", "lam3t", "lam4", "lam4t",
      "lam5", "lam5t", "lam6", "lam6t", "Re(lam7)", "Im(lam7)", "Re(lam8)",
      "Im(lam8)", "Re(eta1)", "Im(eta1)", "Re(eta2)", "Im(eta2)", "Re(eta3)",
      "Im(eta3)", "Re(Y10(0,0))", "Im(Y10(0,0))", "Re(Y10(0,1))", "Im(Y10(0,1))",
      "Re(Y10(0,2))", "Im(Y10(0,2))", "Re(Y10(1,0))", "Im(Y10(1,0))", "Re(Y10(1,1))",
      "Im(Y10(1,1))", "Re(Y10(1,2))", "Im(Y10(1,2))", "Re(Y10(2,0))", "Im(Y10(2,0))",
      "Re(Y10(2,1))", "Im(Y10(2,1))", "Re(Y10(2,2))", "Im(Y10(2,2))", "Re(Y10Pr(0,0))",
      "Im(Y10Pr(0,0))", "Re(Y10Pr(0,1))", "Im(Y10Pr(0,1))", "Re(Y10Pr(0,2))",
      "Im(Y10Pr(0,2))", "Re(Y10Pr(1,0))", "Im(Y10Pr(1,0))", "Re(Y10Pr(1,1))", "Im(Y10Pr(1,1))",
      "Re(Y10Pr(1,2))", "Im(Y10Pr(1,2))", "Re(Y10Pr(2,0))", "Im(Y10Pr(2,0))", "Re(Y10Pr(2,1))",
      "Im(Y10Pr(2,1))", "Re(Y10Pr(2,2))", "Im(Y10Pr(2,2))", "Re(Y5b(0,0))", "Im(Y5b(0,0))",
      "Re(Y5b(0,1))", "Im(Y5b(0,1))", "Re(Y5b(0,2))", "Im(Y5b(0,2))", "Re(Y5b(1,0))",
      "Im(Y5b(1,0))", "Re(Y5b(1,1))", "Im(Y5b(1,1))", "Re(Y5b(1,2))", "Im(Y5b(1,2))",
      "Re(Y5b(2,0))", "Im(Y5b(2,0))", "Re(Y5b(2,1))", "Im(Y5b(2,1))", "Re(Y5b(2,2))",
      "Im(Y5b(2,2))", "Re(Y5bPr(0,0))", "Im(Y5bPr(0,0))", "Re(Y5bPr(0,1))", "Im(Y5bPr(0,1))",
      "Re(Y5bPr(0,2))", "Im(Y5bPr(0,2))", "Re(Y5bPr(1,0))", "Im(Y5bPr(1,0))", "Re(Y5bPr(1,1))",
      "Im(Y5bPr(1,1))", "Re(Y5bPr(1,2))", "Im(Y5bPr(1,2))", "Re(Y5bPr(2,0))", "Im(Y5bPr(2,0))",
      "Re(Y5bPr(2,1))", "Im(Y5bPr(2,1))", "Re(Y5bPr(2,2))", "Im(Y5bPr(2,2))", "Re(Y1(0,0))",
      "Im(Y1(0,0))", "Re(Y1(0,1))", "Im(Y1(0,1))", "Re(Y1(0,2))", "Im(Y1(0,2))", "Re(Y1(1,0))",
      "Im(Y1(1,0))", "Re(Y1(1,1))", "Im(Y1(1,1))", "Re(Y1(1,2))", "Im(Y1(1,2))", "Re(Y1(2,0))",
      "Im(Y1(2,0))", "Re(Y1(2,1))", "Im(Y1(2,1))", "Re(Y1(2,2))", "Im(Y1(2,2))", "Re(Y1Pr(0,0))",
      "Im(Y1Pr(0,0))", "Re(Y1Pr(0,1))", "Im(Y1Pr(0,1))", "Re(Y1Pr(0,2))", "Im(Y1Pr(0,2))",
      "Re(Y1Pr(1,0))", "Im(Y1Pr(1,0))", "Re(Y1Pr(1,1))", "Im(Y1Pr(1,1))", "Re(Y1Pr(1,2))",
      "Im(Y1Pr(1,2))", "Re(Y1Pr(2,0))", "Im(Y1Pr(2,0))", "Re(Y1Pr(2,1))", "Im(Y1Pr(2,1))",
      "Re(Y1Pr(2,2))", "Im(Y1Pr(2,2))", "m10sq", "m5sq", "m5Prsq", "Re(m12sq)", "Im(m12sq)",
      "Re(mu)", "Im(mu)", "Re(muPr)", "Im(muPr)", "v", "vPr", "VG"};

   const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names = {
      "Re(Vd(0,0))", "Im(Vd(0,0))", "Re(Vd(0,1))", "Im(Vd(0,1))", "Re(Vd(0,2))",
      "Im(Vd(0,2))", "Re(Vd(1,0))", "Im(Vd(1,0))", "Re(Vd(1,1))", "Im(Vd(1,1))",
      "Re(Vd(1,2))", "Im(Vd(1,2))", "Re(Vd(2,0))", "Im(Vd(2,0))", "Re(Vd(2,1))",
      "Im(Vd(2,1))", "Re(Vd(2,2))", "Im(Vd(2,2))", "Re(Ud(0,0))", "Im(Ud(0,0))",
      "Re(Ud(0,1))", "Im(Ud(0,1))", "Re(Ud(0,2))", "Im(Ud(0,2))", "Re(Ud(1,0))",
      "Im(Ud(1,0))", "Re(Ud(1,1))", "Im(Ud(1,1))", "Re(Ud(1,2))", "Im(Ud(1,2))",
      "Re(Ud(2,0))", "Im(Ud(2,0))", "Re(Ud(2,1))", "Im(Ud(2,1))", "Re(Ud(2,2))",
      "Im(Ud(2,2))", "Re(Vu(0,0))", "Im(Vu(0,0))", "Re(Vu(0,1))", "Im(Vu(0,1))",
      "Re(Vu(0,2))", "Im(Vu(0,2))", "Re(Vu(1,0))", "Im(Vu(1,0))", "Re(Vu(1,1))",
      "Im(Vu(1,1))", "Re(Vu(1,2))", "Im(Vu(1,2))", "Re(Vu(2,0))", "Im(Vu(2,0))",
      "Re(Vu(2,1))", "Im(Vu(2,1))", "Re(Vu(2,2))", "Im(Vu(2,2))", "Re(Uu(0,0))",
      "Im(Uu(0,0))", "Re(Uu(0,1))", "Im(Uu(0,1))", "Re(Uu(0,2))", "Im(Uu(0,2))",
      "Re(Uu(1,0))", "Im(Uu(1,0))", "Re(Uu(1,1))", "Im(Uu(1,1))", "Re(Uu(1,2))",
      "Im(Uu(1,2))", "Re(Uu(2,0))", "Im(Uu(2,0))", "Re(Uu(2,1))", "Im(Uu(2,1))",
      "Re(Uu(2,2))", "Im(Uu(2,2))", "Re(Ve(0,0))", "Im(Ve(0,0))", "Re(Ve(0,1))",
      "Im(Ve(0,1))", "Re(Ve(0,2))", "Im(Ve(0,2))", "Re(Ve(1,0))", "Im(Ve(1,0))",
      "Re(Ve(1,1))", "Im(Ve(1,1))", "Re(Ve(1,2))", "Im(Ve(1,2))", "Re(Ve(2,0))",
      "Im(Ve(2,0))", "Re(Ve(2,1))", "Im(Ve(2,1))", "Re(Ve(2,2))", "Im(Ve(2,2))",
      "Re(Ue(0,0))", "Im(Ue(0,0))", "Re(Ue(0,1))", "Im(Ue(0,1))", "Re(Ue(0,2))",
      "Im(Ue(0,2))", "Re(Ue(1,0))", "Im(Ue(1,0))", "Re(Ue(1,1))", "Im(Ue(1,1))",
      "Re(Ue(1,2))", "Im(Ue(1,2))", "Re(Ue(2,0))", "Im(Ue(2,0))", "Re(Ue(2,1))",
      "Im(Ue(2,1))", "Re(Ue(2,2))", "Im(Ue(2,2))", "Re(UV(0,0))", "Im(UV(0,0))",
      "Re(UV(0,1))", "Im(UV(0,1))", "Re(UV(0,2))", "Im(UV(0,2))", "Re(UV(0,3))",
      "Im(UV(0,3))", "Re(UV(0,4))", "Im(UV(0,4))", "Re(UV(0,5))", "Im(UV(0,5))",
      "Re(UV(1,0))", "Im(UV(1,0))", "Re(UV(1,1))", "Im(UV(1,1))", "Re(UV(1,2))",
      "Im(UV(1,2))", "Re(UV(1,3))", "Im(UV(1,3))", "Re(UV(1,4))", "Im(UV(1,4))",
      "Re(UV(1,5))", "Im(UV(1,5))", "Re(UV(2,0))", "Im(UV(2,0))", "Re(UV(2,1))",
      "Im(UV(2,1))", "Re(UV(2,2))", "Im(UV(2,2))", "Re(UV(2,3))", "Im(UV(2,3))",
      "Re(UV(2,4))", "Im(UV(2,4))", "Re(UV(2,5))", "Im(UV(2,5))", "Re(UV(3,0))",
      "Im(UV(3,0))", "Re(UV(3,1))", "Im(UV(3,1))", "Re(UV(3,2))", "Im(UV(3,2))",
      "Re(UV(3,3))", "Im(UV(3,3))", "Re(UV(3,4))", "Im(UV(3,4))", "Re(UV(3,5))",
      "Im(UV(3,5))", "Re(UV(4,0))", "Im(UV(4,0))", "Re(UV(4,1))", "Im(UV(4,1))",
      "Re(UV(4,2))", "Im(UV(4,2))", "Re(UV(4,3))", "Im(UV(4,3))", "Re(UV(4,4))",
      "Im(UV(4,4))", "Re(UV(4,5))", "Im(UV(4,5))", "Re(UV(5,0))", "Im(UV(5,0))",
      "Re(UV(5,1))", "Im(UV(5,1))", "Re(UV(5,2))", "Im(UV(5,2))", "Re(UV(5,3))",
      "Im(UV(5,3))", "Re(UV(5,4))", "Im(UV(5,4))", "Re(UV(5,5))", "Im(UV(5,5))"
   };

   const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names
       = {"lam1IN", "lam2IN", "lam3IN", "lam3tIN", "lam4IN",
          "lam4tIN", "lam5IN", "lam5tIN", "lam6IN", "lam6tIN",
          "Re(lam7IN)", "Im(lam7IN)", "Re(lam8IN)", "Im(lam8IN)",
          "Re(eta1IN)", "Im(eta1IN)", "Re(eta2IN)", "Im(eta2IN)",
          "Re(eta3IN)", "Im(eta3IN)"};

   const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names
       = {};

   const std::string model_name = "NMFSU5";

void print(std::ostream& ostr)
{
   ostr
      << "Model information\n"
      << "=================\n"
      << "Model name:                  " << model_name << '\n'
      << "Is a low-energy model:       "
      << (is_low_energy_model ? "yes" : "no") << '\n'
      << "Is a supersymmetric model:   "
      << (is_supersymmetric_model ? "yes" : "no") << '\n'
      << "Is a FlexibleEFTHiggs model: "
      << (is_FlexibleEFTHiggs ? "yes" : "no") << '\n'
      << "Number of multiplets:        " << NUMBER_OF_PARTICLES << '\n'
      << "Number of parameters:        " << NUMBER_OF_PARAMETERS << '\n'
      ;

   ostr << "\n"
      "Multiplets:                  ";
   for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
      ostr << particle_names[i]
           << '[' << particle_multiplicities[i] << ']';
      if (i + 1 < NUMBER_OF_PARTICLES)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Parameters:                  ";
   for (int i = 0; i < NUMBER_OF_PARAMETERS; i++) {
      ostr << parameter_names[i];
      if (i + 1 < NUMBER_OF_PARAMETERS)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Input parameters:            ";
   for (int i = 0; i < NUMBER_OF_INPUT_PARAMETERS; i++) {
      ostr << input_parameter_names[i];
      if (i + 1 < NUMBER_OF_INPUT_PARAMETERS)
         ostr << ", ";
   }

   ostr << '\n';
}

} // namespace NMFSU5_info

} // namespace flexiblesusy

