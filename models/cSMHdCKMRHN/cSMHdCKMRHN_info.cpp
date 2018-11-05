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

// File generated at Mon 5 Nov 2018 12:48:53

#include "cSMHdCKMRHN_info.hpp"

#include <iostream>

namespace flexiblesusy {

namespace cSMHdCKMRHN_info {
   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;

   const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities = {1, 1, 1,
      1, 3, 3, 3, 6, 1, 1, 1};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_names = {"VG", "Hm"
      , "Ah", "hh", "Fd", "Fu", "Fe", "Fv", "VWp", "VP", "VZ"};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names = {
      "g", "H^-", "A^0", "h", "d", "u", "e", "\\nu", "W^+", "\\gamma", "Z"};

   const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names = {"g1",
      "g2", "g3", "Lambdax", "Re(Yd(0,0))", "Im(Yd(0,0))", "Re(Yd(0,1))",
      "Im(Yd(0,1))", "Re(Yd(0,2))", "Im(Yd(0,2))", "Re(Yd(1,0))", "Im(Yd(1,0))",
      "Re(Yd(1,1))", "Im(Yd(1,1))", "Re(Yd(1,2))", "Im(Yd(1,2))", "Re(Yd(2,0))",
      "Im(Yd(2,0))", "Re(Yd(2,1))", "Im(Yd(2,1))", "Re(Yd(2,2))", "Im(Yd(2,2))",
      "Re(Ye(0,0))", "Im(Ye(0,0))", "Re(Ye(0,1))", "Im(Ye(0,1))", "Re(Ye(0,2))",
      "Im(Ye(0,2))", "Re(Ye(1,0))", "Im(Ye(1,0))", "Re(Ye(1,1))", "Im(Ye(1,1))",
      "Re(Ye(1,2))", "Im(Ye(1,2))", "Re(Ye(2,0))", "Im(Ye(2,0))", "Re(Ye(2,1))",
      "Im(Ye(2,1))", "Re(Ye(2,2))", "Im(Ye(2,2))", "Re(Yv(0,0))", "Im(Yv(0,0))",
      "Re(Yv(0,1))", "Im(Yv(0,1))", "Re(Yv(0,2))", "Im(Yv(0,2))", "Re(Yv(1,0))",
      "Im(Yv(1,0))", "Re(Yv(1,1))", "Im(Yv(1,1))", "Re(Yv(1,2))", "Im(Yv(1,2))",
      "Re(Yv(2,0))", "Im(Yv(2,0))", "Re(Yv(2,1))", "Im(Yv(2,1))", "Re(Yv(2,2))",
      "Im(Yv(2,2))", "Re(Yu(0,0))", "Im(Yu(0,0))", "Re(Yu(0,1))", "Im(Yu(0,1))",
      "Re(Yu(0,2))", "Im(Yu(0,2))", "Re(Yu(1,0))", "Im(Yu(1,0))", "Re(Yu(1,1))",
      "Im(Yu(1,1))", "Re(Yu(1,2))", "Im(Yu(1,2))", "Re(Yu(2,0))", "Im(Yu(2,0))",
      "Re(Yu(2,1))", "Im(Yu(2,1))", "Re(Yu(2,2))", "Im(Yu(2,2))", "Re(Mv(0,0))",
      "Im(Mv(0,0))", "Re(Mv(0,1))", "Im(Mv(0,1))", "Re(Mv(0,2))", "Im(Mv(0,2))",
      "Re(Mv(1,0))", "Im(Mv(1,0))", "Re(Mv(1,1))", "Im(Mv(1,1))", "Re(Mv(1,2))",
      "Im(Mv(1,2))", "Re(Mv(2,0))", "Im(Mv(2,0))", "Re(Mv(2,1))", "Im(Mv(2,1))",
      "Re(Mv(2,2))", "Im(Mv(2,2))", "mu2", "v"};

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
      "Im(UV(5,3))", "Re(UV(5,4))", "Im(UV(5,4))", "Re(UV(5,5))", "Im(UV(5,5))",
      "ZZ(0,0)", "ZZ(0,1)", "ZZ(1,0)", "ZZ(1,1)"};

   const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names
       = {"LambdaIN", "Qin", "QEWSB", "MvInput(0,0)", "MvInput(0,1)",
      "MvInput(0,2)", "MvInput(1,0)", "MvInput(1,1)", "MvInput(1,2)",
      "MvInput(2,0)", "MvInput(2,1)", "MvInput(2,2)"};

   const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names
       = {};

   const std::string model_name = "cSMHdCKMRHN";

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

} // namespace cSMHdCKMRHN_info

} // namespace flexiblesusy

