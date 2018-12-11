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

#ifndef NMFSU5_INFO_H
#define NMFSU5_INFO_H

#include "problems.hpp"

#include <array>
#include <iosfwd>
#include <string>

namespace flexiblesusy {

namespace NMFSU5_info {
   enum Particles : int { Fv, Fd, Fu, Fe,
      NUMBER_OF_PARTICLES };

   enum Masses : int { MFv_1, MFv_2, MFv_3, MFv_4, MFv_5, MFv_6, MFd_1, MFd_2,
      MFd_3, MFu_1, MFu_2, MFu_3, MFe_1, MFe_2, MFe_3,
      NUMBER_OF_MASSES };

   enum Parameters : int { g5, gX, lam1, lam2, lam3, lam3t, lam4, lam4t,
      lam5, lam5t, lam6, lam6t, Relam7, Imlam7, Relam8, Imlam8, Reeta1, Imeta1,
      Reeta2, Imeta2, Reeta3, Imeta3, ReY100_0, ImY100_0, ReY100_1, ImY100_1,
      ReY100_2, ImY100_2, ReY101_0, ImY101_0, ReY101_1, ImY101_1,
      ReY101_2, ImY101_2, ReY102_0, ImY102_0, ReY102_1, ImY102_1,
      ReY102_2, ImY102_2, ReY10Pr0_0, ImY10Pr0_0, ReY10Pr0_1, ImY10Pr0_1,
      ReY10Pr0_2, ImY10Pr0_2, ReY10Pr1_0, ImY10Pr1_0, ReY10Pr1_1, ImY10Pr1_1,
      ReY10Pr1_2, ImY10Pr1_2, ReY10Pr2_0, ImY10Pr2_0, ReY10Pr2_1, ImY10Pr2_1,
      ReY10Pr2_2, ImY10Pr2_2, ReY5b0_0, ImY5b0_0, ReY5b0_1, ImY5b0_1,
      ReY5b0_2, ImY5b0_2, ReY5b1_0, ImY5b1_0, ReY5b1_1, ImY5b1_1,
      ReY5b1_2, ImY5b1_2, ReY5b2_0, ImY5b2_0, ReY5b2_1, ImY5b2_1,
      ReY5b2_2, ImY5b2_2, ReY5bPr0_0, ImY5bPr0_0, ReY5bPr0_1, ImY5bPr0_1,
      ReY5bPr0_2, ImY5bPr0_2, ReY5bPr1_0, ImY5bPr1_0, ReY5bPr1_1, ImY5bPr1_1,
      ReY5bPr1_2, ImY5bPr1_2, ReY5bPr2_0, ImY5bPr2_0, ReY5bPr2_1, ImY5bPr2_1,
      ReY5bPr2_2, ImY5bPr2_2, ReY10_0, ImY10_0, ReY10_1, ImY10_1,
      ReY10_2, ImY10_2, ReY11_0, ImY11_0, ReY11_1, ImY11_1,
      ReY11_2, ImY11_2, ReY12_0, ImY12_0, ReY12_1, ImY12_1,
      ReY12_2, ImY12_2, ReY1Pr0_0, ImY1Pr0_0, ReY1Pr0_1, ImY1Pr0_1,
      ReY1Pr0_2, ImY1Pr0_2, ReY1Pr1_0, ImY1Pr1_0, ReY1Pr1_1, ImY1Pr1_1,
      ReY1Pr1_2, ImY1Pr1_2, ReY1Pr2_0, ImY1Pr2_0, ReY1Pr2_1, ImY1Pr2_1,
      ReY1Pr2_2, ImY1Pr2_2, m10sq, m5sq, m5Prsq, Rem12sq, Imm12sq, Remu,
      Immu, RemuPr, ImmuPr, v, vPr, VG, NUMBER_OF_PARAMETERS };

   enum Mixings : int { ReVd0_0, ImVd0_0, ReVd0_1, ImVd0_1, ReVd0_2, ImVd0_2,
      ReVd1_0, ImVd1_0, ReVd1_1, ImVd1_1, ReVd1_2, ImVd1_2, ReVd2_0, ImVd2_0,
      ReVd2_1, ImVd2_1, ReVd2_2, ImVd2_2, ReUd0_0, ImUd0_0, ReUd0_1, ImUd0_1,
      ReUd0_2, ImUd0_2, ReUd1_0, ImUd1_0, ReUd1_1, ImUd1_1, ReUd1_2, ImUd1_2,
      ReUd2_0, ImUd2_0, ReUd2_1, ImUd2_1, ReUd2_2, ImUd2_2, ReVu0_0, ImVu0_0,
      ReVu0_1, ImVu0_1, ReVu0_2, ImVu0_2, ReVu1_0, ImVu1_0, ReVu1_1, ImVu1_1,
      ReVu1_2, ImVu1_2, ReVu2_0, ImVu2_0, ReVu2_1, ImVu2_1, ReVu2_2, ImVu2_2,
      ReUu0_0, ImUu0_0, ReUu0_1, ImUu0_1, ReUu0_2, ImUu0_2, ReUu1_0, ImUu1_0,
      ReUu1_1, ImUu1_1, ReUu1_2, ImUu1_2, ReUu2_0, ImUu2_0, ReUu2_1, ImUu2_1,
      ReUu2_2, ImUu2_2, ReVe0_0, ImVe0_0, ReVe0_1, ImVe0_1, ReVe0_2, ImVe0_2,
      ReVe1_0, ImVe1_0, ReVe1_1, ImVe1_1, ReVe1_2, ImVe1_2, ReVe2_0, ImVe2_0,
      ReVe2_1, ImVe2_1, ReVe2_2, ImVe2_2, ReUe0_0, ImUe0_0, ReUe0_1, ImUe0_1,
      ReUe0_2, ImUe0_2, ReUe1_0, ImUe1_0, ReUe1_1, ImUe1_1, ReUe1_2, ImUe1_2,
      ReUe2_0, ImUe2_0, ReUe2_1, ImUe2_1, ReUe2_2, ImUe2_2, ReUV0_0, ImUV0_0,
      ReUV0_1, ImUV0_1, ReUV0_2, ImUV0_2, ReUV0_3, ImUV0_3, ReUV0_4, ImUV0_4,
      ReUV0_5, ImUV0_5, ReUV1_0, ImUV1_0, ReUV1_1, ImUV1_1, ReUV1_2, ImUV1_2,
      ReUV1_3, ImUV1_3, ReUV1_4, ImUV1_4, ReUV1_5, ImUV1_5, ReUV2_0, ImUV2_0,
      ReUV2_1, ImUV2_1, ReUV2_2, ImUV2_2, ReUV2_3, ImUV2_3, ReUV2_4, ImUV2_4,
      ReUV2_5, ImUV2_5, ReUV3_0, ImUV3_0, ReUV3_1, ImUV3_1, ReUV3_2, ImUV3_2,
      ReUV3_3, ImUV3_3, ReUV3_4, ImUV3_4, ReUV3_5, ImUV3_5, ReUV4_0, ImUV4_0,
      ReUV4_1, ImUV4_1, ReUV4_2, ImUV4_2, ReUV4_3, ImUV4_3, ReUV4_4, ImUV4_4,
      ReUV4_5, ImUV4_5, ReUV5_0, ImUV5_0, ReUV5_1, ImUV5_1, ReUV5_2, ImUV5_2,
      ReUV5_3, ImUV5_3, ReUV5_4, ImUV5_4, ReUV5_5, ImUV5_5,
      NUMBER_OF_MIXINGS };

   enum Input_parameters : int { lam1IN, lam2IN, lam3IN, lam3tIN, lam4IN,
      lam4tIN, lam5IN, lam5tIN, lam6IN, lam6tIN, Relam7IN, Imlam7IN, Relam8IN,
      Imlam8IN, Reeta1IN, Imeta1IN, Reeta2IN, Imeta2IN, Reeta3IN, Imeta3IN,
      NUMBER_OF_INPUT_PARAMETERS
      };

   enum Extra_parameters : int { NUMBER_OF_EXTRA_PARAMETERS };

   extern const double normalization_g5;
   extern const double normalization_gX;

   extern const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_names;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names;
   extern const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names;
   extern const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names;
   extern const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names;
   extern const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names;
   extern const std::string model_name;
   constexpr bool is_low_energy_model = false;
   constexpr bool is_supersymmetric_model = false;
   constexpr bool is_FlexibleEFTHiggs = false;

   void print(std::ostream&);

   class NMFSU5_particle_names : public Names {
   public:
      virtual ~NMFSU5_particle_names() = default;
      virtual const std::string& get(int index) const override {
         return particle_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARTICLES;
      }
   };

   class NMFSU5_parameter_names : public Names {
   public:
      virtual ~NMFSU5_parameter_names() = default;
      virtual const std::string& get(int index) const override {
         return parameter_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARAMETERS;
      }
   };

   const NMFSU5_particle_names  particle_names_getter{};
   const NMFSU5_parameter_names parameter_names_getter{};

} // namespace NMFSU5_info

} // namespace flexiblesusy

#endif
