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

#ifndef cSMHdCKMRHN_INFO_H
#define cSMHdCKMRHN_INFO_H

#include "problems.hpp"

#include <array>
#include <iosfwd>
#include <string>

namespace flexiblesusy {

namespace cSMHdCKMRHN_info {
   enum Particles : int { VG, Hm, Ah, hh, Fd, Fu, Fe, Fv, VWp, VP, VZ,
      NUMBER_OF_PARTICLES };

   enum Masses : int { MVG, MHm, MAh, Mhh, MFd_1, MFd_2, MFd_3, MFu_1, MFu_2,
      MFu_3, MFe_1, MFe_2, MFe_3, MFv_1, MFv_2, MFv_3, MFv_4, MFv_5, MFv_6, MVWp,
      MVP, MVZ, NUMBER_OF_MASSES };

   enum Parameters : int { g1, g2, g3, Lambdax, ReYd0_0, ImYd0_0, ReYd0_1, ImYd0_1
      , ReYd0_2, ImYd0_2, ReYd1_0, ImYd1_0, ReYd1_1, ImYd1_1, ReYd1_2, ImYd1_2,
      ReYd2_0, ImYd2_0, ReYd2_1, ImYd2_1, ReYd2_2, ImYd2_2, ReYe0_0, ImYe0_0,
      ReYe0_1, ImYe0_1, ReYe0_2, ImYe0_2, ReYe1_0, ImYe1_0, ReYe1_1, ImYe1_1,
      ReYe1_2, ImYe1_2, ReYe2_0, ImYe2_0, ReYe2_1, ImYe2_1, ReYe2_2, ImYe2_2,
      ReYv0_0, ImYv0_0, ReYv0_1, ImYv0_1, ReYv0_2, ImYv0_2, ReYv1_0, ImYv1_0,
      ReYv1_1, ImYv1_1, ReYv1_2, ImYv1_2, ReYv2_0, ImYv2_0, ReYv2_1, ImYv2_1,
      ReYv2_2, ImYv2_2, ReYu0_0, ImYu0_0, ReYu0_1, ImYu0_1, ReYu0_2, ImYu0_2,
      ReYu1_0, ImYu1_0, ReYu1_1, ImYu1_1, ReYu1_2, ImYu1_2, ReYu2_0, ImYu2_0,
      ReYu2_1, ImYu2_1, ReYu2_2, ImYu2_2, ReMv0_0, ImMv0_0, ReMv0_1, ImMv0_1,
      ReMv0_2, ImMv0_2, ReMv1_0, ImMv1_0, ReMv1_1, ImMv1_1, ReMv1_2, ImMv1_2,
      ReMv2_0, ImMv2_0, ReMv2_1, ImMv2_1, ReMv2_2, ImMv2_2, mu2, v,
      NUMBER_OF_PARAMETERS };

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
      ReUV5_3, ImUV5_3, ReUV5_4, ImUV5_4, ReUV5_5, ImUV5_5, ZZ0_0, ZZ0_1, ZZ1_0,
      ZZ1_1, NUMBER_OF_MIXINGS };

   enum Input_parameters : int { LambdaIN, Qin, QEWSB, MvInput0_0, MvInput0_1,
      MvInput0_2, MvInput1_0, MvInput1_1, MvInput1_2, MvInput2_0, MvInput2_1,
      MvInput2_2, NUMBER_OF_INPUT_PARAMETERS };

   enum Extra_parameters : int { NUMBER_OF_EXTRA_PARAMETERS };

   extern const double normalization_g1;
   extern const double normalization_g2;
   extern const double normalization_g3;

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

   class cSMHdCKMRHN_particle_names : public Names {
   public:
      virtual ~cSMHdCKMRHN_particle_names() = default;
      virtual const std::string& get(int index) const override {
         return particle_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARTICLES;
      }
   };

   class cSMHdCKMRHN_parameter_names : public Names {
   public:
      virtual ~cSMHdCKMRHN_parameter_names() = default;
      virtual const std::string& get(int index) const override {
         return parameter_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARAMETERS;
      }
   };

   const cSMHdCKMRHN_particle_names  particle_names_getter{};
   const cSMHdCKMRHN_parameter_names parameter_names_getter{};

} // namespace cSMHdCKMRHN_info

} // namespace flexiblesusy

#endif
