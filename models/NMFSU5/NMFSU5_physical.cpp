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

#include "NMFSU5_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

void NMFSU5_physical::clear()
{
   MFv = Eigen::Matrix<double,6,1>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   UV = Eigen::Matrix<std::complex<double>,6,6>::Zero();
}

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void NMFSU5_physical::convert_to_hk()
{

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void NMFSU5_physical::convert_to_slha()
{

}

Eigen::ArrayXd NMFSU5_physical::get() const
{
   Eigen::ArrayXd pars(get_masses());

   pars.conservativeResize(195);

   pars(15) = Re(Vd(0,0));
   pars(16) = Im(Vd(0,0));
   pars(17) = Re(Vd(0,1));
   pars(18) = Im(Vd(0,1));
   pars(19) = Re(Vd(0,2));
   pars(20) = Im(Vd(0,2));
   pars(21) = Re(Vd(1,0));
   pars(22) = Im(Vd(1,0));
   pars(23) = Re(Vd(1,1));
   pars(24) = Im(Vd(1,1));
   pars(25) = Re(Vd(1,2));
   pars(26) = Im(Vd(1,2));
   pars(27) = Re(Vd(2,0));
   pars(28) = Im(Vd(2,0));
   pars(29) = Re(Vd(2,1));
   pars(30) = Im(Vd(2,1));
   pars(31) = Re(Vd(2,2));
   pars(32) = Im(Vd(2,2));
   pars(33) = Re(Ud(0,0));
   pars(34) = Im(Ud(0,0));
   pars(35) = Re(Ud(0,1));
   pars(36) = Im(Ud(0,1));
   pars(37) = Re(Ud(0,2));
   pars(38) = Im(Ud(0,2));
   pars(39) = Re(Ud(1,0));
   pars(40) = Im(Ud(1,0));
   pars(41) = Re(Ud(1,1));
   pars(42) = Im(Ud(1,1));
   pars(43) = Re(Ud(1,2));
   pars(44) = Im(Ud(1,2));
   pars(45) = Re(Ud(2,0));
   pars(46) = Im(Ud(2,0));
   pars(47) = Re(Ud(2,1));
   pars(48) = Im(Ud(2,1));
   pars(49) = Re(Ud(2,2));
   pars(50) = Im(Ud(2,2));
   pars(51) = Re(Vu(0,0));
   pars(52) = Im(Vu(0,0));
   pars(53) = Re(Vu(0,1));
   pars(54) = Im(Vu(0,1));
   pars(55) = Re(Vu(0,2));
   pars(56) = Im(Vu(0,2));
   pars(57) = Re(Vu(1,0));
   pars(58) = Im(Vu(1,0));
   pars(59) = Re(Vu(1,1));
   pars(60) = Im(Vu(1,1));
   pars(61) = Re(Vu(1,2));
   pars(62) = Im(Vu(1,2));
   pars(63) = Re(Vu(2,0));
   pars(64) = Im(Vu(2,0));
   pars(65) = Re(Vu(2,1));
   pars(66) = Im(Vu(2,1));
   pars(67) = Re(Vu(2,2));
   pars(68) = Im(Vu(2,2));
   pars(69) = Re(Uu(0,0));
   pars(70) = Im(Uu(0,0));
   pars(71) = Re(Uu(0,1));
   pars(72) = Im(Uu(0,1));
   pars(73) = Re(Uu(0,2));
   pars(74) = Im(Uu(0,2));
   pars(75) = Re(Uu(1,0));
   pars(76) = Im(Uu(1,0));
   pars(77) = Re(Uu(1,1));
   pars(78) = Im(Uu(1,1));
   pars(79) = Re(Uu(1,2));
   pars(80) = Im(Uu(1,2));
   pars(81) = Re(Uu(2,0));
   pars(82) = Im(Uu(2,0));
   pars(83) = Re(Uu(2,1));
   pars(84) = Im(Uu(2,1));
   pars(85) = Re(Uu(2,2));
   pars(86) = Im(Uu(2,2));
   pars(87) = Re(Ve(0,0));
   pars(88) = Im(Ve(0,0));
   pars(89) = Re(Ve(0,1));
   pars(90) = Im(Ve(0,1));
   pars(91) = Re(Ve(0,2));
   pars(92) = Im(Ve(0,2));
   pars(93) = Re(Ve(1,0));
   pars(94) = Im(Ve(1,0));
   pars(95) = Re(Ve(1,1));
   pars(96) = Im(Ve(1,1));
   pars(97) = Re(Ve(1,2));
   pars(98) = Im(Ve(1,2));
   pars(99) = Re(Ve(2,0));
   pars(100) = Im(Ve(2,0));
   pars(101) = Re(Ve(2,1));
   pars(102) = Im(Ve(2,1));
   pars(103) = Re(Ve(2,2));
   pars(104) = Im(Ve(2,2));
   pars(105) = Re(Ue(0,0));
   pars(106) = Im(Ue(0,0));
   pars(107) = Re(Ue(0,1));
   pars(108) = Im(Ue(0,1));
   pars(109) = Re(Ue(0,2));
   pars(110) = Im(Ue(0,2));
   pars(111) = Re(Ue(1,0));
   pars(112) = Im(Ue(1,0));
   pars(113) = Re(Ue(1,1));
   pars(114) = Im(Ue(1,1));
   pars(115) = Re(Ue(1,2));
   pars(116) = Im(Ue(1,2));
   pars(117) = Re(Ue(2,0));
   pars(118) = Im(Ue(2,0));
   pars(119) = Re(Ue(2,1));
   pars(120) = Im(Ue(2,1));
   pars(121) = Re(Ue(2,2));
   pars(122) = Im(Ue(2,2));
   pars(123) = Re(UV(0,0));
   pars(124) = Im(UV(0,0));
   pars(125) = Re(UV(0,1));
   pars(126) = Im(UV(0,1));
   pars(127) = Re(UV(0,2));
   pars(128) = Im(UV(0,2));
   pars(129) = Re(UV(0,3));
   pars(130) = Im(UV(0,3));
   pars(131) = Re(UV(0,4));
   pars(132) = Im(UV(0,4));
   pars(133) = Re(UV(0,5));
   pars(134) = Im(UV(0,5));
   pars(135) = Re(UV(1,0));
   pars(136) = Im(UV(1,0));
   pars(137) = Re(UV(1,1));
   pars(138) = Im(UV(1,1));
   pars(139) = Re(UV(1,2));
   pars(140) = Im(UV(1,2));
   pars(141) = Re(UV(1,3));
   pars(142) = Im(UV(1,3));
   pars(143) = Re(UV(1,4));
   pars(144) = Im(UV(1,4));
   pars(145) = Re(UV(1,5));
   pars(146) = Im(UV(1,5));
   pars(147) = Re(UV(2,0));
   pars(148) = Im(UV(2,0));
   pars(149) = Re(UV(2,1));
   pars(150) = Im(UV(2,1));
   pars(151) = Re(UV(2,2));
   pars(152) = Im(UV(2,2));
   pars(153) = Re(UV(2,3));
   pars(154) = Im(UV(2,3));
   pars(155) = Re(UV(2,4));
   pars(156) = Im(UV(2,4));
   pars(157) = Re(UV(2,5));
   pars(158) = Im(UV(2,5));
   pars(159) = Re(UV(3,0));
   pars(160) = Im(UV(3,0));
   pars(161) = Re(UV(3,1));
   pars(162) = Im(UV(3,1));
   pars(163) = Re(UV(3,2));
   pars(164) = Im(UV(3,2));
   pars(165) = Re(UV(3,3));
   pars(166) = Im(UV(3,3));
   pars(167) = Re(UV(3,4));
   pars(168) = Im(UV(3,4));
   pars(169) = Re(UV(3,5));
   pars(170) = Im(UV(3,5));
   pars(171) = Re(UV(4,0));
   pars(172) = Im(UV(4,0));
   pars(173) = Re(UV(4,1));
   pars(174) = Im(UV(4,1));
   pars(175) = Re(UV(4,2));
   pars(176) = Im(UV(4,2));
   pars(177) = Re(UV(4,3));
   pars(178) = Im(UV(4,3));
   pars(179) = Re(UV(4,4));
   pars(180) = Im(UV(4,4));
   pars(181) = Re(UV(4,5));
   pars(182) = Im(UV(4,5));
   pars(183) = Re(UV(5,0));
   pars(184) = Im(UV(5,0));
   pars(185) = Re(UV(5,1));
   pars(186) = Im(UV(5,1));
   pars(187) = Re(UV(5,2));
   pars(188) = Im(UV(5,2));
   pars(189) = Re(UV(5,3));
   pars(190) = Im(UV(5,3));
   pars(191) = Re(UV(5,4));
   pars(192) = Im(UV(5,4));
   pars(193) = Re(UV(5,5));
   pars(194) = Im(UV(5,5));

   return pars;
}

void NMFSU5_physical::set(const Eigen::ArrayXd& pars)
{
   set_masses(pars);

   Vd(0,0) = std::complex<double>(pars(15), pars(16));
   Vd(0,1) = std::complex<double>(pars(17), pars(18));
   Vd(0,2) = std::complex<double>(pars(19), pars(20));
   Vd(1,0) = std::complex<double>(pars(21), pars(22));
   Vd(1,1) = std::complex<double>(pars(23), pars(24));
   Vd(1,2) = std::complex<double>(pars(25), pars(26));
   Vd(2,0) = std::complex<double>(pars(27), pars(28));
   Vd(2,1) = std::complex<double>(pars(29), pars(30));
   Vd(2,2) = std::complex<double>(pars(31), pars(32));
   Ud(0,0) = std::complex<double>(pars(33), pars(34));
   Ud(0,1) = std::complex<double>(pars(35), pars(36));
   Ud(0,2) = std::complex<double>(pars(37), pars(38));
   Ud(1,0) = std::complex<double>(pars(39), pars(40));
   Ud(1,1) = std::complex<double>(pars(41), pars(42));
   Ud(1,2) = std::complex<double>(pars(43), pars(44));
   Ud(2,0) = std::complex<double>(pars(45), pars(46));
   Ud(2,1) = std::complex<double>(pars(47), pars(48));
   Ud(2,2) = std::complex<double>(pars(49), pars(50));
   Vu(0,0) = std::complex<double>(pars(51), pars(52));
   Vu(0,1) = std::complex<double>(pars(53), pars(54));
   Vu(0,2) = std::complex<double>(pars(55), pars(56));
   Vu(1,0) = std::complex<double>(pars(57), pars(58));
   Vu(1,1) = std::complex<double>(pars(59), pars(60));
   Vu(1,2) = std::complex<double>(pars(61), pars(62));
   Vu(2,0) = std::complex<double>(pars(63), pars(64));
   Vu(2,1) = std::complex<double>(pars(65), pars(66));
   Vu(2,2) = std::complex<double>(pars(67), pars(68));
   Uu(0,0) = std::complex<double>(pars(69), pars(70));
   Uu(0,1) = std::complex<double>(pars(71), pars(72));
   Uu(0,2) = std::complex<double>(pars(73), pars(74));
   Uu(1,0) = std::complex<double>(pars(75), pars(76));
   Uu(1,1) = std::complex<double>(pars(77), pars(78));
   Uu(1,2) = std::complex<double>(pars(79), pars(80));
   Uu(2,0) = std::complex<double>(pars(81), pars(82));
   Uu(2,1) = std::complex<double>(pars(83), pars(84));
   Uu(2,2) = std::complex<double>(pars(85), pars(86));
   Ve(0,0) = std::complex<double>(pars(87), pars(88));
   Ve(0,1) = std::complex<double>(pars(89), pars(90));
   Ve(0,2) = std::complex<double>(pars(91), pars(92));
   Ve(1,0) = std::complex<double>(pars(93), pars(94));
   Ve(1,1) = std::complex<double>(pars(95), pars(96));
   Ve(1,2) = std::complex<double>(pars(97), pars(98));
   Ve(2,0) = std::complex<double>(pars(99), pars(100));
   Ve(2,1) = std::complex<double>(pars(101), pars(102));
   Ve(2,2) = std::complex<double>(pars(103), pars(104));
   Ue(0,0) = std::complex<double>(pars(105), pars(106));
   Ue(0,1) = std::complex<double>(pars(107), pars(108));
   Ue(0,2) = std::complex<double>(pars(109), pars(110));
   Ue(1,0) = std::complex<double>(pars(111), pars(112));
   Ue(1,1) = std::complex<double>(pars(113), pars(114));
   Ue(1,2) = std::complex<double>(pars(115), pars(116));
   Ue(2,0) = std::complex<double>(pars(117), pars(118));
   Ue(2,1) = std::complex<double>(pars(119), pars(120));
   Ue(2,2) = std::complex<double>(pars(121), pars(122));
   UV(0,0) = std::complex<double>(pars(123), pars(124));
   UV(0,1) = std::complex<double>(pars(125), pars(126));
   UV(0,2) = std::complex<double>(pars(127), pars(128));
   UV(0,3) = std::complex<double>(pars(129), pars(130));
   UV(0,4) = std::complex<double>(pars(131), pars(132));
   UV(0,5) = std::complex<double>(pars(133), pars(134));
   UV(1,0) = std::complex<double>(pars(135), pars(136));
   UV(1,1) = std::complex<double>(pars(137), pars(138));
   UV(1,2) = std::complex<double>(pars(139), pars(140));
   UV(1,3) = std::complex<double>(pars(141), pars(142));
   UV(1,4) = std::complex<double>(pars(143), pars(144));
   UV(1,5) = std::complex<double>(pars(145), pars(146));
   UV(2,0) = std::complex<double>(pars(147), pars(148));
   UV(2,1) = std::complex<double>(pars(149), pars(150));
   UV(2,2) = std::complex<double>(pars(151), pars(152));
   UV(2,3) = std::complex<double>(pars(153), pars(154));
   UV(2,4) = std::complex<double>(pars(155), pars(156));
   UV(2,5) = std::complex<double>(pars(157), pars(158));
   UV(3,0) = std::complex<double>(pars(159), pars(160));
   UV(3,1) = std::complex<double>(pars(161), pars(162));
   UV(3,2) = std::complex<double>(pars(163), pars(164));
   UV(3,3) = std::complex<double>(pars(165), pars(166));
   UV(3,4) = std::complex<double>(pars(167), pars(168));
   UV(3,5) = std::complex<double>(pars(169), pars(170));
   UV(4,0) = std::complex<double>(pars(171), pars(172));
   UV(4,1) = std::complex<double>(pars(173), pars(174));
   UV(4,2) = std::complex<double>(pars(175), pars(176));
   UV(4,3) = std::complex<double>(pars(177), pars(178));
   UV(4,4) = std::complex<double>(pars(179), pars(180));
   UV(4,5) = std::complex<double>(pars(181), pars(182));
   UV(5,0) = std::complex<double>(pars(183), pars(184));
   UV(5,1) = std::complex<double>(pars(185), pars(186));
   UV(5,2) = std::complex<double>(pars(187), pars(188));
   UV(5,3) = std::complex<double>(pars(189), pars(190));
   UV(5,4) = std::complex<double>(pars(191), pars(192));
   UV(5,5) = std::complex<double>(pars(193), pars(194));

}

Eigen::ArrayXd NMFSU5_physical::get_masses() const
{
   Eigen::ArrayXd pars(15);

   pars(0) = MFv(0);
   pars(1) = MFv(1);
   pars(2) = MFv(2);
   pars(3) = MFv(3);
   pars(4) = MFv(4);
   pars(5) = MFv(5);
   pars(6) = MFd(0);
   pars(7) = MFd(1);
   pars(8) = MFd(2);
   pars(9) = MFu(0);
   pars(10) = MFu(1);
   pars(11) = MFu(2);
   pars(12) = MFe(0);
   pars(13) = MFe(1);
   pars(14) = MFe(2);

   return pars;
}

void NMFSU5_physical::set_masses(const Eigen::ArrayXd& pars)
{
   MFv(0) = pars(0);
   MFv(1) = pars(1);
   MFv(2) = pars(2);
   MFv(3) = pars(3);
   MFv(4) = pars(4);
   MFv(5) = pars(5);
   MFd(0) = pars(6);
   MFd(1) = pars(7);
   MFd(2) = pars(8);
   MFu(0)= pars(9);
   MFu(1) = pars(10);
   MFu(2) = pars(11);
   MFe(0) = pars(12);
   MFe(1) = pars(13);
   MFe(2) = pars(14);

}

void NMFSU5_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';

   ostr << "----------------------------------------\n"
           "pole mass mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';
   ostr << "UV = " << UV << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const NMFSU5_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
