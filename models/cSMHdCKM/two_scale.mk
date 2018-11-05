#  ====================================================================
#  This file is part of FlexibleSUSY.
#
#  FlexibleSUSY is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  FlexibleSUSY is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FlexibleSUSY.  If not, see
#  <http://www.gnu.org/licenses/>.
#  ====================================================================

cSMHdCKM_INCLUDE_MK += $(DIR)/two_scale.mk

LIBcSMHdCKM_SRC += \
		$(DIR)/cSMHdCKM_two_scale_convergence_tester.cpp \
		$(DIR)/cSMHdCKM_two_scale_ewsb_solver.cpp \
		$(DIR)/cSMHdCKM_two_scale_high_scale_constraint.cpp \
		$(DIR)/cSMHdCKM_two_scale_initial_guesser.cpp \
		$(DIR)/cSMHdCKM_two_scale_low_scale_constraint.cpp \
		$(DIR)/cSMHdCKM_two_scale_model.cpp \
		$(DIR)/cSMHdCKM_two_scale_spectrum_generator.cpp \
		$(DIR)/cSMHdCKM_two_scale_susy_scale_constraint.cpp
LIBcSMHdCKM_HDR += \
		$(DIR)/cSMHdCKM_two_scale_convergence_tester.hpp \
		$(DIR)/cSMHdCKM_two_scale_ewsb_solver.hpp \
		$(DIR)/cSMHdCKM_two_scale_high_scale_constraint.hpp \
		$(DIR)/cSMHdCKM_two_scale_initial_guesser.hpp \
		$(DIR)/cSMHdCKM_two_scale_low_scale_constraint.hpp \
		$(DIR)/cSMHdCKM_two_scale_model.hpp \
		$(DIR)/cSMHdCKM_two_scale_spectrum_generator.hpp \
		$(DIR)/cSMHdCKM_two_scale_susy_scale_constraint.hpp

LIBcSMHdCKM_SRC += \


LIBcSMHdCKM_HDR += \

