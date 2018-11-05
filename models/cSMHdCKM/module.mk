DIR          := models/cSMHdCKM
MODNAME      := cSMHdCKM
SARAH_MODEL  := SMHd
WITH_$(MODNAME) := yes
MODcSMHdCKM_MOD := SM
MODcSMHdCKM_DEP := $(patsubst %,model_specific/%,$(MODcSMHdCKM_MOD))
MODcSMHdCKM_INC := $(patsubst %,-Imodel_specific/%,$(MODcSMHdCKM_MOD))
MODcSMHdCKM_LIB := $(foreach M,$(MODcSMHdCKM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODcSMHdCKM_SUBMOD  := $(DIR)/cxx_qft
MODcSMHdCKM_SUBMOD_INC := $(patsubst %,-I%,$(MODcSMHdCKM_SUBMOD))

cSMHdCKM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
cSMHdCKM_INSTALL_CXXQFT_DIR := \
		$(cSMHdCKM_INSTALL_DIR)/cxx_qft

cSMHdCKM_MK     := \
		$(DIR)/module.mk

cSMHdCKM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

cSMHdCKM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

cSMHdCKM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

cSMHdCKM_INCLUDE_MK := \
		$(cSMHdCKM_SUSY_BETAS_MK) \
		$(cSMHdCKM_SOFT_BETAS_MK)

cSMHdCKM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.cSMHdCKM_generated \
		$(DIR)/LesHouches.in.cSMHdCKM

cSMHdCKM_REFERENCES := \
		$(DIR)/cSMHdCKM_references.tex

cSMHdCKM_GNUPLOT := \
		$(DIR)/cSMHdCKM_plot_rgflow.gnuplot \
		$(DIR)/cSMHdCKM_plot_spectrum.gnuplot

cSMHdCKM_TARBALL := \
		$(MODNAME).tar.gz

LIBcSMHdCKM_SRC := \
		$(DIR)/cSMHdCKM_a_muon.cpp \
		$(DIR)/cSMHdCKM_edm.cpp \
		$(DIR)/cSMHdCKM_effective_couplings.cpp \
		$(DIR)/cSMHdCKM_info.cpp \
		$(DIR)/cSMHdCKM_input_parameters.cpp \
		$(DIR)/cSMHdCKM_mass_eigenstates.cpp \
		$(DIR)/cSMHdCKM_observables.cpp \
		$(DIR)/cSMHdCKM_physical.cpp \
		$(DIR)/cSMHdCKM_slha_io.cpp \
		$(DIR)/cSMHdCKM_soft_parameters.cpp \
		$(DIR)/cSMHdCKM_susy_parameters.cpp \
		$(DIR)/cSMHdCKM_utilities.cpp \
		$(DIR)/cSMHdCKM_weinberg_angle.cpp

EXEcSMHdCKM_SRC := \
		$(DIR)/run_cSMHdCKM.cpp \
		$(DIR)/run_cmd_line_cSMHdCKM.cpp \
		$(DIR)/scan_cSMHdCKM.cpp
LLcSMHdCKM_LIB  :=
LLcSMHdCKM_OBJ  :=
LLcSMHdCKM_SRC  := \
		$(DIR)/cSMHdCKM_librarylink.cpp

LLcSMHdCKM_MMA  := \
		$(DIR)/cSMHdCKM_librarylink.m \
		$(DIR)/run_cSMHdCKM.m

LIBcSMHdCKM_HDR := \
		$(DIR)/cSMHdCKM_a_muon.hpp \
		$(DIR)/cSMHdCKM_convergence_tester.hpp \
		$(DIR)/cSMHdCKM_edm.hpp \
		$(DIR)/cSMHdCKM_effective_couplings.hpp \
		$(DIR)/cSMHdCKM_ewsb_solver.hpp \
		$(DIR)/cSMHdCKM_ewsb_solver_interface.hpp \
		$(DIR)/cSMHdCKM_high_scale_constraint.hpp \
		$(DIR)/cSMHdCKM_info.hpp \
		$(DIR)/cSMHdCKM_initial_guesser.hpp \
		$(DIR)/cSMHdCKM_input_parameters.hpp \
		$(DIR)/cSMHdCKM_low_scale_constraint.hpp \
		$(DIR)/cSMHdCKM_mass_eigenstates.hpp \
		$(DIR)/cSMHdCKM_model.hpp \
		$(DIR)/cSMHdCKM_model_slha.hpp \
		$(DIR)/cSMHdCKM_observables.hpp \
		$(DIR)/cSMHdCKM_physical.hpp \
		$(DIR)/cSMHdCKM_slha_io.hpp \
		$(DIR)/cSMHdCKM_spectrum_generator.hpp \
		$(DIR)/cSMHdCKM_spectrum_generator_interface.hpp \
		$(DIR)/cSMHdCKM_soft_parameters.hpp \
		$(DIR)/cSMHdCKM_susy_parameters.hpp \
		$(DIR)/cSMHdCKM_susy_scale_constraint.hpp \
		$(DIR)/cSMHdCKM_utilities.hpp \
		$(DIR)/cSMHdCKM_weinberg_angle.hpp

LIBcSMHdCKM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/cSMHdCKM_qft.hpp \
		$(DIR)/cxx_qft/cSMHdCKM_fields.hpp \
		$(DIR)/cxx_qft/cSMHdCKM_vertices.hpp \
		$(DIR)/cxx_qft/cSMHdCKM_context_base.hpp

ifneq ($(findstring two_scale,$(SOLVERS)),)
-include $(DIR)/two_scale.mk
endif
ifneq ($(findstring lattice,$(SOLVERS)),)
-include $(DIR)/lattice.mk
endif
ifneq ($(findstring semi_analytic,$(SOLVERS)),)
-include $(DIR)/semi_analytic.mk
endif

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(cSMHdCKM_SUSY_BETAS_MK)
-include $(cSMHdCKM_SOFT_BETAS_MK)
-include $(cSMHdCKM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(cSMHdCKM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(cSMHdCKM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(cSMHdCKM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

# remove duplicates in case all solvers are used
LIBcSMHdCKM_SRC := $(sort $(LIBcSMHdCKM_SRC))
EXEcSMHdCKM_SRC := $(sort $(EXEcSMHdCKM_SRC))

LIBcSMHdCKM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBcSMHdCKM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBcSMHdCKM_SRC)))

EXEcSMHdCKM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEcSMHdCKM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEcSMHdCKM_SRC)))

EXEcSMHdCKM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEcSMHdCKM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEcSMHdCKM_SRC)))

LIBcSMHdCKM_DEP := \
		$(LIBcSMHdCKM_OBJ:.o=.d)

EXEcSMHdCKM_DEP := \
		$(EXEcSMHdCKM_OBJ:.o=.d)

LLcSMHdCKM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLcSMHdCKM_SRC)))

LLcSMHdCKM_OBJ  := $(LLcSMHdCKM_SRC:.cpp=.o)
LLcSMHdCKM_LIB  := $(LLcSMHdCKM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBcSMHdCKM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_cSMHdCKM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_cSMHdCKM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBcSMHdCKM) $(EXEcSMHdCKM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(cSMHdCKM_INSTALL_DIR)
		install -d $(cSMHdCKM_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(LIBcSMHdCKM_SRC) $(cSMHdCKM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBcSMHdCKM_HDR) $(cSMHdCKM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBcSMHdCKM_CXXQFT_HDR) $(cSMHdCKM_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(EXEcSMHdCKM_SRC) $(cSMHdCKM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLcSMHdCKM_SRC) $(cSMHdCKM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLcSMHdCKM_MMA) $(cSMHdCKM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(cSMHdCKM_MK) $(cSMHdCKM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(cSMHdCKM_INCLUDE_MK) $(cSMHdCKM_INSTALL_DIR)
ifneq ($(cSMHdCKM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(cSMHdCKM_SLHA_INPUT) $(cSMHdCKM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(cSMHdCKM_REFERENCES) $(cSMHdCKM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(cSMHdCKM_GNUPLOT) $(cSMHdCKM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBcSMHdCKM_DEP)
		-rm -f $(EXEcSMHdCKM_DEP)
		-rm -f $(LLcSMHdCKM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBcSMHdCKM)
		-rm -f $(LLcSMHdCKM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBcSMHdCKM_OBJ)
		-rm -f $(EXEcSMHdCKM_OBJ)
		-rm -f $(LLcSMHdCKM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEcSMHdCKM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(cSMHdCKM_TARBALL) \
		$(LIBcSMHdCKM_SRC) $(LIBcSMHdCKM_HDR) $(LIBcSMHdCKM_CXXQFT_HDR) \
		$(EXEcSMHdCKM_SRC) \
		$(LLcSMHdCKM_SRC) $(LLcSMHdCKM_MMA) \
		$(cSMHdCKM_MK) $(cSMHdCKM_INCLUDE_MK) \
		$(cSMHdCKM_SLHA_INPUT) $(cSMHdCKM_REFERENCES) \
		$(cSMHdCKM_GNUPLOT)

$(LIBcSMHdCKM_SRC) $(LIBcSMHdCKM_HDR) $(LIBcSMHdCKM_CXXQFT_HDR) $(EXEcSMHdCKM_SRC) $(LLcSMHdCKM_SRC) $(LLcSMHdCKM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_cSMHdCKM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_cSMHdCKM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_cSMHdCKM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_cSMHdCKM)"
		@echo "Note: to regenerate cSMHdCKM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_cSMHdCKM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_cSMHdCKM):
		@true
endif

$(LIBcSMHdCKM_DEP) $(EXEcSMHdCKM_DEP) $(LLcSMHdCKM_DEP) $(LIBcSMHdCKM_OBJ) $(EXEcSMHdCKM_OBJ) $(LLcSMHdCKM_OBJ) $(LLcSMHdCKM_LIB): \
	CPPFLAGS += $(MODcSMHdCKM_SUBMOD_INC) $(MODcSMHdCKM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBcSMHdCKM_DEP) $(EXEcSMHdCKM_DEP) $(LLcSMHdCKM_DEP) $(LIBcSMHdCKM_OBJ) $(EXEcSMHdCKM_OBJ) $(LLcSMHdCKM_OBJ) $(LLcSMHdCKM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLcSMHdCKM_OBJ) $(LLcSMHdCKM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBcSMHdCKM): $(LIBcSMHdCKM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBcSMHdCKM) $(MODcSMHdCKM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLcSMHdCKM_LIB): $(LLcSMHdCKM_OBJ) $(LIBcSMHdCKM) $(MODcSMHdCKM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBcSMHdCKM_DEP) $(EXEcSMHdCKM_DEP)
ALLSRC += $(LIBcSMHdCKM_SRC) $(EXEcSMHdCKM_SRC)
ALLLIB += $(LIBcSMHdCKM)
ALLEXE += $(EXEcSMHdCKM_EXE)
ALLMODDEP += $(MODcSMHdCKM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLcSMHdCKM_DEP)
ALLSRC += $(LLcSMHdCKM_SRC)
ALLLL  += $(LLcSMHdCKM_LIB)
endif
