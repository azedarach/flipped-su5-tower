DIR          := models/cSMHdCKMRHN
MODNAME      := cSMHdCKMRHN
SARAH_MODEL  := SMHd
WITH_$(MODNAME) := yes
MODcSMHdCKMRHN_MOD := SM
MODcSMHdCKMRHN_DEP := $(patsubst %,model_specific/%,$(MODcSMHdCKMRHN_MOD))
MODcSMHdCKMRHN_INC := $(patsubst %,-Imodel_specific/%,$(MODcSMHdCKMRHN_MOD))
MODcSMHdCKMRHN_LIB := $(foreach M,$(MODcSMHdCKMRHN_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODcSMHdCKMRHN_SUBMOD  := $(DIR)/cxx_qft
MODcSMHdCKMRHN_SUBMOD_INC := $(patsubst %,-I%,$(MODcSMHdCKMRHN_SUBMOD))

cSMHdCKMRHN_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
cSMHdCKMRHN_INSTALL_CXXQFT_DIR := \
		$(cSMHdCKMRHN_INSTALL_DIR)/cxx_qft

cSMHdCKMRHN_MK     := \
		$(DIR)/module.mk

cSMHdCKMRHN_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

cSMHdCKMRHN_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

cSMHdCKMRHN_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

cSMHdCKMRHN_INCLUDE_MK := \
		$(cSMHdCKMRHN_SUSY_BETAS_MK) \
		$(cSMHdCKMRHN_SOFT_BETAS_MK)

cSMHdCKMRHN_SLHA_INPUT := \
		$(DIR)/LesHouches.in.cSMHdCKMRHN_generated \


cSMHdCKMRHN_REFERENCES := \
		$(DIR)/cSMHdCKMRHN_references.tex

cSMHdCKMRHN_GNUPLOT := \
		$(DIR)/cSMHdCKMRHN_plot_rgflow.gnuplot \
		$(DIR)/cSMHdCKMRHN_plot_spectrum.gnuplot

cSMHdCKMRHN_TARBALL := \
		$(MODNAME).tar.gz

LIBcSMHdCKMRHN_SRC := \
		$(DIR)/cSMHdCKMRHN_a_muon.cpp \
		$(DIR)/cSMHdCKMRHN_edm.cpp \
		$(DIR)/cSMHdCKMRHN_effective_couplings.cpp \
		$(DIR)/cSMHdCKMRHN_info.cpp \
		$(DIR)/cSMHdCKMRHN_input_parameters.cpp \
		$(DIR)/cSMHdCKMRHN_mass_eigenstates.cpp \
		$(DIR)/cSMHdCKMRHN_observables.cpp \
		$(DIR)/cSMHdCKMRHN_physical.cpp \
		$(DIR)/cSMHdCKMRHN_slha_io.cpp \
		$(DIR)/cSMHdCKMRHN_soft_parameters.cpp \
		$(DIR)/cSMHdCKMRHN_susy_parameters.cpp \
		$(DIR)/cSMHdCKMRHN_utilities.cpp \
		$(DIR)/cSMHdCKMRHN_weinberg_angle.cpp

EXEcSMHdCKMRHN_SRC := \
		$(DIR)/run_cSMHdCKMRHN.cpp \
		$(DIR)/run_cmd_line_cSMHdCKMRHN.cpp \
		$(DIR)/scan_cSMHdCKMRHN.cpp
LLcSMHdCKMRHN_LIB  :=
LLcSMHdCKMRHN_OBJ  :=
LLcSMHdCKMRHN_SRC  := \
		$(DIR)/cSMHdCKMRHN_librarylink.cpp

LLcSMHdCKMRHN_MMA  := \
		$(DIR)/cSMHdCKMRHN_librarylink.m \
		$(DIR)/run_cSMHdCKMRHN.m

LIBcSMHdCKMRHN_HDR := \
		$(DIR)/cSMHdCKMRHN_a_muon.hpp \
		$(DIR)/cSMHdCKMRHN_convergence_tester.hpp \
		$(DIR)/cSMHdCKMRHN_edm.hpp \
		$(DIR)/cSMHdCKMRHN_effective_couplings.hpp \
		$(DIR)/cSMHdCKMRHN_ewsb_solver.hpp \
		$(DIR)/cSMHdCKMRHN_ewsb_solver_interface.hpp \
		$(DIR)/cSMHdCKMRHN_high_scale_constraint.hpp \
		$(DIR)/cSMHdCKMRHN_info.hpp \
		$(DIR)/cSMHdCKMRHN_initial_guesser.hpp \
		$(DIR)/cSMHdCKMRHN_input_parameters.hpp \
		$(DIR)/cSMHdCKMRHN_low_scale_constraint.hpp \
		$(DIR)/cSMHdCKMRHN_mass_eigenstates.hpp \
		$(DIR)/cSMHdCKMRHN_model.hpp \
		$(DIR)/cSMHdCKMRHN_model_slha.hpp \
		$(DIR)/cSMHdCKMRHN_observables.hpp \
		$(DIR)/cSMHdCKMRHN_physical.hpp \
		$(DIR)/cSMHdCKMRHN_slha_io.hpp \
		$(DIR)/cSMHdCKMRHN_spectrum_generator.hpp \
		$(DIR)/cSMHdCKMRHN_spectrum_generator_interface.hpp \
		$(DIR)/cSMHdCKMRHN_soft_parameters.hpp \
		$(DIR)/cSMHdCKMRHN_susy_parameters.hpp \
		$(DIR)/cSMHdCKMRHN_susy_scale_constraint.hpp \
		$(DIR)/cSMHdCKMRHN_utilities.hpp \
		$(DIR)/cSMHdCKMRHN_weinberg_angle.hpp

LIBcSMHdCKMRHN_CXXQFT_HDR := \
		$(DIR)/cxx_qft/cSMHdCKMRHN_qft.hpp \
		$(DIR)/cxx_qft/cSMHdCKMRHN_fields.hpp \
		$(DIR)/cxx_qft/cSMHdCKMRHN_vertices.hpp \
		$(DIR)/cxx_qft/cSMHdCKMRHN_context_base.hpp

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
-include $(cSMHdCKMRHN_SUSY_BETAS_MK)
-include $(cSMHdCKMRHN_SOFT_BETAS_MK)
-include $(cSMHdCKMRHN_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(cSMHdCKMRHN_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(cSMHdCKMRHN_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(cSMHdCKMRHN_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBcSMHdCKMRHN_SRC := $(sort $(LIBcSMHdCKMRHN_SRC))
EXEcSMHdCKMRHN_SRC := $(sort $(EXEcSMHdCKMRHN_SRC))

LIBcSMHdCKMRHN_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBcSMHdCKMRHN_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBcSMHdCKMRHN_SRC)))

EXEcSMHdCKMRHN_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEcSMHdCKMRHN_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEcSMHdCKMRHN_SRC)))

EXEcSMHdCKMRHN_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEcSMHdCKMRHN_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEcSMHdCKMRHN_SRC)))

LIBcSMHdCKMRHN_DEP := \
		$(LIBcSMHdCKMRHN_OBJ:.o=.d)

EXEcSMHdCKMRHN_DEP := \
		$(EXEcSMHdCKMRHN_OBJ:.o=.d)

LLcSMHdCKMRHN_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLcSMHdCKMRHN_SRC)))

LLcSMHdCKMRHN_OBJ  := $(LLcSMHdCKMRHN_SRC:.cpp=.o)
LLcSMHdCKMRHN_LIB  := $(LLcSMHdCKMRHN_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBcSMHdCKMRHN     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_cSMHdCKMRHN := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_cSMHdCKMRHN := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBcSMHdCKMRHN) $(EXEcSMHdCKMRHN_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(cSMHdCKMRHN_INSTALL_DIR)
		install -d $(cSMHdCKMRHN_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(LIBcSMHdCKMRHN_SRC) $(cSMHdCKMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBcSMHdCKMRHN_HDR) $(cSMHdCKMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBcSMHdCKMRHN_CXXQFT_HDR) $(cSMHdCKMRHN_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(EXEcSMHdCKMRHN_SRC) $(cSMHdCKMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLcSMHdCKMRHN_SRC) $(cSMHdCKMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLcSMHdCKMRHN_MMA) $(cSMHdCKMRHN_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(cSMHdCKMRHN_MK) $(cSMHdCKMRHN_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(cSMHdCKMRHN_INCLUDE_MK) $(cSMHdCKMRHN_INSTALL_DIR)
ifneq ($(cSMHdCKMRHN_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(cSMHdCKMRHN_SLHA_INPUT) $(cSMHdCKMRHN_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(cSMHdCKMRHN_REFERENCES) $(cSMHdCKMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(cSMHdCKMRHN_GNUPLOT) $(cSMHdCKMRHN_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBcSMHdCKMRHN_DEP)
		-rm -f $(EXEcSMHdCKMRHN_DEP)
		-rm -f $(LLcSMHdCKMRHN_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBcSMHdCKMRHN)
		-rm -f $(LLcSMHdCKMRHN_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBcSMHdCKMRHN_OBJ)
		-rm -f $(EXEcSMHdCKMRHN_OBJ)
		-rm -f $(LLcSMHdCKMRHN_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEcSMHdCKMRHN_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(cSMHdCKMRHN_TARBALL) \
		$(LIBcSMHdCKMRHN_SRC) $(LIBcSMHdCKMRHN_HDR) $(LIBcSMHdCKMRHN_CXXQFT_HDR) \
		$(EXEcSMHdCKMRHN_SRC) \
		$(LLcSMHdCKMRHN_SRC) $(LLcSMHdCKMRHN_MMA) \
		$(cSMHdCKMRHN_MK) $(cSMHdCKMRHN_INCLUDE_MK) \
		$(cSMHdCKMRHN_SLHA_INPUT) $(cSMHdCKMRHN_REFERENCES) \
		$(cSMHdCKMRHN_GNUPLOT)

$(LIBcSMHdCKMRHN_SRC) $(LIBcSMHdCKMRHN_HDR) $(LIBcSMHdCKMRHN_CXXQFT_HDR) $(EXEcSMHdCKMRHN_SRC) $(LLcSMHdCKMRHN_SRC) $(LLcSMHdCKMRHN_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_cSMHdCKMRHN)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_cSMHdCKMRHN): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_cSMHdCKMRHN)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_cSMHdCKMRHN)"
		@echo "Note: to regenerate cSMHdCKMRHN source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_cSMHdCKMRHN)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_cSMHdCKMRHN):
		@true
endif

$(LIBcSMHdCKMRHN_DEP) $(EXEcSMHdCKMRHN_DEP) $(LLcSMHdCKMRHN_DEP) $(LIBcSMHdCKMRHN_OBJ) $(EXEcSMHdCKMRHN_OBJ) $(LLcSMHdCKMRHN_OBJ) $(LLcSMHdCKMRHN_LIB): \
	CPPFLAGS += $(MODcSMHdCKMRHN_SUBMOD_INC) $(MODcSMHdCKMRHN_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBcSMHdCKMRHN_DEP) $(EXEcSMHdCKMRHN_DEP) $(LLcSMHdCKMRHN_DEP) $(LIBcSMHdCKMRHN_OBJ) $(EXEcSMHdCKMRHN_OBJ) $(LLcSMHdCKMRHN_OBJ) $(LLcSMHdCKMRHN_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLcSMHdCKMRHN_OBJ) $(LLcSMHdCKMRHN_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBcSMHdCKMRHN): $(LIBcSMHdCKMRHN_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBcSMHdCKMRHN) $(MODcSMHdCKMRHN_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLcSMHdCKMRHN_LIB): $(LLcSMHdCKMRHN_OBJ) $(LIBcSMHdCKMRHN) $(MODcSMHdCKMRHN_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBcSMHdCKMRHN_DEP) $(EXEcSMHdCKMRHN_DEP)
ALLSRC += $(LIBcSMHdCKMRHN_SRC) $(EXEcSMHdCKMRHN_SRC)
ALLLIB += $(LIBcSMHdCKMRHN)
ALLEXE += $(EXEcSMHdCKMRHN_EXE)
ALLMODDEP += $(MODcSMHdCKMRHN_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLcSMHdCKMRHN_DEP)
ALLSRC += $(LLcSMHdCKMRHN_SRC)
ALLLL  += $(LLcSMHdCKMRHN_LIB)
endif
