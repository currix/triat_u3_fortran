######################################################################
#	Makefile to compile, print or check the U(3) programs that
#	compute bending vibrational spectra for linear/bent
#	triatomic molecules using the two-dimensional limit of the
#	Vibron Model.
#
#       by Currix TM.
#######################################################################
#
# 
BINPATH = ../bin
#
VERSION = "2.2.3"
#
# General 2-3-4 body Hamiltonian
#
#
# PROGRAM THAT CARRY OUT FITTING TO EXPERIMENTAL ENERGY SPECTRA USING MINUIT
TRIMINSRCT = chisqre_u3_min_gen.f readenerg_gfort.f \
          hbldu3_gen.f assign_gen.f scaleh_gen.f \
          maxc.f maxc_exp.f cpexp.f disdata.f \
          selectvl.f minuit-cern.f intrac.f \
          so3casimirW2.f  so3casimirW4.f so3casimirW2WB2.f \
          fcn_sub_gen.f minuit_u3.f
#
TRIMINOBJT = $(TRIMINSRCT:.f=.o) 
###############################################################
# ALGEBRAIC MODEL ENERGIES (GENERAL HAMILTONIAN)
ENRGENSRCT = enrgy_gen_U3.f hbldu3_gen.f  \
          so3casimirW2.f  so3casimirW4.f so3casimirW2WB2.f
#
ENRGENOBJT = $(ENRGENSRCT:.f=.o)
###############################################################
# ALGEBRAIC MODEL ENERGIES AND EIGENVALUES (GENERAL HAMILTONIAN)
AVECGENSRCT = avalavec_gen_U3.f hbldu3_gen.f  \
          so3casimirW2.f so3casimirW4.f so3casimirW2WB2.f change_basis.f check_phase.f
#
AVECGENOBJT = $(AVECGENSRCT:.f=.o)
#########################################################################
# ALGEBRAIC MODEL IPR, ENERGIES, AND EIGENVALUES (GENERAL HAMILTONIAN)
AVECIPRGENSRCT = avalavec_gen_U3_ipr.f hbldu3_gen.f  \
          so3casimirW2.f so3casimirW4.f so3casimirW2WB2.f change_basis.f check_phase.f
#
AVECIPRGENOBJT = $(AVECIPRGENSRCT:.f=.o)
###############################################################
# ALGEBRAIC MODEL ENERGIES (MODEL HAMILTONIAN)
#
ENRGSRCTRI = enrgy_modelH_U3_tri.f hbldu3_modelH_tri.f
#
ENRGOBJTRI = $(ENRGSRCTRI:.f=.o)
##
ENRGSRC = enrgy_modelH_U3.f hbldu3_modelH.f so3casimirW2.f 
#
ENRGOBJ = $(ENRGSRC:.f=.o)
####################################################################
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (MODEL HAMILTONIAN)
#
AVECSRC = avalavec_modelH_U3.f hbldu3_modelH.f  \
          so3casimirW2.f
AVECOBJ = $(AVECSRC:.f=.o)
#
####################################################################
# ALGEBRAIC MODEL INTENSITIES (INFRARED CASE)
#
IRINOBJ = ir_intensity.o 
IRINSRC = ir_intensity.f90
#
####################################################################
#
FC = gfortran
#
FOPT	= -c 
#FOPT   += -Wall
#
LOPT     =  -L /usr/lib/atlas -lblas -llapack
LSTATMIN =  -L /usr/lib/atlas -L /usr/lib/atlas-base/atlas -llapack -lblas -lpthread
LSTAT   =   -L /usr/lib/atlas -L /usr/lib/atlas-base/atlas -llapack -lblas -lpthread
all: triat_min_T energenh_T avecgenh_T energmodh energmodhtrid avecmodh IR_inten avecgenh_T_ipr
#
stats : triat_min_Tstat energenh_Tstat avecgenh_Tstat energmodh_stat energmodhtrid_stat avecmodh_stat IR_inten_stat
#
.PHONY : all
# Standard compilation rules
#
%.o: %.f90
	$(FC) $(FOPT) $<
#
%.o: %.F90
	$(FC) $(FOPT) $<
#
%.o: %.f
	$(FC) $(FOPT) $<
#########################################################
# U3 general hamiltonian program with minuit minimization
triat_min_T: $(TRIMINOBJT) 
	$(FC) -o $(BINPATH)/tri_min_T $(TRIMINOBJT) $(LOPT) 
triat_min_Tstat: $(TRIMINOBJT) 
	$(FC) -static -o $(BINPATH)/tri_min_Tstat $(TRIMINOBJT) $(LSTATMIN)
# HAMILTONIAN EIGENVALUES (GENERAL HAMILTONIAN)
energenh_T: $(ENRGENOBJT)
	$(FC) -o $(BINPATH)/en_U3_gen_T $(ENRGENOBJT) $(LOPT)
energenh_Tstat: $(ENRGENOBJT)
	$(FC) -static -o $(BINPATH)/en_U3_gen_Tstat $(ENRGENOBJT) $(LSTAT)
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (GENERAL HAMILTONIAN)
avecgenh_T: $(AVECGENOBJT)
	$(FC) -o $(BINPATH)/avec_U3_gen_T $(AVECGENOBJT) $(LOPT)
avecgenh_T_ipr: $(AVECIPRGENOBJT)
	$(FC) -o $(BINPATH)/avec_U3_gen_T_ipr $(AVECIPRGENOBJT) $(LOPT)
avecgenh_Tstat: $(AVECGENOBJT)
	$(FC) -static -o $(BINPATH)/avec_U3_gen_Tstat $(AVECGENOBJT) $(LSTAT)
# HAMILTONIAN EIGENVALUES (MODEL HAMILTONIAN)
energmodh: $(ENRGOBJ)
	$(FC) -o $(BINPATH)/en_U3_mh $(ENRGOBJ) $(LOPT)
energmodhtrid: $(ENRGOBJTRI)
	$(FC) -o $(BINPATH)/en_U3_mh_trid $(ENRGOBJTRI) $(LOPT)
energmodh_stat: $(ENRGOBJ)
	$(FC) -static -o $(BINPATH)/en_U3_mh_stat $(ENRGOBJ) $(LSTAT)
energmodhtrid_stat: $(ENRGOBJTRI)
	$(FC) -static -o $(BINPATH)/en_U3_mh_trid_stat $(ENRGOBJTRI) $(LSTAT)
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN)
avecmodh: $(AVECOBJ)
	$(FC) -o $(BINPATH)/avec_U3_mh $(AVECOBJ) $(LOPT)
avecmodh_stat: $(AVECOBJ)
	$(FC) -static -o $(BINPATH)/avec_U3_mh_stat $(AVECOBJ) $(LSTAT)
# INFRARED INTENSITY (F90 CODE)
IR_inten: $(IRINOBJ)
	$(FC) -o $(BINPATH)/ir_intensity $(IRINOBJ)
IR_inten_stat: $(IRINOBJ)
	$(FC) -static -o $(BINPATH)/ir_intensity_stat $(IRINOBJ)
#################################################################
install:
	@cp $(BINPATH)/alg_U3_energy  /usr/local/bin/alg_U3_energy_$(VERSION)
	@ln -fs /usr/local/bin/alg_U3_energy_$(VERSION) /usr/local/bin/alg_U3_energy
#
	@cp $(BINPATH)/avec_U3_gen_T  /usr/local/bin/avec_U3_gen_T_$(VERSION)
	@ln -fs /usr/local/bin/avec_U3_gen_T_$(VERSION) /usr/local/bin/avec_U3_gen_T
	@cp $(BINPATH)/avec_U3_gen_Tstat /usr/local/bin/avec_U3_gen_Tstat_$(VERSION)
	@ln -fs /usr/local/bin/avec_U3_gen_Tstat_$(VERSION) /usr/local/bin/avec_U3_gen_Tstat
#
	@cp $(BINPATH)/avec_U3_mh /usr/local/bin/avec_U3_mh_$(VERSION)
	@ln -fs /usr/local/bin/avec_U3_mh_$(VERSION) /usr/local/bin/avec_U3_mh
	@cp $(BINPATH)/avec_U3_mh_stat /usr/local/bin/avec_U3_mh_stat_$(VERSION)
	@ln -fs /usr/local/bin/avec_U3_mh_stat_$(VERSION) /usr/local/bin/avec_U3_mh_stat
#
	@cp $(BINPATH)/en_U3_gen_T  /usr/local/bin/en_U3_gen_T_$(VERSION)
	@ln -fs /usr/local/bin/en_U3_gen_T_$(VERSION) /usr/local/bin/en_U3_gen_T
	@cp $(BINPATH)/en_U3_gen_Tstat /usr/local/bin/en_U3_gen_Tstat_$(VERSION)
	@ln -fs /usr/local/bin/en_U3_gen_Tstat_$(VERSION) /usr/local/bin/en_U3_gen_Tstat
#
	@cp $(BINPATH)/en_U3_mh /usr/local/bin/en_U3_mh_$(VERSION)
	@ln -fs /usr/local/bin/en_U3_mh_$(VERSION) /usr/local/bin/en_U3_mh
	@cp $(BINPATH)/en_U3_mh_stat /usr/local/bin/en_U3_mh_stat_$(VERSION)
	@ln -fs /usr/local/bin/en_U3_mh_stat_$(VERSION) /usr/local/bin/en_U3_mh_stat
#
	@cp $(BINPATH)/en_U3_mh_trid /usr/local/bin/en_U3_mh_trid_$(VERSION)
	@ln -fs /usr/local/bin/en_U3_mh_trid_$(VERSION) /usr/local/bin/en_U3_mh_trid
	@cp $(BINPATH)/en_U3_mh_trid_stat /usr/local/bin/en_U3_mh_trid_stat_$(VERSION)
	@ln -fs /usr/local/bin/en_U3_mh_trid_stat_$(VERSION) /usr/local/bin/en_U3_mh_trid_stat
#
	@cp $(BINPATH)/exc_energ_diag_U3_mh  /usr/local/bin/exc_energ_diag_U3_mh_$(VERSION)
	@ln -fs /usr/local/bin/exc_energ_diag_U3_mh_$(VERSION) /usr/local/bin/exc_energ_diag_U3_mh
#
	@cp $(BINPATH)/expected_val_n /usr/local/bin/expected_val_n_$(VERSION) 
	@ln -fs /usr/local/bin/expected_val_n_$(VERSION) /usr/local/bin/expected_val_n
#
	@cp $(BINPATH)/intensity_IR  /usr/local/bin/intensity_IR_$(VERSION)
	@ln -fs /usr/local/bin/intensity_IR_$(VERSION) /usr/local/bin/intensity_IR
#
	@cp $(BINPATH)/ir_intensity /usr/local/bin/ir_intensity_$(VERSION)
	@ln -fs /usr/local/bin/ir_intensity_$(VERSION) /usr/local/bin/ir_intensity
	@cp $(BINPATH)/ir_intensity_stat /usr/local/bin/ir_intensity_stat_$(VERSION)
	@ln -fs /usr/local/bin/ir_intensity_stat_$(VERSION) /usr/local/bin/ir_intensity_stat
#
	@cp $(BINPATH)/minuit_U3_genH  /usr/local/bin/minuit_U3_genH_$(VERSION)
	@ln -fs /usr/local/bin/minuit_U3_genH_$(VERSION) /usr/local/bin/minuit_U3_genH
#
	@cp $(BINPATH)/tri_min_T  /usr/local/bin/tri_min_T_$(VERSION)
	@ln -fs /usr/local/bin/tri_min_T_$(VERSION) /usr/local/bin/tri_min_T
	@cp $(BINPATH)/tri_min_Tstat /usr/local/bin/tri_min_Tstat_$(VERSION)
	@ln -fs /usr/local/bin/tri_min_Tstat_$(VERSION) /usr/local/bin/tri_min_Tstat
#######################################################
uninstall:
	@rm -f  /usr/local/bin/alg_U3_energy /usr/local/bin/avec_U3_gen_T /usr/local/bin/avec_U3_gen_Tstat /usr/local/bin/avec_U3_mh /usr/local/bin/avec_U3_mh_stat /usr/local/bin/en_U3_gen_T /usr/local/bin/en_U3_gen_Tstat /usr/local/bin/en_U3_mh /usr/local/bin/en_U3_mh_stat /usr/local/bin/en_U3_mh_trid /usr/local/bin/en_U3_mh_trid_stat /usr/local/bin/exc_energ_diag_U3_mh /usr/local/bin/expected_val_n /usr/local/bin/intensity_IR /usr/local/bin/ir_intensity /usr/local/bin/ir_intensity_stat /usr/local/bin/minuit_U3_genH  /usr/local/bin/tri_min_T  /usr/local/bin/tri_min_Tstat  
#
	@rm -f  /usr/local/bin/alg_U3_energy_$(VERSION) /usr/local/bin/avec_U3_gen_T_$(VERSION) /usr/local/bin/avec_U3_gen_Tstat_$(VERSION) /usr/local/bin/avec_U3_mh_$(VERSION) /usr/local/bin/avec_U3_mh_stat_$(VERSION) /usr/local/bin/en_U3_gen_T_$(VERSION) /usr/local/bin/en_U3_gen_Tstat_$(VERSION) /usr/local/bin/en_U3_mh_$(VERSION) /usr/local/bin/en_U3_mh_stat_$(VERSION) /usr/local/bin/en_U3_mh_trid_$(VERSION) /usr/local/bin/en_U3_mh_trid_stat_$(VERSION) /usr/local/bin/exc_energ_diag_U3_mh_$(VERSION) /usr/local/bin/expected_val_n_$(VERSION) /usr/local/bin/intensity_IR_$(VERSION) /usr/local/bin/ir_intensity_$(VERSION) /usr/local/bin/ir_intensity_stat_$(VERSION) /usr/local/bin/minuit_U3_genH_$(VERSION) /usr/local/bin/tri_min_T_$(VERSION) /usr/local/bin/tri_min_Tstat_$(VERSION)  
###############
clean:
	@echo "Cleaning object files"
	@rm -f $(TRIMINOBJT) $(ENRGENOBJT) $(AVECGENOBJT) $(ENRGOBJ) $(ENRGOBJTRI) $(AVECOBJ) $(IRINOBJ) $(AVECIPRGENOBJT)





