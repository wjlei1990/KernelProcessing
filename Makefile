SRCDIR=src
OBJDIR=obj
BINDIR=bin

FC=gfortran
MPIFC=mpif90
#MPIFC=ftn

FCFLAGS=-O3 -Wall -J $(OBJDIR) -I $(OBJDIR)

adios_link=$(shell adios_config -lf)
adios_inc=$(shell adios_config -cf)

#adios_link=$(shell /ccs/home/ccui/adios-gcc/adios_config -lf)
#adios_inc=$(shell /ccs/home/ccui/adios-gcc/adios_config -cf)

objects= $(OBJDIR)/adios_helpers_definitions.o $(OBJDIR)/adios_helpers_writers.o $(OBJDIR)/adios_helpers.o $(OBJDIR)/gll_library.o $(OBJDIR)/global_var.o $(OBJDIR)/AdiosIO.o

all: $(BINDIR)/xsteepDescent \
	$(BINDIR)/xcg_direction \
	$(BINDIR)/xlbfgs \
	$(BINDIR)/xsum_kernels \
	$(BINDIR)/xsum_azi_kernels \
	$(BINDIR)/xprecond_kernels \
	$(BINDIR)/xprecond_vp_vs_kernels \
	$(BINDIR)/xprecond_azi_kernels \
	$(BINDIR)/xmerge_kernels \
	$(BINDIR)/xmerge_azi_kernels \
	$(BINDIR)/xupdate_model \
	$(BINDIR)/xmodel_perturb_ref \
	$(BINDIR)/xblend_model \
	$(BINDIR)/xgauss_single \
	$(BINDIR)/xgauss_multiple \
	$(BINDIR)/xbp2binary \
	$(BINDIR)/xascii2bp \
	$(BINDIR)/xabs_kernel \
	$(BINDIR)/xcompute_azi_params \
	$(BINDIR)/xsrc_mask \
	$(BINDIR)/xcompute_vp_vs_hess \
	$(BINDIR)/xrandom_probe_lbfgs


# ###################
# Compile
# ###################
$(OBJDIR)/random.o: $(SRCDIR)/random.f90
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJDIR)/global_var.o: $(SRCDIR)/global_var.f90 $(OBJDIR)/gll_library.o
	$(MPIFC) $(FCFLAGS) -c $< -o $@

$(OBJDIR)/AdiosIO.o: $(SRCDIR)/AdiosIO.f90 $(OBJDIR)/adios_helpers.o
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/adios_helpers.o: $(SRCDIR)/adios_helpers.f90 $(OBJDIR)/adios_helpers_definitions.o $(OBJDIR)/adios_helpers_writers.o
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/sum_kernels.o: $(SRCDIR)/sum_kernels.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/sum_azi_kernels.o: $(SRCDIR)/sum_azi_kernels.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/merge_kernels.o: $(SRCDIR)/merge_kernels.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/merge_azi_kernels.o: $(SRCDIR)/merge_azi_kernels.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/precond_kernels.o: $(SRCDIR)/precond_kernels.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/precond_azi_kernels.o: $(SRCDIR)/precond_azi_kernels.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/precond_vp_vs_kernels_subs.o: $(SRCDIR)/precond_vp_vs_kernels_subs.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/precond_vp_vs_kernels.o: $(SRCDIR)/precond_vp_vs_kernels.f90 $(OBJDIR)/precond_vp_vs_kernels_subs.o $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/convert_adios_subs.o: $(SRCDIR)/convert_adios_subs.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/convert_adios_to_binary.o: $(SRCDIR)/convert_adios_to_binary.f90 $(OBJDIR)/convert_adios_subs.o $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/convert_ascii_to_adios.o: $(SRCDIR)/convert_ascii_to_adios.f90 $(OBJDIR)/convert_adios_subs.o $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/steepDescent.o: $(SRCDIR)/steepDescent.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/conjugateGradient.o: $(SRCDIR)/conjugateGradient.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/lbfgs_subs.o: $(SRCDIR)/lbfgs_subs.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/lbfgs.o: $(SRCDIR)/lbfgs.f90 $(OBJDIR)/lbfgs_subs.o $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/update_model.o: $(SRCDIR)/update_model.f90  $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/model_perturb_ref.o: $(SRCDIR)/model_perturb_ref.f90  $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/blend_model.o: $(SRCDIR)/blend_model.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/gpp_utils.o: $(SRCDIR)/gpp_utils.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/gaussian_perturb_single.o: $(SRCDIR)/gaussian_perturb_single.f90 $(objects) $(OBJDIR)/gpp_utils.o
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/gaussian_perturb_multiple.o: $(SRCDIR)/gaussian_perturb_multiple.f90 $(objects) $(OBJDIR)/gpp_utils.o
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/abs_kernel.o: $(SRCDIR)/abs_kernel.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/compute_azi_params.o: $(SRCDIR)/compute_azi_params.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/apply_source_mask.o: $(SRCDIR)/apply_source_mask.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/compute_vp_vs_hess.o: $(SRCDIR)/compute_vp_vs_hess.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/random_probe_lbfgs.o: $(SRCDIR)/random_probing.f90 $(OBJDIR)/lbfgs.o $(OBJDIR)/compute_vp_vs_hess.o $(OBJDIR)/random.o $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

# ######################
# Link
# ######################

$(BINDIR)/xsum_kernels: $(OBJDIR)/sum_kernels.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xsum_azi_kernels: $(OBJDIR)/sum_azi_kernels.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xmerge_kernels: $(OBJDIR)/merge_kernels.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xmerge_azi_kernels: $(OBJDIR)/merge_azi_kernels.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xprecond_kernels: $(OBJDIR)/precond_kernels.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xprecond_azi_kernels: $(OBJDIR)/precond_azi_kernels.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xprecond_vp_vs_kernels: $(OBJDIR)/precond_vp_vs_kernels.o $(OBJDIR)/precond_vp_vs_kernels_subs.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xbp2binary: $(OBJDIR)/convert_adios_to_binary.o $(OBJDIR)/convert_adios_subs.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xascii2bp: $(OBJDIR)/convert_ascii_to_adios.o $(OBJDIR)/convert_adios_subs.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xsteepDescent: $(OBJDIR)/steepDescent.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xcg_direction: $(OBJDIR)/conjugateGradient.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xlbfgs: $(OBJDIR)/lbfgs.o $(OBJDIR)/lbfgs_subs.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xupdate_model: $(OBJDIR)/update_model.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xmodel_perturb_ref: $(OBJDIR)/model_perturb_ref.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xblend_model: $(OBJDIR)/blend_model.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xgauss_single: $(OBJDIR)/gaussian_perturb_single.o $(objects) $(OBJDIR)/gpp_utils.o
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xgauss_multiple: $(OBJDIR)/gaussian_perturb_multiple.o $(objects) $(OBJDIR)/gpp_utils.o
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xcompute_azi_params: $(OBJDIR)/compute_azi_params.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xabs_kernel: $(OBJDIR)/abs_kernel.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xsrc_mask: $(OBJDIR)/apply_source_mask.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xcompute_vp_vs_hess_subs: $(OBJDIR)/compute_vp_vs_hess_subs.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xcompute_vp_vs_hess: $(OBJDIR)/compute_vp_vs_hess.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xrandom_probe_lbfgs: $(OBJDIR)/random_probe_lbfgs.o $(OBJDIR)/random.o $(OBJDIR)/lbfgs_subs.o $(OBJDIR)/precond_vp_vs_kernels_subs.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc) -I$(OBJDIR) -L$(OBJDIR)

clean:
	rm obj/* bin/*
