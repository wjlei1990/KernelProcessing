SRCDIR=src
OBJDIR=obj
BINDIR=bin
FC=gfortran
MPIFC=mpif90
#MPIFC=ftn
FCFLAGS=-O3 -Wall -J $(OBJDIR) -I $(OBJDIR)

#adios_link=$(shell adios_config -lf)
#adios_inc=$(shell adios_config -cf)

adios_link=$(shell /ccs/home/ccui/adios-gcc/adios_config -lf)
adios_inc=$(shell /ccs/home/ccui/adios-gcc/adios_config -cf)

objects= $(OBJDIR)/adios_helpers_definitions.o $(OBJDIR)/adios_helpers_writers.o $(OBJDIR)/adios_helpers.o $(OBJDIR)/gll_library.o $(OBJDIR)/global_var.o $(OBJDIR)/AdiosIO.o

all: $(BINDIR)/xsteepDescent $(BINDIR)/xcg_direction $(BINDIR)/xlbfgs $(BINDIR)/xsum_kernels $(BINDIR)/xprecond_kernels $(BINDIR)/xmerge_kernels $(BINDIR)/xupdate_model $(BINDIR)/xmodel_perturb_ref $(BINDIR)/xblend_model $(BINDIR)/xgauss_single $(BINDIR)/xgauss_multiple

$(OBJDIR)/global_var.o: $(SRCDIR)/global_var.f90 $(OBJDIR)/gll_library.o
	$(MPIFC) $(FCFLAGS) -c $< -o $@

$(OBJDIR)/AdiosIO.o: $(SRCDIR)/AdiosIO.f90 $(OBJDIR)/adios_helpers.o
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/adios_helpers.o: $(SRCDIR)/adios_helpers.f90 $(OBJDIR)/adios_helpers_definitions.o $(OBJDIR)/adios_helpers_writers.o
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/sum_kernels.o: $(SRCDIR)/sum_kernels.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/merge_kernels.o: $(SRCDIR)/merge_kernels.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/precond_kernels.o: $(SRCDIR)/precond_kernels.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_link) $(adios_inc)

$(OBJDIR)/steepDescent.o: $(SRCDIR)/steepDescent.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/conjugateGradient.o: $(SRCDIR)/conjugateGradient.f90 $(objects)
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(OBJDIR)/lbfgs.o: $(SRCDIR)/lbfgs.f90 $(objects)
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

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(MPIFC) $(FCFLAGS) -c $< -o $@ $(adios_inc)

$(BINDIR)/xsum_kernels: $(OBJDIR)/sum_kernels.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xmerge_kernels: $(OBJDIR)/merge_kernels.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xprecond_kernels: $(OBJDIR)/precond_kernels.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xsteepDescent: $(OBJDIR)/steepDescent.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xcg_direction: $(OBJDIR)/conjugateGradient.o $(objects)
	$(MPIFC) $(FCFLAGS) -o $@ $^ $(adios_link) $(adios_inc)

$(BINDIR)/xlbfgs: $(OBJDIR)/lbfgs.o $(objects)
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

clean:
	rm obj/* bin/*
