subroutine get_sys_args(ref_model_file, new_model_file, outputfn)
  use global, only : myrank
  use global, only : exit_mpi
  implicit none
  character(len=500), intent(in) :: ref_model_file, new_model_file, outputfn

  call getarg(1, ref_model_file)
  call getarg(2, new_model_file)
  call getarg(3, outputfn)

  if(trim(ref_model_file) == '' .or. trim(new_model_file) == '' .or. &
      trim(outputfn) == "") then
    call exit_mpi('Usage: xmodel_perturbs ref_model_file new_model_file outputfn')
  endif

  if (myrank == 0) then
    print*, "ref model file: ", trim(ref_model_file)
    print*, "new model file: ", trim(new_model_file)
    print*, "output file: ", trim(outputfn)
  endif

end subroutine get_sys_args

program main

  use mpi
  use global, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank, init_mpi
  use AdiosIO
  implicit none

  character(len=500) :: ref_model_file, new_model_file, outputfn

  integer, parameter :: nvars = 6
  character(len=500), dimension(nvars), parameter :: model_names = &
    (/character(len=500) :: "reg1/vpv", "reg1/vph", "reg1/vsv", &
                            "reg1/vsh", "reg1/eta", "reg1/rho"/)

  character(len=500), dimension(nvars), parameter :: perturb_names = &
    (/character(len=500) :: "reg1/dvpvvpv", "reg1/dvphvph", "reg1/dvsvvsv", &
                            "reg1/dvshvsh", "reg1/detaeta", "reg1/drhorho"/)

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, nvars) :: ref_model, &
                                                                          new_model, &
                                                                          perturb_model
  character(len=500) :: group_name = "KERNELS_GROUP"

  integer :: ier

  call init_mpi()
  !call MPI_INIT(ier)
  !call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)
  !call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ier)

  if(myrank == 0) print*, "mpi done"

  call get_sys_args(ref_model_file, new_model_file, outputfn)

  if(myrank == 0) print*, "Read ref model file"
  call read_bp_file_real(ref_model_file, model_names, ref_model)
  if(myrank == 0) print*, "Read new model file"
  call read_bp_file_real(new_model_file, model_names, new_model)
  if(myrank == 0) print*, "Done reading"

  perturb_model = log(new_model / ref_model)

  if(myrank == 0) print*, "Write perturb model file"
  call write_bp_file(perturb_model, perturb_names, group_name, outputfn)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

  if(myrank == 0) print*, "Job finished"

end program main
