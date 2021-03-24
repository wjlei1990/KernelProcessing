! Convert ADIOS file to binary
program convert_adios_to_binary
  use mpi
  use adios_read_mod
  use AdiosIO
  use global_var, only : NGLLX, NGLLY, NGLLZ, NSPEC, myrank, CUSTOM_REAL
  use global_var, only : init_mpi
  use convert_adios_subs

  implicit none

  !integer, parameter :: NKERNELS = 6    !bulk_betah, bulk_betav, bulk_c, eta
  !character(len=500), parameter :: kernel_names(NKERNELS) = &
  !  (/character(len=500) :: "hess_kl_crust_mantle", "bulk_betah_kl_crust_mantle", &
  !                          "bulk_betav_kl_crust_mantle", "bulk_c_kl_crust_mantle", &
  !                          "eta_kl_crust_mantle", "rho_kl_crust_mantle"/)
  !integer, parameter :: NKERNELS = 4
  !character(len=500), parameter :: kernel_names(NKERNELS) = &
  !  (/character(len=500) :: "bulk_betah_kl_crust_mantle", &
  !                          "bulk_betav_kl_crust_mantle", &
  !                          "bulk_c_kl_crust_mantle", &
  !                          "eta_kl_crust_mantle"/)

  character(len=512) :: input_file, varname, output_dir

  integer, parameter :: iregion_code = 1
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC):: vdata = 0.0

  integer:: ier
  character(len=512) :: kname, prname, outputfn, var_path, reg_name
  logical :: is_kernel

  call init_mpi()

  call get_sys_args(varname, input_file, output_dir)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  call get_short_name(varname, kname, is_kernel)

  if (is_kernel) then
    var_path = varname
  else
    write(reg_name,"('reg',i1,'/')") iregion_code
    var_path = trim(reg_name) // trim(varname)
  endif

  if(myrank == 0) write(*, *) "varname, varpath, bin prefix: ", &
    trim(varname), ", ", trim(var_path), ", ", trim(kname)

  call read_bp_file_real(input_file, var_path, vdata)

  ! write binary file out
  call create_name_database(prname, myrank, iregion_code, output_dir)

  outputfn = trim(prname) // trim(kname) // ".bin"
  if (myrank == 0) then
    write(*, *) "[Rank=", myrank, "] Save kernel ", trim(kname), " to: ", trim(outputfn)
  endif
  call write_binary_file(vdata, outputfn)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program convert_adios_to_binary
