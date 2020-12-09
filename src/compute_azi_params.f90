module compute_azi_sub
  use global_var, only : NGLLX, NGLLY, NGLLZ, NSPEC, myrank, CUSTOM_REAL, exit_mpi, myrank

  implicit none

  contains

  subroutine get_sys_args(model_gll_file, kernel_file, output_file)
    character(len=*), intent(inout) :: model_gll_file, kernel_file, output_file

    call getarg(1, model_gll_file)
    call getarg(2, kernel_file)
    call getarg(3, output_file)

    if( model_gll_file == '' .or. kernel_file == '' .or. output_file == '') then
      call exit_mpi("Usage: xcompute_azi_params model_gll kernel output" )
    endif

    if(myrank == 0) then
      write(*, *) "Model GLL file: ", trim(model_gll_file)
      write(*, *) "Input kernel file: ", trim(kernel_file)
      write(*, *) "output file:", trim(output_file)
    endif
  end subroutine get_sys_args

end module compute_azi_sub

program main

  use mpi
  use adios_read_mod
  use AdiosIO
  use global_var, only : NGLLX, NGLLY, NGLLZ, NSPEC, myrank, CUSTOM_REAL
  use global_var, only: init_mpi
  use compute_azi_sub

  implicit none

  character(len=512) :: model_gll_file, kernel_file, output_file

  character(len=512), parameter :: kernel_names(2) = &
    (/character(len=512) :: "Gs_prime_kl_crust_mantle", "Gc_prime_kl_crust_mantle"/)

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 2) :: kernels

  character(len=512), parameter :: mu0_str = "reg1/mu0"
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: mu0

  character(len=512), parameter :: output_names(4) = &
    (/character(len=512) :: "reg1/Gs", "reg1/Gc", "reg1/G0", "reg1/eta"/)
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 4) :: results

  integer :: ier

  call init_mpi()

  call get_sys_args(model_gll_file, kernel_file, output_file)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  !call read_bp_file_real(model_gll_file, mu0_str, mu0)

  call read_bp_file_real(kernel_file, kernel_names, kernels)

  kernels = kernels * 1.0e10

  results(:, :, :, :, 1) = kernels(:, :, :, :, 1)
  results(:, :, :, :, 2) = kernels(:, :, :, :, 2)

  ! G0 = sqrt(Gs**2 + Gc**2)
  results(:, :, :, :, 3) = sqrt(kernels(:, :, :, :, 1)**2 + kernels(:, :, :, :,  2)**2)
  ! eta = 0.5 * atan2(Gs, Gc)
  results(:, :, :, :, 4) = 0.5 * atan2(kernels(:, :, :, :, 1), kernels(:, :, :, :, 2))

  call write_bp_file(results, output_names, "KERNEL_GROUPS", output_file)

end program main

