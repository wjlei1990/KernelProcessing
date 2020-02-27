module merge_kernels_sub
  use global_var, only : exit_mpi, myrank
  implicit none

  contains

  subroutine get_sys_args(kernel_smooth_dir, output_file)
    character(len=500), intent(inout) :: kernel_smooth_dir
    character(len=500), intent(inout) :: output_file

    call getarg(1, kernel_smooth_dir)
    call getarg(2, output_file)

    if(trim(kernel_smooth_dir) == "" .or. trim(output_file) == "") then
      call exit_mpi("Usage: xprecond_kernels kernel_smooth_dir output_file")
    endif

    if(myrank == 0) then
      write(*, *) "Kerenel smooth dir(input): ", trim(kernel_smooth_dir)
      write(*, *) "output file: ", trim(output_file)
    endif

  end subroutine get_sys_args

end module merge_kernels_sub

program merge_kernels
  use mpi
  use adios_read_mod
  use AdiosIO
  use global_var, only : NSPEC, NGLLX, NGLLY, NGLLZ, CUSTOM_REAL
  use global_var, only : init_mpi
  use merge_kernels_sub

  implicit none

  integer, parameter:: NKERNELS = 6
  character(len=500), parameter :: kernel_names(NKERNELS) = &
    (/character(len=500) :: "bulk_betah_kl_crust_mantle", "bulk_betav_kl_crust_mantle", &
                            "bulk_c_kl_crust_mantle", "eta_kl_crust_mantle",&
                            "rho_kl_crust_mantle", "hess_kl_crust_mantle"/)

  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: total_kernels

  character(len=500) :: tmp_kernel_name
  character(len=500) :: tmp_kernel_file
  character(len=500) :: tmp_kernel_list(1)
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 1):: tmp_kernel

  character(len=500) :: kernel_smooth_dir
  character(len=500) :: output_file
  integer:: ier, iker

  call init_mpi()

  call get_sys_args(kernel_smooth_dir, output_file)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  do iker=1, NKERNELS
    tmp_kernel_name = trim(kernel_names(iker))
    tmp_kernel_file = trim(adjustl(kernel_smooth_dir))//"/kernels_smooth_"//trim(tmp_kernel_name)//".bp"

    if(myrank==0) then
      write(*, *) "----------- ", iker, "/", NKERNELS, " ------------"
      write(*, *) 'Reading in kernel: ', trim(tmp_kernel_name)
      write(*, *) 'Kernel file: ', trim(tmp_kernel_file)
    endif

    tmp_kernel_list(1) = tmp_kernel_name
    call read_bp_file_real(tmp_kernel_file, tmp_kernel_list, tmp_kernel)

    total_kernels(:, :, :, :, iker) = tmp_kernel(:, :, :, :, 1)
  end do

  if(myrank == 0) write(*, *) "Output file: ", trim(output_file)
  call write_bp_file(total_kernels, kernel_names, "KERNEL_GROUPS", output_file)

  call adios_finalize (myrank, ier)
  call MPI_FINALIZE(ier)

end program merge_kernels
