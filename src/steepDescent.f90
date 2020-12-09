subroutine get_sys_args(grad_file, direction_file)

  use global_var, only : myrank, exit_mpi

  character(len=500), intent(inout):: grad_file, direction_file

  call getarg(1, grad_file)
  call getarg(2, direction_file)

  if(trim(grad_file) == '' .or. trim(direction_file) == '') then
        call exit_mpi('Usage: xcg_direction gradient_file direction_file')
  endif

  if(myrank == 0) then
    write(*, *) "Grad 0 file   (input): ", trim(grad_file)
    write(*, *) "Direct 1 file (output): ", trim(direction_file)
  endif

end subroutine get_sys_args

program main
  use mpi
  use adios_read_mod
  use global_var
  use AdiosIO

  integer, parameter :: NKERNELS = 4
  character(len=500), parameter :: kernel_names(NKERNELS) = &
    (/character(len=500) :: "bulk_betah_kl_crust_mantle",&
                            "bulk_betav_kl_crust_mantle",&
                            "bulk_c_kl_crust_mantle",&
                            "eta_kl_crust_mantle"/)

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: gradient
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: direction

  character(len=500) :: grad_file, direction_file

  integer:: ier

  call init_mpi()

  call get_sys_args(grad_file, direction_file)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  if (myrank == 0) write(*, *) "Start Reading: ", trim(grad_file)
  call read_bp_file_real(grad_file, kernel_names, gradient)

  ! steep descent(just reverse the gradient as the search direction)
  direction = - gradient

  if (myrank == 0) write(*, *) "Start Write: ", trim(direction_file)
  call write_bp_file(direction, kernel_names, "KERNEL_GROUPS", direction_file)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program
