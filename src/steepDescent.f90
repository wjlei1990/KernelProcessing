subroutine get_sys_args(grad_file, precond_file, direction_file)

  use global_var, only : myrank, exit_mpi

  character(len=512), intent(inout):: grad_file, precond_file, direction_file

  call getarg(1, grad_file)
  call getarg(2, precond_file)
  call getarg(3, direction_file)

  if(trim(grad_file) == '' .or. trim(precond_file) == '' .or. trim(direction_file) == '') then
        call exit_mpi('Usage: xcg_direction gradient_file precond_file direction_file')
  endif

  if(myrank == 0) then
    write(*, *) "Gradient file   (input): ", trim(grad_file)
    write(*, *) "Precond file    (input): ", trim(precond_file)
    write(*, *) "Direction file (output): ", trim(direction_file)
  endif

end subroutine get_sys_args

program main
  use mpi
  use adios_read_mod
  use global_var
  use AdiosIO

  integer, parameter :: NKERNELS = 4
  character(len=512), parameter :: kernel_names(NKERNELS) = &
    (/character(len=512) :: "bulk_c_kl_crust_mantle",&
                            "bulk_betav_kl_crust_mantle",&
                            "bulk_betah_kl_crust_mantle",&
                            "eta_kl_crust_mantle"/)

  character(len=512), parameter :: precond_names(NKERNELS) = &
    (/character(len=512) :: "precond_bulk_c",&
                            "precond_bulk_betav",&
                            "precond_bulk_betah",&
                            "precond_eta"/)

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: gradient
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: precond
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: direction

  character(len=512) :: grad_file, precond_file, direction_file

  integer:: ier

  call init_mpi()

  call get_sys_args(grad_file, precond_file, direction_file)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  if (myrank == 0) write(*, *) "Reading Gradient: ", trim(grad_file)
  call read_bp_file_real(grad_file, kernel_names, gradient)

  if (myrank == 0) write(*, *) "Reading Preconditioner: ", trim(grad_file)
  call read_bp_file_real(precond_file, precond_names, precond)

  ! steep descent method with preconditioner applied
  direction = - precond * gradient

  if (myrank == 0) write(*, *) "Write search direction: ", trim(direction_file)
  call write_bp_file(direction, kernel_names, "KERNEL_GROUPS", direction_file)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program
