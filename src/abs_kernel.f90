module abs_kernel_subs

  use global_var, only : CUSTOM_REAL, myrank, exit_mpi
  implicit none

  contains

  subroutine get_sys_args(kernel_name, inputfn, outputfn)
    character(len=*), intent(inout) :: kernel_name, inputfn, outputfn

    call getarg(1, kernel_name)
    call getarg(2, inputfn)
    call getarg(3, outputfn)

    if(trim(inputfn) == '' .or. trim(kernel_name) == '' &
        .or. trim(outputfn) == '') then
      call exit_mpi("Usage: xabs_kernel kernel_name inputfn outputfn")
    endif

    if(myrank == 0) then
      write(*, *) "Kernel name (in): ", trim(kernel_name)
      write(*, *) "Input  file (in): ", trim(inputfn)
      write(*, *) "Output file (out, kernel sums): ", trim(outputfn)
    endif
  end subroutine

end module abs_kernel_subs


program abs_kernel
! Adios implementation: Matthieu & Ebru
! Princeton, August 2013
  use mpi
  use adios_read_mod
  use global_var, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank
  use global_var, only : init_mpi, exit_mpi
  use AdiosIO
  use abs_kernel_subs

  implicit none

  integer:: ier

  character(len=512) :: kernel_name, inputfn, outputfn
  character(len=512) :: kernel_names(1)

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: kernel
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 1) :: kernel_out

  call init_mpi()

  call get_sys_args(kernel_name, inputfn, outputfn)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  kernel = 0.0
  call read_bp_file_real(inputfn, kernel_name, kernel)

  ! only keep the abs value of the hess kernel
  kernel_out(:, :, :, :, 1) = abs(kernel)

  kernel_names(1) = kernel_name
  call write_bp_file(kernel_out, kernel_names, "KERNEL_GROUPS", outputfn)

  call adios_finalize (myrank, ier)
  call MPI_FINALIZE(ier)

end program abs_kernel
