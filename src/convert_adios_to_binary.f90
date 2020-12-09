! Ebru1: The choice of THRESHOLD value is somewhat subjective. It is not trivial to set it like the 20% of max value
! which may be OK for smaller scale studies but global scale needs a few trial&error to adjust this parameter for
! every iteration. Needs some more investigation..

! Ebru2: I find the preconditioner behave better after changing the order of smoothing and preconditioning in
! post-processing upon the suggestion by Ryan & Yanhua.
! However, I am still not convinced by Ryan's latest suggestion that preconditioner should be smoothed more than the
! gradients of other parameters that the preconditioner to be applied. I currently smooth the preconditioner and
! the other gradients in the same way.

module convert_adios_to_bin_sub
  use mpi
  use global_var, only : CUSTOM_REAL, exit_mpi, myrank
  implicit none

  contains

  subroutine get_sys_args(varname, input_file, output_dir)
    character(len=*), intent(inout) :: varname, input_file, output_dir

    call getarg(1, varname)
    call getarg(2, input_file)
    call getarg(3, output_dir)

    if(varname == '' .or. input_file == '' .or. output_dir == '') then
      call exit_mpi("Usage: xprecond_kernels param input_kernel output_kernel")
    endif

    if(myrank == 0) then
      write(*, *) "Variable: ", trim(varname)
      write(*, *) "Input ADIOS kernel: ", trim(input_file)
      write(*, *) "Output dir: ", trim(output_dir)
    endif
  end subroutine get_sys_args

  subroutine write_binary_file(values, filename)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: values
    character(len=*), intent(in) :: filename
    integer :: IOUT = 321

    open(unit=IOUT, file=filename, status='unknown', form='unformatted', action='write')
    write(IOUT) values
    close(IOUT)
  end subroutine write_binary_file

  ! create the name of the database for the mesher and the solver
  subroutine create_name_database(prname, iproc, iregion_code, LOCAL_PATH)
    character(len=*), intent(inout) :: prname
    integer, intent(in) :: iproc, iregion_code
    character(len=*), intent(in) :: LOCAL_PATH

    character(len=500) :: procname

    ! create the name for the database of the current slide and region
    write(procname, "('/proc',i6.6,'_reg',i1,'_')") iproc,iregion_code
    ! create full name with path
    prname = trim(LOCAL_PATH) // trim(procname)
  end subroutine create_name_database

  subroutine get_short_name(full_name, short_name, is_kernel)
    character(len=*), intent(in) :: full_name
    character(len=*), intent(out) :: short_name
    logical :: is_kernel

    integer :: slen

    slen = len_trim(full_name)

    if ( slen > 16 .and. full_name((slen-15):slen) == 'kl_crust_mantle' ) then
      short_name = full_name(1:(slen-16)) // "_kernel"
      is_kernel = .true.
    else
      short_name = full_name
      is_kernel = .false.
    endif

  end subroutine get_short_name

end module convert_adios_to_bin_sub

! Convert ADIOS file to binary
program convert_adios_to_binary
  use mpi
  use adios_read_mod
  use AdiosIO
  use global_var, only : NGLLX, NGLLY, NGLLZ, NSPEC, myrank, CUSTOM_REAL
  use global_var, only : init_mpi
  use convert_adios_to_bin_sub

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
