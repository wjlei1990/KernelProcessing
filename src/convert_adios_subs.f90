module convert_adios_subs
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

end module convert_adios_subs
