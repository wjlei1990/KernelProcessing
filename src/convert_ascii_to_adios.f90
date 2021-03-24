program main

  use mpi
  use AdiosIO
  use global_var, only : CUSTOM_REAL, init_mpi, exit_mpi, myrank
  use global_var, only : NGLLX, NGLLY, NGLLZ, NSPEC
  use convert_adios_subs, only : create_name_database

  integer :: iregion_code = 1
  character(len=512) :: inputdir, outputfile
  character(len=512) :: prname, filename

  integer, parameter :: NKERNELS = 4
  character(len=512), parameter :: kernel_names(NKERNELS) = &
        (/character(len=512) :: "bulk_c_kl_crust_mantle", &
                                "bulk_betav_kl_crust_mantle", &
                                "bulk_betah_kl_crust_mantle", &
                                "eta_kl_crust_mantle"/)

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS) :: kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: ldata

  integer :: i
  integer :: ioint = 111
  logical :: file_exist

  call init_mpi()

  call getarg(1, inputdir)
  call getarg(2, outputfile)

  if(trim(inputdir) == "" .or. trim(outputfile) == "") then
    call exit_mpi("Usage: xascii2bp inputdir outputfile")
  endif

  call create_name_database(prname, myrank, iregion_code, inputdir)

  do i =1, NKERNELS
    filename = trim(prname) // trim(kernel_names(i)) // ".dat"
    print*, "rank: ", myrank, " | ", trim(filename)
    inquire( file=trim(filename), exist=file_exist)
    if(.not. file_exist) call exit_mpi("input file not found")

    open(unit=ioint, file=filename, status='old',    &
         access='sequential', form='formatted', action='read')
    read(ioint, *) ldata
    kernels(:, :, :, :, i) = ldata
    close(unit=ioint)
  enddo

  call write_bp_file(kernels, kernel_names, "kernel_groups", outputfile)

  call adios_finalize(myrank, ier)
  call mpi_finalize(ier)

end program main
