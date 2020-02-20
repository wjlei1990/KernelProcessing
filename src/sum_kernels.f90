module sum_kernels_subs

  use global, only : CUSTOM_REAL, myrank, exit_mpi
  implicit none

  contains

  subroutine read_event_file(eventfile, nevent, event_list, weight_list)
    character(len=*), intent(in) :: eventfile
    integer, intent(inout) :: nevent
    character(len=*), dimension(:), allocatable, intent(inout) :: event_list
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout) :: weight_list

    ! local variables
    character(len=500) :: eventname
    real(kind=CUSTOM_REAL) :: weight
    integer :: ios, i
    integer, parameter :: fh = 1001

    open(unit=fh, file=trim(eventfile), status='old', iostat=ios)
    if ( ios /= 0 ) then
      print*, 'ERROR OPENING', trim(eventfile)
      stop
    end if

    if(myrank == 0) write(*, *) "Reading events and weights from: ", trim(eventfile)
    read(fh, *) nevent
    if(myrank == 0) write(*, *) "Number of events: ", nevent
    allocate(event_list(nevent))
    allocate(weight_list(nevent))

    do i=1, nevent
      read(fh, *) eventname, weight
      if (myrank == 0) write(*, '(A, I5, A, I5, A, A, A, ES16.8)') &
        "[", i, "/", nevent, "]", trim(eventname), ": ", weight
      event_list(i) = trim(eventname)
      weight_list(i) = weight
    enddo
    close(fh)
  end subroutine read_event_file

  subroutine get_sys_args(eventfile, kernel_dir, outputfn)
    character(len=*), intent(inout) :: eventfile, kernel_dir, outputfn

    call getarg(1, eventfile)
    call getarg(2, kernel_dir)
    call getarg(3, outputfn)

    if(trim(eventfile) == '' .or. trim(kernel_dir) == '' &
        .or. trim(outputfn) == '') then
      call exit_mpi("Usage: xsum_kernels eventfile kernel_dir outputfn")
    endif

    if(myrank == 0) then
      write(*, *) "Event file(in): ", trim(eventfile)
      write(*, *) "Kernel dir(in): ", trim(kernel_dir)
      write(*, *) "Output file(out, kernel sums): ", trim(outputfn)
    endif
  end subroutine

end module sum_kernels_subs


program sum_kernels
! Adios implementation: Matthieu & Ebru
! Princeton, August 2013
  use mpi
  use adios_read_mod
  use global, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank
  use global, only : init_mpi, exit_mpi
  use AdiosIO
  use sum_kernels_subs

  implicit none

  ! for transverse anisotropy
  ! for isotropy
  !integer,parameter:: NKERNEL=13    !bulk_betah, bulk_betav, bulk_c, eta

  integer:: nevent, ievent, ier
  character(len=500):: eventfile, kernel_dir, outputfn, kernel_file

  character(len=500):: eventname
  character(len=500), dimension(:), allocatable :: event_list
  real(kind=CUSTOM_REAL):: weight
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: weight_list

  integer, parameter :: NKERNELS = 6    !bulk_betah, bulk_betav, bulk_c, eta, rho, hess
  integer, parameter :: hess_idx = 6
  character(len=500), parameter :: kernel_names(NKERNELS) = &
    (/character(len=150) :: "bulk_betah_kl_crust_mantle", "bulk_betav_kl_crust_mantle", &
                            "bulk_c_kl_crust_mantle",     "eta_kl_crust_mantle",        &
                            "rho_kl_crust_mantle",        "hess_kl_crust_mantle"/)

  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: total_kernel
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: kernels

  call init_mpi()

  if (kernel_names(hess_idx) /= "hess_kl_crust_mantle") call exit_mpi("Incorrect hess_idx")

  call get_sys_args(eventfile, kernel_dir, outputfn)
  call read_event_file(eventfile, nevent, event_list, weight_list)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  total_kernel=0.
  do ievent=1, nevent
    eventname = event_list(ievent)
    weight = weight_list(ievent)

    if (myrank==0) write(*,*) 'Reading in kernel for [', ievent, &
        "/", nevent, "]: ", trim(event_list(ievent))
    ! ------------------------------------------------------------
    ! Construct the path of individual kernel files (User-defined)
    ! ------------------------------------------------------------
    kernel_file = trim(kernel_dir)//'/'//trim(eventname)//'.kernels.bp'
    call read_bp_file_real(kernel_file, kernel_names, kernels)

    ! only keep the abs value of the hess kernel
    kernels(:, :, :, :, hess_idx) = abs(kernels(:, :, :, :, hess_idx))

    total_kernel = total_kernel + kernels * weight
  end do

  call write_bp_file(total_kernel, kernel_names, "KERNEL_GROUPS", outputfn)

  if (myrank==0) print*, 'Done summing all the kernels'

  call adios_finalize (myrank, ier)
  call MPI_FINALIZE(ier)

end program sum_kernels
