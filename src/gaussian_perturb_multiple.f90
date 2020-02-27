module psf_file_mod
  use global_var, only : CUSTOM_REAL, myrank, exit_mpi
  use gpp_utils, only : NKERNELS, kernel_names, kernels, ibool, x_glob, y_glob, z_glob, &
    add_gaussian_perturb_hv, NGLLX, NGLLY, NGLLZ, NSPEC

  implicit none

  integer :: npts
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: latitudes, longitudes, depths
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: signs, sigma_hs, sigma_vs

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: total_kernels

  contains

  subroutine read_points_file(fn)
    implicit none

    character(len=*), intent(in) :: fn

    integer, parameter :: handle = 1001
    integer :: i, ios

    open(unit=handle, file=trim(fn), status='old', iostat=ios)
    if ( ios /= 0 ) then
      print *, 'Error opening ', trim(fn)
      !call exit_mpi("Error on opening input file")
      stop
    end if

    read(handle, *) npts
    if (myrank == 0) then
      write(*, '(A, I10)') "npts=", npts
    endif

    allocate(latitudes(npts))
    allocate(longitudes(npts))
    allocate(depths(npts))
    allocate(signs(npts))
    allocate(sigma_hs(npts))
    allocate(sigma_vs(npts))

    read(handle, *)
    do i=1, npts
      read(handle, *) latitudes(i), longitudes(i), depths(i), signs(i), &
        sigma_hs(i), sigma_vs(i)
      if(myrank == 0) write(*, '(I5, F10.2, F10.2, F10.2, F8.1, F10.2, F10.2)') &
          i, latitudes(i), longitudes(i), depths(i), &
          signs(i), sigma_hs(i), sigma_vs(i)
    enddo

    close(handle)
  end subroutine read_points_file

  subroutine get_sys_args(point_file, solver_file, output_file)
    character(len=*), intent(in) :: point_file, solver_file, output_file

    call getarg(1, point_file)
    call getarg(2, solver_file)
    call getarg(3, output_file)

    if (trim(point_file) == '' .or. trim(solver_file) == '' .or. &
      trim(output_file) == '') then
      call exit_mpi('Usage: mpirun -n 384 ./bin/xpsf_file point_file solver_file output_file')
    endif

    if (myrank == 0) then
      print*, "point file: ", trim(point_file)
      print*, "solver file: ", trim(solver_file)
      print*, "output file: ", trim(output_file)
    endif

  end subroutine get_sys_args

  subroutine calculate_perturbations(perturb_idx)
    integer, intent(in) :: perturb_idx

    integer :: i

    total_kernels = 0.0
    do i=1, npts
      if(myrank == 0) write(*, '(A, I6, A, I6, A)') "|||>>> Reading point ", i, " / ", npts, " <<<|||"
      kernels = 0.0
      call add_gaussian_perturb_hv(depths(i), latitudes(i), longitudes(i), signs(i), &
        sigma_hs(i), sigma_vs(i), perturb_idx)
      total_kernels = total_kernels + kernels
    enddo

  end subroutine calculate_perturbations

end module psf_file_mod

program main
  use mpi
  use global_var, only : myrank, nprocs, CUSTOM_REAL, init_mpi, exit_mpi
  use gpp_utils, only : read_model_file, kernel_names, write_bp_file
  use psf_file_mod, only : read_points_file, get_sys_args, calculate_perturbations, total_kernels

  implicit none

  character(len=256) :: input_file, solver_file, output_file
  integer, parameter :: perturb_idx = 2

  integer :: ier

  call init_mpi()
  if (myrank == 0) print*, "Rank=", myrank, "/", nprocs
  call get_sys_args(input_file, solver_file, output_file)

  call read_model_file(solver_file)
  call MPI_Barrier(MPI_COMM_WORLD, ier)

  if (myrank == 0) print*, "Reading input file: ", trim(input_file)

  call read_points_file(input_file)
  call MPI_Barrier(MPI_COMM_WORLD, ier)

  if (myrank == 0) print*, "Finish reading input file"

  call calculate_perturbations(perturb_idx)
  call MPI_Barrier(MPI_COMM_WORLD, ier)

  call write_bp_file(total_kernels, kernel_names, "KERNELS_GROUP", output_file)

  call MPI_Barrier(MPI_COMM_WORLD, ier)
  call MPI_FINALIZE(ier)

end program main
