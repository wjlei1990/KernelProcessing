!-------------------------------------
! Guassian perturbation for psf
!-------------------------------------
program gaussian_perturb_psf
  use mpi
  use adios_read_mod
  use gpp_utils

  implicit none

  integer, parameter :: perturb_idx = 2
  character(len=500):: input_solver_file, output_file
  integer:: ier

  call init_mpi()
  call get_sys_args(input_solver_file, output_file)

  if( perturb_idx /= betav_kl_idx .or. &
      trim(kernel_names(perturb_idx)) /= "bulk_betav_kl_crust_mantle" ) then
    call exit_mpi("Checking perturb_idx is right (not betav perturbation now)!")
  endif

  if(myrank == 0) write(*, '(A)') "|<-------- Read model files -------->|"
  call read_model_file(input_solver_file)

  if(myrank == 0) write(*, '(A)') "|<-------- Add gaussian perturb -------->|"
  kernels(:,:,:,:, :) = 0.0_CUSTOM_REAL

  ! Currently perturbations are hard-coded in the program
  ! AFAR
  !if(myrank == 0) write(*, '(A)') "|<-------- Afar Plume -------->|"
  !!! call add_gaussian_perturb_sphere(1250.0, 0.0, 28.0, 250.0, perturb_idx)
  !!! call add_gaussian_perturb_hv(2000.0, 0.0, 28.0, 400.0, 400.0, perturb_idx)
  !!! call add_gaussian_perturb_elliptic(300.0, -10.0, 28.0, 500.0, 250.0, perturb_idx)
  !call add_gaussian_perturb_hv(2000.0, -18.0, 20.0, 400.0, 400.0, perturb_idx)

  ! Pacific Plume
  !if(myrank == 0) write(*, '(A)') "|<-------- Pacific Plume -------->|"
  !call add_gaussian_perturb_hv(1250.0, -10.0, -170.0, 500.0, 300.0, perturb_idx)
  !call add_gaussian_perturb_hv(2300.0, -10.0, -170.0, 600.0, 400.0, perturb_idx)

  ! Hawaii
  !if(myrank == 0) write(*, '(A)') "|<-------- Hawaii -------->|"
  !call add_gaussian_perturb_hv(2300.0, 19.5, -155.5, 500.0, 250.0, perturb_idx)

  ! Yellowstone
  !if(myrank == 0) write(*, '(A)') "|<-------- Yellowstone -------->|"
  !call add_gaussian_perturb_hv(300.0, 44.5, -110.4, 100.0, 100.0, perturb_idx)

  !if(myrank == 0) write(*, '(A)') "|<-------- Yellowstone Complex -------->|"
  !call add_gaussian_perturb_hv(200.0, 44.5, -110.4, 100.0, 100.0, perturb_idx)
  !call add_gaussian_perturb_hv(500.0, 44.5, -110.4, 200.0, 200.0, perturb_idx)
  !call add_gaussian_perturb_hv(1000.0, 44.5, -110.4, 250.0, 250.0, perturb_idx)
  !call add_gaussian_perturb_hv(1600.0, 44.5, -110.4, 300.0, 300.0, perturb_idx)
  !call add_gaussian_perturb_hv(2200.0, 44.5, -110.4, 300.0, 300.0, perturb_idx)
  !call add_gaussian_perturb_hv(2800.0, 44.5, -110.4, 300.0, 300.0, perturb_idx)

  ! USA Mid
  !if(myrank == 0) write(*, '(A)') "|<-------- USA Mid -------->|"
  !call add_gaussian_perturb_hv(500.0, 38.0, -98.0, 160.0, 160.0, perturb_idx)

  ! USA East
  !if(myrank == 0) write(*, '(A)') "|<-------- USA East -------->|"
  !call add_gaussian_perturb_hv(1600.0, 36.0, -82.0, 250.0, 250.0, perturb_idx)

  !if(myrank == 0) write(*, '(A)') "|<-------- East-Galapagos -------->|"
  !call add_gaussian_perturb_hv(2500.0, -15.0, -100.0, 300.0, 300.0, perturb_idx)

  if(myrank == 0) write(*, '(A)') "|<-------- Subduction Su -------->|"
  call add_gaussian_perturb_hv(900.0, -6.0, 112.0, 1.0, 200.0, 200.0, perturb_idx)

  if(myrank == 0) write(*, '(A)'), "|<-------- Write gaussian perturb -------->|"
  call write_bp_file(kernels, kernel_names, "KERNELS_GROUP", output_file)

  if(myrank == 0) write(*, '(A)'), "|<-------- Done -------->|"
  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program gaussian_perturb_psf
