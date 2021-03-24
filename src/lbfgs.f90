program main
  use mpi
  use adios_read_mod
  use lbfgs_subs
  use AdiosIO, only : write_bp_file, calculate_jacobian_matrix
  use global_var, only : init_mpi
  implicit none

  integer, parameter :: NKERNELS = 4
  character(len=512), parameter :: kernel_names(NKERNELS) = &
        (/character(len=512) :: "bulk_c_kl_crust_mantle", &
                                "bulk_betav_kl_crust_mantle", &
                                "bulk_betah_kl_crust_mantle", &
                                "eta_kl_crust_mantle"/)

  ! Number of previous iteration used in L-BFGS
  integer :: niter
  character(len=512) :: input_path_file, solver_file, outputfn
  character(len=512), dimension(:), allocatable :: gradient_files, model_change_files

  ! jacobian related to the geometry of mesh
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: jacobian
  ! precond kernels (default = 1.0, no preconditioner applied)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: precond = 1.0
  ! most recent gradient (this iteration's gradient)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: gradient
  ! this iterations' search direction
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: direction

  ! According to Numerical Optimization, page 177:
  !     sk = x(k+1) - x(k), which is the model change
  !     yk = g(k+1) - g(k), which is the gradient change
  real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), allocatable :: yks, sks

  integer :: ier

  call init_mpi()

  if(myrank == 0) print*, "|<---- Get System Args ---->|"
  call get_sys_args(input_path_file, solver_file, outputfn)

  if(myrank == 0) print*, "|<---- Parse Input Path File ---->|"
  call parse_input_path(input_path_file, niter, gradient_files, model_change_files)

  allocate(sks(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS, niter))
  allocate(yks(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS, niter))

  call read_all_bp_files(niter, NKERNELS, gradient_files, model_change_files, &
                         kernel_names, gradient, yks, sks)

  if(myrank == 0) print*, "|<---- Calculate Jacobian ---->|"
  call calculate_jacobian_matrix(solver_file, jacobian)

  if(myrank == 0) print*, "|<---- L-BFGS Direction ---->|"
  call calculate_LBFGS_direction(niter, NKERNELS, jacobian, gradient, precond, yks, sks, direction)

  call write_bp_file(direction, kernel_names, "KERNELS_GROUP", outputfn)
  if(myrank == 0) print*, "LBFGS direction saved: ", trim(outputfn)

  call adios_finalize(myrank, ier)
  call MPI_finalize(ier)

end program main
