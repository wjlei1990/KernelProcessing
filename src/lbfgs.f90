module lbfgs_subs

  use mpi
  use global, only : CUSTOM_REAL, myrank, exit_mpi, NGLLX, NGLLY, NGLLZ, &
                        NSPEC
  use global, only : Parallel_ComputeInnerProduct, Parallel_ComputeL2normSquare
  use AdiosIO, only : read_bp_file_real
  implicit none

  real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), allocatable :: stored_gradients
  ! According to Numerical Optimization, page 177:
  !     sk = x(k+1) - x(k), which is the model change
  !     yk = g(k+1) - g(k), which is the gradient change
  real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), allocatable :: yks, sks

  integer, parameter :: NKERNELS = 4
  character(len=500), parameter :: kernel_names(NKERNELS) = &
        (/character(len=150) :: "bulk_c_kl_crust_mantle", &
                                "bulk_betav_kl_crust_mantle", &
                                "bulk_betah_kl_crust_mantle", &
                                "eta_kl_crust_mantle"/)

  contains

  subroutine get_sys_args(input_path_file, solver_file, outputfn)
    character(len=*), intent(inout) :: input_path_file, solver_file, outputfn

    if(myrank == 0) print*, '|<============= Get System Args =============>|'
    call getarg(1, input_path_file)
    call getarg(2, solver_file)
    call getarg(3, outputfn)

    if(trim(input_path_file) == '' .or. trim(solver_file) == '' &
        .or. trim(outputfn) == '') then
      call exit_mpi("Usage: ./xlbfgs input_path_file solver_file outputfn")
    endif

    if(myrank == 0) then
      print*, "Input path file: ", trim(input_path_file)
      print*, "solver file:", trim(solver_file)
      print*, "Output file: ", trim(outputfn)
    endif
  end subroutine get_sys_args

  subroutine parse_input_path(fn, niter, gradient_files, model_change_files)
    character(len=*), intent(in) :: fn
    integer, intent(inout) :: niter
    character(len=500), dimension(:), allocatable, intent(inout) :: gradient_files, &
                                                                    model_change_files

    integer :: fh = 1001
    integer :: i, ier

    if(myrank == 0) print*, '|<============= Parsing Input Path File =============>|'

    open(fh, file=trim(fn), status='old')

    read(fh, *) niter
    if(myrank == 0) print*, "niter: ", niter
    ! The gradient files is for each iteration, gradient of M(n), ..., M(n-k)
    ! So the size is (k + 1)
    allocate(gradient_files(niter+1))
    ! The model change files is the model change between two continues
    ! So the size is (k)
    allocate(model_change_files(niter))

    ! read the last k iteration
    do i=1, niter
      read(fh, '(A)') gradient_files(i)
      read(fh, '(A)') model_change_files(i)
      if(myrank == 0) then
        write(*, '(A, I3, A, A)'), "[iter=", niter, "] grad file:   ", trim(gradient_files(i))
        write(*, '(A, I3, A, A)'), "[iter=", niter, "] dkernel file:", trim(model_change_files(i))
      endif
    enddo

    ! read the gradient from current iteration
    read(fh, '(A)') gradient_files(niter+1)
    if(myrank == 0) then
      print*, "[Current iter] grad file: ", trim(gradient_files(niter+1))
    endif

    close(fh)

    call MPI_Barrier(MPI_COMM_WORLD, ier)
  end subroutine parse_input_path

  subroutine read_all_bp_files(niter, gradient_files, model_change_files)
    integer, intent(in) :: niter
    character(len=*), dimension(:), intent(in) :: gradient_files, model_change_files

    integer :: i

    if(myrank == 0) print*, '|<============= Reading Previous BP files =============>|'
    allocate(stored_gradients(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS, niter+1))
    allocate(sks(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS, niter))
    allocate(yks(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS, niter))

    do i=1, niter+1
      if(myrank == 0) write(*, '(A, I2, A, A)') &
        "Reading [iter=", i, "] Gradient: ", trim(gradient_files(i))
      call read_bp_file_real(gradient_files(i), kernel_names, &
                             stored_gradients(:, :, :, :, :, i))
    enddo

    do i=1, niter
      yks(:, :, :, :, :, i) = &
        stored_gradients(:, :, :, :, :, i+1) - stored_gradients(:, :, :, :, :, i)
    enddo

    do i=1, niter
      if(myrank == 0) write(*, '(A, I2, A, A)') &
        "Reading [iter=", i, "] dkernel: ", trim(model_change_files(i))
      call read_bp_file_real(model_change_files(i), kernel_names, &
                             sks(:, :, :, :, :, i))
    enddo

  end subroutine read_all_bp_files

  subroutine calculate_LBFGS_direction(niter, jacobian, direction)
    integer, intent(in) :: niter
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: jacobian
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: direction

    integer :: i
    real(kind=CUSTOM_REAL) :: tmp, rhok, beta, norm_y
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: pk_store, ak_store

    if(myrank == 0) print*, '|<============= Calculating LBFGS Direction =============>|'
    allocate(pk_store(niter))
    allocate(ak_store(niter))

    direction = stored_gradients(:, :, :, :, :, niter+1)

    ! first round
    do i=niter, 1, -1
      call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, i), yks(:, :, :, :, :, i), &
                                        NKERNELS, jacobian, tmp)
      pk_store(i) = 1.0 / tmp

      call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, i), direction, &
                                        NKERNELS, jacobian, tmp)
      ak_store(i) = pk_store(i) * tmp

      direction = direction - ak_store(i) * yks(:, :, :, :, :, i)
    enddo

    ! Precondition
    call Parallel_ComputeL2normSquare(yks(:, :, :, :, :, niter), NKERNELS, &
                                          jacobian, norm_y)
    rhok = 1.0 / (pk_store(niter) * norm_y)
    direction = rhok * direction

    ! second round
    do i=1, niter, 1
      call Parallel_ComputeInnerProduct(yks(:, :, :, :, :, i),  direction, &
                                        NKERNELS, jacobian, tmp)
      beta = pk_store(i) * tmp

      direction = direction + (ak_store(i) - beta) * sks(:, :, :, :, :, i)
    enddo

    ! reverse to get descent direction
    direction = -1.0 * direction

  end subroutine calculate_LBFGS_direction
end module lbfgs_subs

program main
  use mpi
  use adios_read_mod
  use lbfgs_subs
  use AdiosIO, only : write_bp_file, calculate_jacobian_matrix
  use global, only : init_mpi
  implicit none

  ! Number of previous iteration used in L-BFGS
  integer :: niter
  character(len=500) :: input_path_file, solver_file, outputfn
  character(len=500), dimension(:), allocatable :: gradient_files, model_change_files

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: jacobian
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: direction

  integer :: ier

  call init_mpi()

  if(myrank == 0) print*, "|<---- Get System Args ---->|"
  call get_sys_args(input_path_file, solver_file, outputfn)

  if(myrank == 0) print*, "|<---- Parse Input Path File ---->|"
  call parse_input_path(input_path_file, niter, gradient_files, model_change_files)

  call read_all_bp_files(niter, gradient_files, model_change_files)

  if(myrank == 0) print*, "|<---- Calculate Jacobian ---->|"
  call calculate_jacobian_matrix(solver_file, jacobian)

  if(myrank == 0) print*, "|<---- L-BFGS Direction ---->|"
  call calculate_LBFGS_direction(niter, jacobian, direction)

  call write_bp_file(direction, kernel_names, "KERNELS_GROUP", outputfn)
  if(myrank == 0) print*, "LBFGS direction saved: ", trim(outputfn)

  call adios_finalize(myrank, ier)
  call MPI_finalize(ier)

end program main
