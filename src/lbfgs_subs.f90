module lbfgs_subs

  use mpi
  use global_var, only : CUSTOM_REAL, myrank, exit_mpi, NGLLX, NGLLY, NGLLZ, &
                     NSPEC
  use global_var, only : Parallel_ComputeInnerProduct, Parallel_ComputeL2normSquare
  use AdiosIO, only : read_bp_file_real
  implicit none

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
    character(len=512), dimension(:), allocatable, intent(inout) :: gradient_files, &
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
        write(*, '(A, I3, A, A)'), "[iter=", i, "] grad file:   ", trim(gradient_files(i))
        write(*, '(A, I3, A, A)'), "[iter=", i, "] dkernel file:", trim(model_change_files(i))
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

  subroutine read_all_bp_files(niter, nkernels, gradient_files, model_change_files, &
                               kernel_names, gradient, yks, sks)

    integer, intent(in) :: niter
    character(len=*), dimension(:), intent(in) :: gradient_files, model_change_files
    integer, intent(in) :: nkernels
    character(len=*), dimension(:), intent(in) :: kernel_names

    ! (NGLLX, NGLLY, NGLLZ, NSPEC, nkernels)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :) :: gradient

    ! (NGLLX, NGLLY, NGLLZ, NSPEC, nkernels, niter)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), allocatable :: yks, sks

    ! local variable
    ! (NGLLX, NGLLY, NGLLZ, NSPEC, nkernels, niter)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), allocatable :: stored_gradients

    integer :: i

    if(myrank == 0) print*, '|<============= Reading Previous BP files =============>|'
    allocate(stored_gradients(NGLLX, NGLLY, NGLLZ, NSPEC, nkernels, niter+1))

    ! Read gradient
    do i=1, niter+1
      if(myrank == 0) write(*, '(A, I2, A, A)') &
        "Reading [iter=", i, "] Gradient: ", trim(gradient_files(i))
      call read_bp_file_real(gradient_files(i), kernel_names, &
                             stored_gradients(:, :, :, :, :, i))
    enddo

    ! store the last iteration gradient
    gradient = stored_gradients(:, :, :, :, :, niter+1)

    ! Calculate the gradient change
    do i=1, niter
      yks(:, :, :, :, :, i) = &
        stored_gradients(:, :, :, :, :, i+1) - stored_gradients(:, :, :, :, :, i)
    enddo

    ! Read Model Change
    do i=1, niter
      if(myrank == 0) write(*, '(A, I2, A, A)') &
        "Reading [iter=", i, "] dkernel: ", trim(model_change_files(i))
      call read_bp_file_real(model_change_files(i), kernel_names, &
                             sks(:, :, :, :, :, i))
    enddo

    deallocate(stored_gradients)

  end subroutine read_all_bp_files

  subroutine calculate_LBFGS_direction(niter, nkernels, jacobian, gradient, &
                                       yks, sks, direction)
    integer, intent(in) :: niter
    integer, intent(in) :: nkernels
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: jacobian
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: gradient
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), intent(in) :: yks, sks
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: direction

    integer :: i
    real(kind=CUSTOM_REAL) :: tmp, rhok, beta, norm_y
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: pk_store, ak_store

    if(myrank == 0) print*, '|<============= Calculating LBFGS Direction =============>|'
    allocate(pk_store(niter))
    allocate(ak_store(niter))

    direction = gradient

    ! first round
    do i=niter, 1, -1
      call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, i), &
                                        yks(:, :, :, :, :, i), &
                                        nkernels, jacobian, tmp)
      pk_store(i) = 1.0 / tmp

      call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, i), &
                                        direction, &
                                        nkernels, jacobian, tmp)
      ak_store(i) = pk_store(i) * tmp

      if(myrank == 0) write(*, '(A, I2, A, i2, A, ES18.8, A, ES18.8)') &
        "First Loop[iter=", i, "/", niter, "] pk=", pk_store(i), ", ak=", ak_store(i)

      direction = direction - ak_store(i) * yks(:, :, :, :, :, i)
    enddo

    ! Precondition
    call Parallel_ComputeL2normSquare(yks(:, :, :, :, :, niter), nkernels, &
                                      jacobian, norm_y)
    !rhok = 1.0 / (pk_store(niter) * norm_y)
    rhok = 1.0
    if(myrank == 0) write(*, '(A, ES18.8)') "Precondition coef: ", rhok
    direction = rhok * direction

    ! second round
    do i=1, niter, 1
      call Parallel_ComputeInnerProduct(yks(:, :, :, :, :, i),  direction, &
                                        nkernels, jacobian, tmp)
      beta = pk_store(i) * tmp

      if(myrank == 0) write(*, '(A, I3, A, I3, A, ES18.8)') &
        "Second Loop[iter=", i, "/", niter, "] beta=", beta
      direction = direction + (ak_store(i) - beta) * sks(:, :, :, :, :, i)
    enddo

    ! reverse to get descent direction
    direction = -1.0 * direction

  end subroutine calculate_LBFGS_direction

end module lbfgs_subs
