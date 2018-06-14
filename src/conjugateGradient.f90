! -------------------------------------------------------------------
! A parallel implementation of Conjugate Gradient method
! -------------------------------------------------------------------
module ConjugateGradient

  use mpi
  use global, only : CUSTOM_REAL, CUSTOM_MPI_TYPE, myrank, &
    Parallel_ComputeL2normSquare, Parallel_ComputeInnerProduct, &
    sum_all_all_cr

  implicit none

  contains

  subroutine get_beta_old(gradient_0, gradient_1, beta)
    ! This is the old version from Ebru, which is WRONG
    ! The inner product in GLL space should considering the weighting and
    ! jacobian matrix.
    ! This function is just for the use of benchmark.
    real(kind=CUSTOM_REAL),dimension(:, :, :, :, :), intent(in):: gradient_0, gradient_1
    real(kind=CUSTOM_REAL), intent(inout) :: beta

    ! local variables
    real(kind=CUSTOM_REAL),dimension(:), allocatable:: beta_upper_all, beta_down_all
    real(kind=CUSTOM_REAL)::beta_upper, beta_down, beta_upper_all_tmp, beta_down_all_tmp
    integer :: ier, iker
    integer :: nkernels

    nkernels = size(gradient_0, 5)
    if (myrank == 0) write(*, *) "Number of kerenels: ", nkernels

    allocate(beta_upper_all(nkernels))
    allocate(beta_down_all(nkernels))

    do iker = 1, nkernels
      beta_upper = sum( gradient_1(:, :, :, :, iker) * &
        (gradient_1(:, :, :, :, iker) - gradient_0(:, :, :, :, iker)))

      beta_down = sum(gradient_0(:, :, :, :, iker) * gradient_0(:, :, :, :, iker))

      call mpi_barrier(MPI_COMM_WORLD, ier)
      call sum_all_all_cr(beta_upper, beta_upper_all_tmp)
      call sum_all_all_cr(beta_down, beta_down_all_tmp)
      beta_upper_all(iker) = beta_upper_all_tmp
      beta_down_all(iker) = beta_down_all_tmp
    end do

    beta = sum(beta_upper_all) / sum(beta_down_all)

    if (myrank == 0) then
      write(*, *) "beta(old version, wrong!): ", beta
    endif
  end subroutine get_beta_old

  subroutine get_beta(gradient_0, gradient_1, jacobian, beta)
    real(kind=CUSTOM_REAL),dimension(:, :, :, :, :), intent(in):: gradient_0, gradient_1
    real(kind=CUSTOM_REAL),dimension(:, :, :, :), intent(in):: jacobian
    real(kind=CUSTOM_REAL), intent(inout) :: beta

    real(kind=CUSTOM_REAL) :: beta_up, beta_down
    real(kind=CUSTOM_REAL) :: orth, orth_up, orth_down
    integer :: nkernels

    nkernels = size(gradient_0, 5)
    if (myrank == 0) write(*, *) "Number of kerenels: ", nkernels

    call Parallel_ComputeInnerProduct(gradient_1, gradient_1 - gradient_0, &
                                      nkernels, jacobian, beta_up)
    call Parallel_ComputeL2normSquare(gradient_0, nkernels, jacobian, beta_down)

    beta = beta_up / beta_down
    write(*, *) "Initial beta: ", beta
    ! Restart condition 1: beta must be >= 0
    if (beta < 0.0) then
      if (myrank == 0) write(*, *) "Beta change by restart condition(beta>=0): ", beta, "-> 0.0"
      beta = 0.0
      return
    endif

    ! Restart condition 2: g0 and g1 must be orthogonal enough
    ! Eq(5.52) on Page 125 on Numerical Optimization
    call Parallel_ComputeInnerProduct(gradient_1, gradient_0, nkernels, jacobian, orth_up)
    call Parallel_ComputeL2normSquare(gradient_1, nkernels, jacobian, orth_down)
    orth = abs(orth_up / orth_down)
    if(myrank == 0) write(*, *) "Orthogonal coef: ", orth
    if (orth > 0.1) then
      if (myrank == 0) write(*, *) "Beta change by restart condition(Orth<0.1): ", beta, "-> 0.0"
      beta = 0.0
    endif
  end subroutine get_beta

  subroutine compute_search_direction(gradient_0, gradient_1, direction_0, jacobian, &
                                      direction_1)
    ! Dimension of gradient and direction would be
    !     (NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: gradient_0, gradient_1
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: direction_0
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in)    :: jacobian
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: direction_1

    real(kind=CUSTOM_REAL) :: beta

    call get_beta_old(gradient_0, gradient_1, beta)
    !call get_beta(gradient_0, gradient_1, jacobian, beta)

    if(myrank == 0) write(*, *) "Final beta used: ", beta

    direction_1 = -gradient_1 + beta * direction_0

  end subroutine
end module ConjugateGradient

subroutine get_sys_args(grad_0_file, grad_1_file, &
                        direction_0_file, direction_1_file, solver_file)

  use global, only : myrank, exit_mpi

  character(len=*), intent(inout):: grad_0_file, grad_1_file
  character(len=*), intent(inout):: direction_0_file, direction_1_file
  character(len=*), intent(inout):: solver_file

  call getarg(1, grad_0_file)
  call getarg(2, grad_1_file)
  call getarg(3, direction_0_file)
  call getarg(4, solver_file)
  call getarg(5, direction_1_file)

  if(trim(grad_0_file) == '' .or. trim(grad_1_file) == '' &
     .or. trim(direction_0_file) == '' .or. trim(direction_1_file) == '' &
     .or. trim(solver_file) == '') then
        call exit_mpi('Usage: xcg_direction gradient_0_file '//&
          'gradient_1_file direction_0_file solver_bp_file outputfn')
  endif

  if(myrank == 0) then
    write(*, *) "Grad 0 file   (input): ", trim(grad_0_file)
    write(*, *) "Grad 1 file   (input): ", trim(grad_1_file)
    write(*, *) "Direct 0 file (input): ", trim(direction_0_file)
    write(*, *) "solver bp file(input): ", trim(solver_file)
    write(*, *) "Direct 1 file (output): ", trim(direction_1_file)
  endif

end subroutine get_sys_args


program main
  use mpi
  use global
  use AdiosIO
  use ConjugateGradient

  integer, parameter:: NKERNEL=4
  character(len=500), dimension(NKERNEL):: kernel_names = &
    (/character(len=500) :: "bulk_betah_kl_crust_mantle",&
                            "bulk_betav_kl_crust_mantle",&
                            "bulk_c_kl_crust_mantle",&
                            "eta_kl_crust_mantle"/)

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL):: gradient_0, gradient_1
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL):: direction_0, direction_1
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: jacobian

  character(len=500) :: solver_file
  character(len=500) :: grad_0_file, grad_1_file
  character(len=500) :: direction_0_file
  character(len=500) :: direction_1_file ! outputfn

  integer:: ier

  call init_mpi()

  call get_sys_args(grad_0_file, grad_1_file, &
                    direction_0_file, direction_1_file, solver_file)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  if (myrank == 0) write(*, *) "|<----- Start Reading ----->|"
  call read_bp_file_real(grad_0_file, kernel_names, gradient_0)
  call read_bp_file_real(grad_1_file, kernel_names, gradient_1)
  call read_bp_file_real(direction_0_file, kernel_names, direction_0)

  if (myrank == 0) write(*, *) "|<----- Calculate Jacobian ----->|"
  call calculate_jacobian_matrix(solver_file, jacobian)

  if (myrank == 0) write(*, *) "|<----- Compute Search Direction ----->|"
  call compute_search_direction(gradient_0, gradient_1, direction_0, jacobian, direction_1)

  if (myrank == 0) write(*, *) "|<----- Start Writing ----->|"
  call write_bp_file(direction_1, kernel_names, "KERNEL_GROUPS", direction_1_file)

  if (myrank == 0) write(*, *) "|<----- Done Writing ----->|"

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program main
