! Ebru1: The choice of THRESHOLD value is somewhat subjective. It is not trivial to set it like the 20% of max value
! which may be OK for smaller scale studies but global scale needs a few trial&error to adjust this parameter for
! every iteration. Needs some more investigation..

! Ebru2: I find the preconditioner behave better after changing the order of smoothing and preconditioning in
! post-processing upon the suggestion by Ryan & Yanhua.
! However, I am still not convinced by Ryan's latest suggestion that preconditioner should be smoothed more than the
! gradients of other parameters that the preconditioner to be applied. I currently smooth the preconditioner and
! the other gradients in the same way.


module preconditioner_subs
  use mpi
  use global_var, only : max_all_all_cr, min_all_all_cr, CUSTOM_REAL, exit_mpi, &
                     myrank
  implicit none

  contains

  subroutine get_sys_args(hess_file, output_file, threshold_hess)
    character(len=*), intent(inout) :: output_file, hess_file
    real(kind=CUSTOM_REAL), intent(inout) :: threshold_hess

    character(len=20) :: threshold_str

    call getarg(1, hess_file)
    call getarg(2, output_file)
    call getarg(3, threshold_str)

    if(hess_file == '' .or. output_file == '' .or. threshold_str == '' ) then
      call exit_mpi("Usage: xprepare_vp_vs_precond hess_file "// &
        "output_kernel threshold")
    endif

    read(threshold_str, *) threshold_hess

    if (myrank == 0) then
      write(*, *) "Hess kernel: ", trim(hess_file)
      write(*, *) "Preconditioner(output): ", trim(output_file)
      write(*, *) "Dampening threshold: ", threshold_hess
    endif

  end subroutine get_sys_args

  ! Take the inverse of Hessian with dampening (as preconditioner kernels)
  ! This method will apply different lambda for invidual model parameters
  ! but using the same threshold ratio regarding to its own max value
  subroutine prepare_invHess(hess, hess_names, threshold, invHess)

    real(CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: hess, invHess
    character(len=*), dimension(:) :: hess_names
    real(CUSTOM_REAL), intent(in) :: threshold

    real(kind=CUSTOM_REAL):: maxh_all, minh_all, dampening
    integer :: ik
    integer :: nhess

    nhess = size(hess_names)
    do ik = 1, nhess
      if (myrank == 0) then
        write(*, '(/,A,A)') "====> Prepare invHess for parameter: ", trim(hess_names(ik))
      endif

      ! stats before
      call max_all_all_cr(maxval(hess(:, :, :, :, ik)), maxh_all)
      call min_all_all_cr(minval(hess(:, :, :, :, ik)), minh_all)

      if ( maxh_all < 1.e-18 ) then
        call exit_mpi("hess max value < 1.e-18")
      end if

      ! damping value
      dampening = maxh_all * threshold
      if (myrank == 0) then
        write(*, '(A, ES12.2, 5X, ES12.2, 5X, A, ES12.2)') &
          "Max and Min of hess: ", maxh_all, minh_all, &
          " || Condition Number:", maxh_all / minh_all
        write(*, '(A, E12.2, 5X, E12.2)') "Dampen threshold and value: ", threshold, dampening
      endif

      invHess(:, :, :, :, ik) = 1.0 / (hess(:, :, :, :, ik) + dampening)

      ! stats afterwards
      call max_all_all_cr(maxval(invHess(:, :, :, :, ik)), maxh_all)
      call min_all_all_cr(minval(invHess(:, :, :, :, ik)), minh_all)

      if (myrank == 0) then
        write(*, *) "Apply dampening and inverse..."
        write(*, '(A, ES12.2, 5X, ES12.2, 5X, A, ES12.2)') &
          "Max and Min of hess: ", maxh_all, minh_all, &
          " || Condition Number:", maxh_all / minh_all
      endif

    enddo
  end subroutine prepare_invHess

end module preconditioner_subs

program prepare_vp_vs_preconditioner
  use mpi
  use adios_read_mod
  use AdiosIO
  use global_var, only : NGLLX, NGLLY, NGLLZ, NSPEC, myrank, CUSTOM_REAL
  use global_var, only : init_mpi
  use preconditioner_subs

  implicit none

  integer, parameter :: NHESS = 3
  character(len=500), parameter :: hess_names(NHESS) = &
    (/character(len=500) :: "hess_vp_kl_crust_mantle", &
                            "hess_vs_kl_crust_mantle", &
                            "hess_eta_kl_crust_mantle"/)
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NHESS) :: hess = 0.0, invHess = 0.0

  integer, parameter :: NPRECOND = 4
  character(len=500), parameter :: precond_names(NPRECOND) = &
    (/character(len=500) :: "precond_bulk_c_kl_crust_mantle", &
                            "precond_bulk_betav_kl_crust_mantle", &
                            "precond_bulk_betah_kl_crust_mantle", &
                            "precond_eta_kl_crust_mantle"/)
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NPRECOND) :: precond = 0.0

  character(len=500) :: hess_file, output_file
  real(kind=CUSTOM_REAL) :: threshold_hess
  integer:: ier

  call init_mpi()

  if (myrank == 0) print*, "Prepare Preconditions for c, Vp, Vs, eta"

  call get_sys_args(hess_file, output_file, threshold_hess)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)
  call read_bp_file_real(hess_file, hess_names, hess)

  call prepare_invHess(hess, hess_names, threshold_hess, invHess)

  ! bulk_c precond
  precond(:, :, :, :, 1) = invHess(:, :, :, :, 1)
  ! betav
  precond(:, :, :, :, 2) = invHess(:, :, :, :, 2)
  ! betah
  precond(:, :, :, :, 3) = invHess(:, :, :, :, 2)
  ! eta
  precond(:, :, :, :, 4) = invHess(:, :, :, :, 3)

  call write_bp_file(precond, precond_names, "KERNEL_GOURPS", output_file)
  if(myrank == 0) then
    write(*, '(/, A, A)') "output preconditioner: ", trim(output_file)
  endif

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program prepare_vp_vs_preconditioner
