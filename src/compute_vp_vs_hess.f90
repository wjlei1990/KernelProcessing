! Ebru1: The choice of THRESHOLD value is somewhat subjective. It is not trivial to set it like the 20% of max value
! which may be OK for smaller scale studies but global scale needs a few trial&error to adjust this parameter for
! every iteration. Needs some more investigation..

! Ebru2: I find the preconditioner behave better after changing the order of smoothing and preconditioning in
! post-processing upon the suggestion by Ryan & Yanhua.
! However, I am still not convinced by Ryan's latest suggestion that preconditioner should be smoothed more than the
! gradients of other parameters that the preconditioner to be applied. I currently smooth the preconditioner and
! the other gradients in the same way.

module compute_vp_vs_hess_sub
  use mpi
  use global_var, only : max_all_all_cr, min_all_all_cr, CUSTOM_REAL, exit_mpi, &
                     myrank
  implicit none

  contains

  subroutine get_sys_args(input_file, model_file, output_file)
    character(len=*), intent(inout) :: input_file, model_file, output_file

    call getarg(1, input_file)
    call getarg(2, model_file)
    call getarg(3, output_file)

    if(input_file == '' .or. model_file == '' .or. output_file == '') then
      call exit_mpi("Usage: xcompute_vp_vs_hess input_kernel solver_file output_kernel")
    endif

    if(myrank == 0) then
      write(*, *) "Input kernel: ", trim(input_file)
      write(*, *) "Model file: ", trim(model_file)
      write(*, *) "Output kernel: ", trim(output_file)
    endif

  end subroutine get_sys_args

  subroutine normalize_hess(hess, hess_names)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: hess
    character(len=*), dimension(:), intent(in) :: hess_names

    real(kind=CUSTOM_REAL) :: maxh_all, minh_all
    integer :: i

    call max_all_all_cr(maxval(hess), maxh_all)
    call min_all_all_cr(minval(hess), minh_all)

    if (myrank==0) then
      write(*, '(A, ES12.6, ES12.6, F12.2)') "Max and Min of hess: ", maxh_all, minh_all, maxh_all/minh_all
      write(*, *) 'Normalize factor(max hess) for all processors ', maxh_all
    endif

    ! normalized hess
    hess = hess / maxh_all

    call max_all_all_cr(maxval(hess), maxh_all)
    call min_all_all_cr(minval(hess), minh_all)

    if (myrank==0) then
      write(*, '(A, ES12.6)') 'min and max hess after norm: ', minh_all, maxh_all
      write(*, *) "Number of hess kernels: ", size(hess_names)
    endif

    do i=1, size(hess_names)
      call max_all_all_cr(maxval(hess(:, :, :, :, i)), maxh_all)
      call min_all_all_cr(minval(hess(:, :, :, :, i)), minh_all)
      if (myrank == 0) then
        write(*, '(A, A, A, ES12.6, 5X, ES12.6, 5X, ES12.6)') '[', trim(hess_names(i)), &
          '] min, max: ', minh_all, maxh_all, maxh_all / minh_all
      endif
    enddo

  end subroutine normalize_hess


end module compute_vp_vs_hess_sub

program compute_vp_vs_hess
  use mpi
  use adios_read_mod
  use AdiosIO
  use global_var, only : NGLLX, NGLLY, NGLLZ, NSPEC, myrank, CUSTOM_REAL
  use global_var, only : init_mpi
  use compute_vp_vs_hess_sub

  implicit none

  integer, parameter :: NKERNELS = 2
  character(len=500), parameter :: kernel_names(NKERNELS) = &
    (/character(len=500) :: "hess_kappa_kl_crust_mantle", "hess_mu_kl_crust_mantle"/)

  integer, parameter :: NKERNELS_out = 3
  character(len=500), parameter :: kernel_out_names(NKERNELS_out) = &
    (/character(len=500) :: "hess_vp_kl_crust_mantle", "hess_vs_kl_crust_mantle", &
                            "hess_eta_kl_crust_mantle"/)

  integer, parameter :: NMODELS = 3
  character(len=500), parameter :: model_names(NMODELS) = &
    (/character(len=500) :: "reg1/kappavstore", "reg1/muvstore", "reg1/muhstore"/)


  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS) :: hess_in
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: hess_kappa, hess_mu

  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NMODELS) :: models
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: kappa, mu

  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS_out) :: hess_out
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: hess_vp, hess_vs, hess_eta

  character(len=500) :: input_file, model_file, output_file
  integer:: ier

  call init_mpi()

  if (myrank == 0) then
    write(*, *) "Compute Vp and Vs Hessian kernels"
  endif

  call get_sys_args(input_file, model_file, output_file)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  call read_bp_file_real(input_file, kernel_names, hess_in)
  call read_bp_file_real(model_file, model_names, models)

  kappa = models(:, :, :, :, 1)
  mu = sqrt((2.0*models(:, :, :, :, 2)**2 + models(:, :, :, :, 3)**2) / 3.0)

  ! WARNING: tmp modification. only for first iteration since error in SPECFEM
  ! REMOVE LATER !!!
  !hess_in(:, :, :, :, 1) = 3 * hess_in(:, :, :, :, 1)

  hess_kappa = hess_in(:, :, :, :, 1)
  hess_mu = hess_in(:, :, :, :, 2)

  ! compute for hess_vp, hess_vs, and hess_eta
  hess_vp = 4.0 * (1.0 + 4.0/3.0*(mu/kappa))**2 * hess_kappa
  hess_vs = 4.0 * hess_mu + 64.0 / 9.0 * (mu/kappa)**2 * hess_kappa
  hess_eta = hess_kappa + 4.0 / 9.0 * hess_mu

  hess_out(:, :, :, :, 1) = hess_vp
  hess_out(:, :, :, :, 2) = hess_vs
  hess_out(:, :, :, :, 3) = hess_eta

  call normalize_hess(hess_out, kernel_out_names)

  call write_bp_file(hess_out, kernel_out_names, "KERNEL_GOURPS", output_file)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program compute_vp_vs_hess
