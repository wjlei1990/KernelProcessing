! Ebru1: The choice of THRESHOLD value is somewhat subjective. It is not trivial to set it like the 20% of max value
! which may be OK for smaller scale studies but global scale needs a few trial&error to adjust this parameter for
! every iteration. Needs some more investigation..

! Ebru2: I find the preconditioner behave better after changing the order of smoothing and preconditioning in
! post-processing upon the suggestion by Ryan & Yanhua.
! However, I am still not convinced by Ryan's latest suggestion that preconditioner should be smoothed more than the
! gradients of other parameters that the preconditioner to be applied. I currently smooth the preconditioner and
! the other gradients in the same way.


module precond_vp_vs_kernels_sub
  use mpi
  use global_var, only : max_all_all_cr, min_all_all_cr, CUSTOM_REAL, exit_mpi, &
                     myrank
  implicit none

  contains

  subroutine get_sys_args(input_file, hess_file, output_file, threshold_hess)
    character(len=*), intent(inout) :: input_file, output_file, hess_file
    real(kind=CUSTOM_REAL), intent(inout) :: threshold_hess

    character(len=20) :: threshold_str

    call getarg(1, input_file)
    call getarg(2, hess_file)
    call getarg(3, output_file)
    call getarg(4, threshold_str)

    if(input_file == '' .or. hess_file == '' .or. output_file == '' .or. &
       threshold_str == '' ) then
      call exit_mpi("Usage: xprecond_vp_vs_kernels input_kernel hess_kernel "// &
        "output_kernel threshold")
    endif

    read(threshold_str, *) threshold_hess

    if (myrank == 0) then
      write(*, *) "Input kernel: ", trim(input_file)
      write(*, *) "Hess kernel: ", trim(hess_file)
      write(*, *) "Output kernel: ", trim(output_file)
      write(*, *) "Threshold hessian: ", threshold_hess
    endif

  end subroutine get_sys_args

  subroutine prepare_hessian(hess, hess_names, threshold, invHess)
    real(CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: hess, invHess
    character(len=*), dimension(:) :: hess_names
    real(CUSTOM_REAL), intent(in) :: threshold

    real(kind=CUSTOM_REAL):: maxh_all, minh_all
    integer :: i

    call max_all_all_cr(maxval(hess), maxh_all)
    call min_all_all_cr(minval(hess), minh_all)

    if ( maxh_all < 1.e-18 ) then
      call exit_mpi("hess max value < 1.e-18")
    end if

    if (myrank==0) then
      write(*, *) "Max and Min of hess: ", maxh_all, minh_all
      write(*, *) 'Normalize factor(max hess) for all processors ', maxh_all
    endif

    ! normalized hess
    hess = hess / maxh_all

    call max_all_all_cr(maxval(hess), maxh_all)
    call min_all_all_cr(minval(hess), minh_all)

    if (myrank==0) then
      write(*, '(A, E12.6)') 'min and max hess after norm', minh_all, maxh_all
      write(*, '(A, E12.6)') "Hessian Threshold: ", threshold
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

    where(hess > threshold )
      invHess = 1.0_CUSTOM_REAL / hess
    elsewhere
      invHess = 1.0_CUSTOM_REAL / threshold
    endwhere

  end subroutine prepare_hessian

end module precond_vp_vs_kernels_sub

program precond_vp_vs_kernels
  use mpi
  use adios_read_mod
  use AdiosIO
  use global_var, only : NGLLX, NGLLY, NGLLZ, NSPEC, myrank, CUSTOM_REAL
  use global_var, only : init_mpi
  use precond_vp_vs_kernels_sub

  implicit none

  integer, parameter :: NKERNELS = 4    !bulk_betah, bulk_betav, bulk_c, eta
  character(len=500), parameter :: kernel_names(NKERNELS) = &
    (/character(len=500) :: "bulk_c_kl_crust_mantle", "bulk_betav_kl_crust_mantle", &
                            "bulk_betah_kl_crust_mantle", "eta_kl_crust_mantle"/)
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS) :: kernels = 0.0
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS) :: kernels_precond = 0.0

  integer, parameter :: NHESS = 3
  character(len=500), parameter :: hess_names(NHESS) = &
    (/character(len=500) :: "hess_vp_kl_crust_mantle", "hess_vs_kl_crust_mantle", &
                            "hess_eta_kl_crust_mantle"/)
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NHESS) :: hess = 0.0, invHess = 0.0

  character(len=500) :: input_file, output_file, hess_file
  real(kind=CUSTOM_REAL) :: threshold_hess
  integer:: ier

  call init_mpi()

  if (myrank == 0) print*, "Precondition using Vp and Vs Hessian Kernels"

  call get_sys_args(input_file, hess_file, output_file, threshold_hess)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  call read_bp_file_real(input_file, kernel_names, kernels)
  call read_bp_file_real(hess_file, hess_names, hess)

  call prepare_hessian(hess, hess_names, threshold_hess, invHess)

  ! bulk_c
  kernels_precond(:, :, :, :, 1) = kernels(:, :, :, :, 1) * invHess(:, :, :, :, 1)
  ! betav and betah
  kernels_precond(:, :, :, :, 2) = kernels(:, :, :, :, 2) * invHess(:, :, :, :, 2)
  kernels_precond(:, :, :, :, 3) = kernels(:, :, :, :, 3) * invHess(:, :, :, :, 2)
  ! eta
  kernels_precond(:, :, :, :, 4) = kernels(:, :, :, :, 4) * invHess(:, :, :, :, 3)

  call write_bp_file(kernels_precond, kernel_names, "KERNEL_GOURPS", output_file)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program precond_vp_vs_kernels
