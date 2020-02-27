! ----------------------------------
! gaussian perturbation utils
! ----------------------------------

module gpp_utils
  use mpi
  use AdiosIO
  use adios_read_mod
  use global_var, only : myrank, nprocs, NGLLX, NGLLY, NGLLZ, NSPEC, NGLOB, CUSTOM_REAL
  use global_var, only : R_EARTH_KM, DEGREES_TO_RADIANS
  use global_var, only : init_mpi, exit_mpi, max_all_all_cr, min_all_all_cr

  implicit None

  integer, parameter :: NKERNELS = 4
  character(len=500), dimension(NKERNELS), parameter :: kernel_names = &
    (/character(len=500) :: "bulk_c_kl_crust_mantle", &
                            "bulk_betav_kl_crust_mantle", &
                            "bulk_betah_kl_crust_mantle", &
                            "eta_kl_crust_mantle"/)
  integer, parameter :: bulk_c_kl_idx = 1, betav_kl_idx = 2, betah_kl_idx = 3, &
                        eta_kl_idx = 4

  ! Kernel array
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: kernels
  ! model array
  integer, dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: x_glob, y_glob, z_glob

  contains

  subroutine get_sys_args(input_solver_file, output_file)
    character(len=*), intent(inout) :: input_solver_file, output_file

    call getarg(1, input_solver_file)
    call getarg(2, output_file)

    if (trim(input_solver_file) == '' .or. trim(output_file) == '') then
      call exit_mpi('Usage: mpiexec -n 384 ./bin/gaussian_perturb input_solver_file output_file')
    endif

    if (myrank == 0) then
      print*, "input solver file: ", trim(input_solver_file)
      print*, "output file: ", trim(output_file)
    endif

  end subroutine get_sys_args

  subroutine convert_to_xyz(dep, lat, lon, x, y, z)
    real(kind=CUSTOM_REAL), intent(in) :: dep, lat, lon
    real(kind=CUSTOM_REAL), intent(out) :: x, y, z

    double precision :: r1, phi, theta

    phi = lon * DEGREES_TO_RADIANS
    theta = (90.0 - lat) * DEGREES_TO_RADIANS
    r1 = (R_EARTH_KM - dep ) / R_EARTH_KM

    x = REAL(r1 * sin(theta) * cos(phi), CUSTOM_REAL)
    y = REAL(r1 * sin(theta) * sin(phi), CUSTOM_REAL)
    z = REAL(r1 * cos(theta), CUSTOM_REAL)

  end subroutine convert_to_xyz

  subroutine add_gaussian_perturb_sphere(r, lat, lon, sigma0, perturb_idx)
    ! version by Ebru Bozdag
    ! Pure spherical Gaussian perturbation
    real(kind=CUSTOM_REAL), intent(in) :: r, lat, lon, sigma0
    integer, intent(in) :: perturb_idx

    real(kind=CUSTOM_REAL) :: xcent, ycent, zcent, sigma
    integer :: ispec, i, j, k, iglob
    real(kind=CUSTOM_REAL) :: gauss, dist
    real(kind=CUSTOM_REAL) :: vmax, vmin

    if (myrank == 0) then
      write(*, '(A, F8.2, A, F8.2, A, F8.2, A, F8.2, A)') &
        "==>  Gaussian Sphere(Depth=", r, "km, lat=", lat, ", lon=", lon, ", sigma=", sigma0, ")"
    endif

    ! convert sigma from km to xyz coordinate
    sigma = REAL(sigma0 / R_EARTH_KM, CUSTOM_REAL)
    call convert_to_xyz(r, lat, lon, xcent, ycent, zcent)
    if (myrank == 0) then
      write(*, '(A, ES18.8, A, ES18.8, A, ES18.8)') &
        "|    xyz: ", xcent, ",", ycent, ",", zcent, 'sigma:', sigma
      write(*, '(A, ES18.8, A, ES18.8, A)') &
        '|    Param(R_EARTH=', R_EARTH_KM, "km, DEGREE_TO_RADIANS=", DEGREES_TO_RADIANS, ")"
    endif

    do ispec = 1,NSPEC
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob=ibool(i,j,k,ispec)

            dist = sqrt( (x_glob(iglob) - xcent) ** 2 + (y_glob(iglob) - ycent) ** 2 + &
                         (z_glob(iglob) - zcent) ** 2 )
            ! postive search direction --> fast anomaly
            !gauss = exp( - (dist/sigma) **2 )
            ! negative search direction --> slow anomaly
            gauss = -exp( - (dist/sigma) **2 )

            kernels(i,j,k,ispec,perturb_idx) = kernels(i,j,k,ispec,perturb_idx) + gauss
          end do
        end do
      end do
    end do

    vmax = 0.0
    call max_all_all_cr(maxval(kernels(:,:,:,:,perturb_idx)), vmax)
    vmin = 0.0
    call min_all_all_cr(minval(kernels(:,:,:,:,perturb_idx)), vmin)
    if(myrank == 0) then
      write(*, '(A, F10.2, A, F10.2)') "Kernel min and max value: ", vmin, ",", vmax
    endif
  end subroutine add_gaussian_perturb_sphere

  subroutine add_gaussian_perturb_elliptic(r, lat, lon, sigmah0, sigmav0, perturb_idx)
    ! Version By Wenjie
    ! Create elliptical
    real(kind=CUSTOM_REAL), intent(in) :: r, lat, lon, sigmah0, sigmav0
    integer, intent(in) :: perturb_idx

    real(kind=CUSTOM_REAL) :: xcent, ycent, zcent, sigmah, sigmav
    integer :: ispec, i, j, k, iglob
    real(kind=CUSTOM_REAL) :: gauss, distv, disth, r0, r1, ratio, theta
    real(kind=CUSTOM_REAL) :: vmax, vmin

    if (myrank == 0) then
      write(*, '(A)') "|--> Gaussian HV"
      write(*, '(A, F8.2, A, F8.2, A, F8.2, A, F8.2, A, F8.2, A)') &
        "|    Depth=", r, "km, lat=", lat, ", lon=", lon, &
        ", sigmav0=", sigmav0, "km, sigmah0=", sigmah0, "km"
    endif

    ! convert sigma from km to xyz coordinate
    sigmah = REAL(sigmah0 / R_EARTH_KM, CUSTOM_REAL)
    sigmav = REAL(sigmav0 / R_EARTH_KM, CUSTOM_REAL)

    call convert_to_xyz(r, lat, lon, xcent, ycent, zcent)
    if (myrank == 0) then
      write(*, '(A, ES18.3, A, ES18.3, A, ES18.3)') &
        "|    xyz: ", xcent, ",", ycent, ",", zcent
      write(*, '(A, ES18.3, A, ES18.3)') &
        '|    sigmah:', sigmah, ", sigmav:", sigmav
      write(*, '(A, ES18.3, A, ES18.3)') &
        '|    R_EARTH=', R_EARTH_KM, "km, DEGREE_TO_RADIANS=", DEGREES_TO_RADIANS
    endif

    r0=sqrt( xcent * xcent + ycent * ycent + zcent * zcent)
    do ispec = 1,NSPEC
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob=ibool(i,j,k,ispec)

            r1=sqrt( x_glob(iglob) * x_glob(iglob) + y_glob(iglob) * y_glob(iglob) &
                    + z_glob(iglob) * z_glob(iglob) )

            ratio = ( x_glob(iglob) * xcent + y_glob(iglob) * ycent + z_glob(iglob) * zcent ) &
                    / (r0*r1)
            ! checks boundaries of ratio (due to numerical inaccuracies)
            if( ratio > 1.0_CUSTOM_REAL ) ratio = 1.0_CUSTOM_REAL
            if( ratio < -1.0_CUSTOM_REAL ) ratio = -1.0_CUSTOM_REAL

            theta = acos( ratio )

            disth = r1 * sin(theta)
            distv = r0 - r1 * cos(theta)

            ! postive search direction --> fast anomaly
            !gauss = exp(-(disth / sigmah) ** 2 - (distv / sigmav) ** 2)
            ! negative search direction --> slow anomaly
            gauss = -exp(- (disth / sigmah) ** 2 - (distv / sigmav) ** 2)

            kernels(i,j,k,ispec,perturb_idx) = kernels(i,j,k,ispec,perturb_idx) + gauss
          end do
        end do
      end do
    end do

    vmax = 0.0
    call max_all_all_cr(maxval(kernels(:,:,:,:,perturb_idx)), vmax)
    vmin = 0.0
    call min_all_all_cr(minval(kernels(:,:,:,:,perturb_idx)), vmin)
    if(myrank == 0) then
      write(*, '(A, ES18.8, A, ES18.8)') "Kernel min and max value: ", vmin, ",", vmax
    endif
  end subroutine add_gaussian_perturb_elliptic

  subroutine add_gaussian_perturb_hv(dep, lat, lon, perturb_sign, sigmah0, sigmav0, perturb_idx)
    ! Version by Hejun
    ! assign different horizontal and vertical length
    ! This one should be used with careful check because the shape of
    ! the preturbations will be doom-liked sphere.(not sysmetircal with
    ! regarding to the top and bottom)
    real(kind=CUSTOM_REAL), intent(in) :: dep, lat, lon, perturb_sign, sigmah0, sigmav0
    integer, intent(in) :: perturb_idx

    real(kind=CUSTOM_REAL) :: xcent, ycent, zcent, sigmah, sigmav
    integer :: ispec, i, j, k, iglob
    real(kind=CUSTOM_REAL) :: gauss, dist, distv, disth, r0, r1, ratio, theta
    real(kind=CUSTOM_REAL) :: vmax, vmin

    if (myrank == 0) then
      write(*, '(A)') "|--> Gaussian HV"
      write(*, '(A, F8.2, A, F8.2, A, F8.2, A, F8.2, A, F8.2, A)') &
        "|    Depth=", dep, "km, lat=", lat, ", lon=", lon, &
        ", sigmav0=", sigmav0, "km, sigmah0=", sigmah0, "km"
    endif

    ! convert sigma from km to xyz coordinate
    sigmah = REAL(sigmah0 / R_EARTH_KM, CUSTOM_REAL)
    sigmav = REAL(sigmav0 / R_EARTH_KM, CUSTOM_REAL)

    call convert_to_xyz(dep, lat, lon, xcent, ycent, zcent)
    if (myrank == 0) then
      write(*, '(A, ES18.3, A, ES18.3, A, ES18.3)') &
        "|    xyz: ", xcent, ",", ycent, ",", zcent
      write(*, '(A, ES18.3, A, ES18.3)') &
        '|    sigmah:', sigmah, ", sigmav:", sigmav
      write(*, '(A, ES18.3, A, ES18.3)') &
        '|    R_EARTH=', R_EARTH_KM, "km, DEGREE_TO_RADIANS=", DEGREES_TO_RADIANS
    endif

    r0=sqrt( xcent * xcent + ycent * ycent + zcent * zcent)
    do ispec = 1,NSPEC
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob=ibool(i,j,k,ispec)
            dist = sqrt( (x_glob(iglob) - xcent) ** 2 + (y_glob(iglob) - ycent) ** 2 + &
                         (z_glob(iglob) - zcent) ** 2 )

            r1=sqrt( x_glob(iglob) * x_glob(iglob) + y_glob(iglob) * y_glob(iglob) &
                    + z_glob(iglob) * z_glob(iglob) )
            distv = r1 - r0

            ratio = ( x_glob(iglob) * xcent + y_glob(iglob) * ycent + z_glob(iglob) * zcent ) &
                    / (r0*r1)
            ! checks boundaries of ratio (due to numerical inaccuracies)
            if( ratio > 1.0_CUSTOM_REAL ) ratio = 1.0_CUSTOM_REAL
            if( ratio < -1.0_CUSTOM_REAL ) ratio = -1.0_CUSTOM_REAL

            theta = acos( ratio )
            disth = r1 * theta

            gauss = perturb_sign * exp(-(disth / sigmah) ** 2 - (distv / sigmav) ** 2)

            kernels(i,j,k,ispec,perturb_idx) = kernels(i,j,k,ispec,perturb_idx) + gauss
          end do
        end do
      end do
    end do

    vmax = 0.0
    call max_all_all_cr(maxval(kernels(:,:,:,:,perturb_idx)), vmax)
    vmin = 0.0
    call min_all_all_cr(minval(kernels(:,:,:,:,perturb_idx)), vmin)
    if(myrank == 0) then
      write(*, '(A, ES18.8, A, ES18.8)') "Kernel min and max value: ", vmin, ",", vmax
    endif
  end subroutine add_gaussian_perturb_hv

  subroutine read_model_file(input_solver_file)
    character(len=*), intent(in) :: input_solver_file

    integer :: ier

    call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, "verbose=1", ier)
    call read_bp_file_int(input_solver_file, "reg1/ibool", ibool)
    call read_bp_file_real(input_solver_file, "reg1/x_global", x_glob)
    call read_bp_file_real(input_solver_file, "reg1/y_global", y_glob)
    call read_bp_file_real(input_solver_file, "reg1/z_global", z_glob)
  end subroutine read_model_file

end module gpp_utils
