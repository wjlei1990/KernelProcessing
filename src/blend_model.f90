module model_blend_subs

  use mpi
  use AdiosIO
  use global_var, only : myrank, nprocs, NGLLX, NGLLY, NGLLZ, NSPEC, NGLOB, CUSTOM_REAL
  use global_var, only : DEGREES_TO_RADIANS
  use global_var, only : init_mpi, exit_mpi, max_all_all_cr, min_all_all_cr

  integer, parameter :: NPARS = 6
  character(len=500), dimension(NPARS), parameter :: parameters = &
    (/character(len=500) :: "reg1/vpv", "reg1/vph", "reg1/vsv", "reg1/vsh", "reg1/eta", "reg1/rho"/)
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NPARS) :: models0, models1, models

  ! model array
  integer, dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: x_glob, y_glob, z_glob

  contains

  subroutine get_distance(x, y, z, x0, y0, z0, dist, distv, disth)
    real(kind=CUSTOM_REAL), intent(in) :: x, y, z, x0, y0, z0
    real(kind=CUSTOM_REAL), intent(inout) :: dist, distv, disth

    r0 = sqrt( x0 * x0 + y0 * y0 + z0 * z0)
    r1 = sqrt( x * x + y * y  + z * z )

    ratio = ( x * x0 + y * y0 + z * z0 ) / (r0 * r1)
    if( ratio > 1.0_CUSTOM_REAL ) ratio = 1.0_CUSTOM_REAL
    if( ratio < -1.0_CUSTOM_REAL ) ratio = -1.0_CUSTOM_REAL

    theta = acos( ratio )
    disth = r1 * sin(theta)
    distv = r0 - r1 * cos(theta)

    dist = sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)

  end subroutine get_distance

  subroutine convert_to_xyz(r, lat, lon, x, y, z)
    real(kind=CUSTOM_REAL), intent(in) :: r, lat, lon
    real(kind=CUSTOM_REAL), intent(out) :: x, y, z

    double precision :: phi, theta

    phi = lon * DEGREES_TO_RADIANS
    theta = (90.0 - lat) * DEGREES_TO_RADIANS

    x = REAL(r * sin(theta) * cos(phi), CUSTOM_REAL)
    y = REAL(r * sin(theta) * sin(phi), CUSTOM_REAL)
    z = REAL(r * cos(theta), CUSTOM_REAL)

  end subroutine convert_to_xyz

  subroutine modify_coef(r_top, r_bot, radius, lat0, lon0, coefs)
    real(kind=CUSTOM_REAL), intent(in) :: r_top, r_bot, radius, lat0, lon0
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(inout) :: coefs

    integer :: i, j, k, ispec
    real(kind=CUSTOM_REAL) :: dist, distv, disth

    call convert_to_xyz(r_top, lat0, lon0, x0_top, y0_top, z0_top)
    call convert_to_xyz(r_bot, lat0, lon0, x0_bot, y0_bot, z0_bot)

    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            gauss = 0.0
            iglob = ibool(i, j, k, ispec)
            x = x_glob(iglob)
            y = y_glob(iglob)
            z = z_glob(iglob)
            r = sqrt( x * x + y * y  + z * z )
            if (r >= r_top) then
              call get_distance(x, y, z, x0_top, y0_top, z0_top, dist, distv, disth)
            else if (r <= r_bot) then
              call get_distance(x, y, z, x0_bot, y0_bot, z0_bot, dist, distv, disth)
            else
              call convert_to_xyz(r, lat0, lon0, x0_tmp, y0_tmp, z0_tmp)
              call get_distance(x, y, z, x0_tmp, y0_tmp, z0_tmp, dist, distv, disth)
            endif

            if(dist > 6 * radius) cycle
            gauss = exp( - (disth / radius) ** 2 - (distv / radius) ** 2)
            coefs(i, j, k, ispec) = coefs(i, j, k, ispec) + gauss
            ! max value of coefs is 1.0
            if (coefs(i, j, k, ispec) > 1.0) coefs(i, j, k, ispec) = 1.0_CUSTOM_REAL
          enddo
        enddo
      enddo
    enddo

  end subroutine modify_coef

  subroutine write_coef_bp_file(coefs, output_file)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: coefs
    character(len=*), intent(in) :: output_file

    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 2) :: coefs_bp

    character(len=500), dimension(2), parameter :: coef_names = &
      (/character(len=500) :: "reg1/coef0", "reg1/coef1"/)

    coefs_bp(:, :, :, :, 1) = coefs
    coefs_bp(:, :, :, :, 2) = 1.0 - coefs

    vmax = 0.0
    call max_all_all_cr(maxval(coefs), vmax)
    vmin = 0.0
    call min_all_all_cr(minval(coefs), vmin)
    if(myrank == 0) then
      write(*, '(A, F10.2, A, F10.2)') "Coef min and max value: ", vmin, ",", vmax
    endif

    if (myrank == 0) print*, "coef file: ", trim(output_file)
    call write_bp_file(coefs_bp, coef_names, "COEF_GROUP", output_file)

  end subroutine write_coef_bp_file


  subroutine blend_model(model0_file, model1_file, coef0, outputfile)
    character(len=*), intent(in) :: model0_file, model1_file
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: coef0
    character(len=*), intent(in) :: outputfile

    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: coef1
    ! read in the model0 and model1
    call read_bp_file_real(model0_file, parameters, models0)
    call read_bp_file_real(model1_file, parameters, models1)

    coef1 = 1.0 - coef0
    do i = 1, NPARS
      models(:,:,:,:,i) = coef0 * models0(:,:,:,:,i) + coef1 * models1(:,:,:,:,i)
    enddo
    if(myrank == 0) write(*, '(A, A)') "output blend model: ", outputfile
    call write_bp_file(models, parameters, "KERNELS_GROUP", outputfile)
  end subroutine blend_model


end module model_blend_subs


program main
  use mpi
  use adios_read_mod
  use model_blend_subs

  implicit none

  ! coef for model that to be blended in
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: coefs
  character(len=500) :: outputdir, outputfile
  character(len=500) :: input_solver_file, model0_file, model1_file

  real(kind=CUSTOM_REAL) :: r_top, r_bot, radius, lat0, lon0
  integer :: ier

  outputdir = "/lustre/atlas2/geo111/proj-shared/Wenjie/model_blend"
  model0_file = "/lustre/atlas2/geo111/proj-shared/Wenjie/model_blend/M15/model_gll.bp"
  model1_file = "/lustre/atlas2/geo111/proj-shared/Wenjie/model_blend/M25/model_gll.bp"
  input_solver_file = &
    "/lustre/atlas2/geo111/proj-shared/Wenjie/model_blend/M25/DATABASES_MPI/solver_data.bp"

  call init_mpi()

  if(myrank == 0) write(*, '(A)') "|<------ Read Model Files ------>|"
  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, "verbose=1", ier)
  if(myrank == 0) write(*, '(A, A)')  "Read from solver file: ", trim(input_solver_file)
  call read_bp_file_int(input_solver_file, "reg1/ibool", ibool)
  call read_bp_file_real(input_solver_file, "reg1/x_global", x_glob)
  call read_bp_file_real(input_solver_file, "reg1/y_global", y_glob)
  call read_bp_file_real(input_solver_file, "reg1/z_global", z_glob)

  ! init coef
  coefs = 0.0
  if(myrank == 0) write(*, '(A)') "|<------ Calculate Coef ------>|"

  r_top = 0.84
  r_bot = 0.77
  radius = 0.05
  lat0 = -14.64
  lon0 = 170.70
  call modify_coef(r_top, r_bot, radius, lat0, lon0, coefs)

  r_top = 0.87
  r_bot = 0.87
  radius = 0.04
  lat0 = 53.00
  lon0 = 152.00
  call modify_coef(r_top, r_bot, radius, lat0, lon0, coefs)

  if(myrank == 0) write(*, '(A)') "|<------ Write Coef ------>|"
  outputfile = trim(outputdir) // "/blend_coef.bp"
  call write_coef_bp_file(coefs, outputfile)

  outputfile = trim(outputdir) // "/model_gll-blend.bp"
  call blend_model(model0_file, model1_file, coefs, outputfile)

  if(myrank == 0) write(*, '(A)') "|<------ Done ------>|"
  call adios_finalize(myrank, ier)
  call MPI_Finalize(ier)

end program main
