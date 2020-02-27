module global_var

  use mpi
  implicit none

  integer, parameter :: CUSTOM_REAL = 4
  integer, parameter :: CUSTOM_MPI_TYPE = MPI_REAL

  ! values from the mesh
  double precision, parameter :: GAUSSALPHA = 0.d0, GAUSSBETA = 0.d0
  integer, parameter :: NGLLX = 5, NGLLY = NGLLX, NGLLZ = NGLLX
  integer, parameter :: NSPEC_CRUST_MANTLE = 10496
  integer, parameter :: NGLOB_CRUST_MANTLE = 699193
  integer, parameter :: NSPEC = NSPEC_CRUST_MANTLE
  integer, parameter :: NGLOB = NGLOB_CRUST_MANTLE

  character(len=*), parameter :: ADIOS_TRANSPORT_METHOD = "MPI"
  integer, parameter :: ADIOS_BUFFER_SIZE_IN_MB = 200

  ! Earth
  double precision, parameter :: R_EARTH = 6371000.d0
  double precision, parameter :: R_EARTH_KM = R_EARTH / 1000.d0

  real(kind=CUSTOM_REAL), parameter :: FOUR_THIRDS = 4._CUSTOM_REAL/3._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO_THIRDS = 2._CUSTOM_REAL/3._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE_THIRD = 1._CUSTOM_REAL/3._CUSTOM_REAL

  ! flags for idoubling
  integer, parameter :: IFLAG_80_MOHO = 2
  integer, parameter :: IFLAG_220_80 = 3
  integer, parameter :: IFLAG_670_220 = 4

  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: DEGREES_TO_RADIANS = PI / 180.d0

  ! MPI variable that need to be init in the main program
  integer :: myrank, nprocs

  contains

  subroutine sum_all_all_cr(sendbuf, recvbuf)
    real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
    integer :: ier

    call MPI_AllReduce(sendbuf, recvbuf, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ier)
  end subroutine

  subroutine max_all_all_cr(sendbuf, recvbuf)

    real(kind=CUSTOM_REAL):: sendbuf, recvbuf
    integer :: ier

    call MPI_ALLREDUCE(sendbuf, recvbuf, 1, CUSTOM_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD, ier)

  end subroutine max_all_all_cr

  subroutine min_all_all_cr(sendbuf, recvbuf)

    real(kind=CUSTOM_REAL):: sendbuf, recvbuf
    integer :: ier

    call MPI_ALLREDUCE(sendbuf, recvbuf, 1, CUSTOM_MPI_TYPE, MPI_MIN, MPI_COMM_WORLD, ier)

  end subroutine min_all_all_cr

  subroutine build_jacobian(ibool, xix, xiy, xiz, etax, etay, etaz, &
                            gammax, gammay, gammaz, jacobian)
  integer, dimension(:, :, :, :), intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: xix, xiy, xiz
  real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: etax, etay, etaz
  real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: gammax, gammay, gammaz
  real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(inout) :: jacobian

  integer :: i, j, k, ispec, iglob
  real(kind=CUSTOM_REAL) :: xixl, xiyl, xizl, etaxl, etayl, etazl
  real(kind=CUSTOM_REAL) :: gammaxl, gammayl, gammazl, jacobianl

  do ispec = 1, NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! build jacobian
          ! get derivatives of ux, uy and uz with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)
          ! compute the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                    - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                    + xizl*(etaxl*gammayl-etayl*gammaxl))
          jacobian(i, j, k, ispec) = jacobianl
        enddo
      enddo
    enddo
  enddo

  end subroutine build_jacobian

  subroutine build_gll_weight(wgll_cube)
    double precision, dimension(:, :, :), intent(inout) :: wgll_cube

    ! Gauss-Lobatto-Legendre points of integration and weights
    double precision, dimension(NGLLX) :: xigll, wxgll
    double precision, dimension(NGLLY) :: yigll, wygll
    double precision, dimension(NGLLZ) :: zigll, wzgll

    integer :: i, j, k

    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
    call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
           wgll_cube(i,j,k) = wxgll(i)*wygll(j)*wzgll(k)
        enddo
      enddo
    enddo
  end subroutine build_gll_weight

  subroutine Parallel_ComputeInnerProduct(vect1, vect2, Niv, jacobian, qp)
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(in) :: vect1, vect2
    integer,                intent(in)    :: Niv
    real(kind=CUSTOM_REAL), dimension(:, :, :, :) :: jacobian
    real(kind=CUSTOM_REAL), intent(inout) :: qp

    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, Niv) :: wks_1n, wks_2n

    real(kind=CUSTOM_REAL)  :: qp_tmp_single
    ! try double precision
    real(kind=8)            :: jacobianl, weight, qp_tmp
    integer                 :: ipar, i, j, k, ispec
    real(kind=CUSTOM_REAL)  :: coeff, coeff_n1, coeff_n2
    real(kind=8)            :: coeff_n1_dp, coeff_n2_dp

    double precision, dimension(NGLLX, NGLLY, NGLLZ) :: wgll_cube

    !! try normalization to avoid numerical errors
    !call Parallel_ComputeL2normSquare(vect1 , Niv, coeff_n1)
    !call Parallel_ComputeL2normSquare(vect2 , Niv, coeff_n2)

    coeff=maxval(abs(vect1(:,:,:,:,:)))
    call max_all_all_cr(coeff, coeff_n1)
    if (coeff_n1 == 0._CUSTOM_REAL) coeff_n1=1._CUSTOM_REAL
    wks_1n(:,:,:,:,:) = vect1(:,:,:,:,:) / coeff_n1

    coeff=maxval(abs(vect2(:,:,:,:,:)))
    call max_all_all_cr(coeff, coeff_n2)
    if (coeff_n2 == 0._CUSTOM_REAL) coeff_n2=1._CUSTOM_REAL
    wks_2n(:,:,:,:,:) = vect2(:,:,:,:,:) / coeff_n2

    coeff_n1_dp = coeff_n1
    coeff_n2_dp = coeff_n2

    call build_gll_weight(wgll_cube)

    qp_tmp=0._CUSTOM_REAL
    do ipar=1, Niv
      do ispec = 1, NSPEC
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              weight = wgll_cube(i, j, k)
              jacobianl = jacobian(i, j, k, ispec)
              qp_tmp = qp_tmp + jacobianl * weight * &
                       wks_1n(i,j,k,ispec,ipar) * wks_2n(i,j,k,ispec,ipar)
              !qp = qp + jacobianl * weight * &
              !          vect1(i,j,k,ispec,ipar) * vect2(i,j,k,ispec,ipar)
            enddo
          enddo
        enddo
      enddo
    enddo

    qp_tmp_single = REAL(qp_tmp * coeff_n1_dp * coeff_n2_dp, CUSTOM_REAL)
    qp=0.
    call sum_all_all_cr(qp_tmp_single, qp)

  end subroutine Parallel_ComputeInnerProduct

  subroutine Parallel_ComputeL2normSquare(vect1 , Niv, jacobian, qp)
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(in) :: vect1
    integer,                intent(in)    :: Niv
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: jacobian
    real(kind=CUSTOM_REAL), intent(inout) :: qp

    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, Niv) :: wks_1n
    real(kind=CUSTOM_REAL) :: coeff, coeff_n1
    real(kind=CUSTOM_REAL) :: qp_tmp
    real(kind=8)           :: jacobianl, weight, qp_dp, coeff_n1_dp
    integer                :: ipar, i, j, k, ispec

    ! Gauss-Lobatto-Legendre points of integration and weights
    double precision, dimension(NGLLX, NGLLY, NGLLZ) :: wgll_cube

    coeff=maxval(abs(vect1(:,:,:,:,:)))
    call max_all_all_cr(coeff, coeff_n1)

    if (coeff_n1 == 0._CUSTOM_REAL) coeff_n1=1._CUSTOM_REAL

    wks_1n(:,:,:,:,:) = vect1(:,:,:,:,:) / coeff_n1
    coeff_n1_dp=coeff_n1

    call build_gll_weight(wgll_cube)

    qp_dp=0.d0
    do ipar=1,Niv
      do ispec = 1, NSPEC
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              weight = wgll_cube(i, j, k)
              jacobianl = jacobian(i, j, k, ispec)
              qp_dp = qp_dp + jacobianl * weight * (wks_1n(i,j,k,ispec,ipar)**2)
            enddo
           enddo
        enddo
      enddo
    enddo

    qp_tmp = REAL(qp_dp * coeff_n1_dp * coeff_n1_dp, CUSTOM_REAL)
    qp = 0.d0
    call sum_all_all_cr(qp_tmp, qp)

  end subroutine Parallel_ComputeL2normSquare

  subroutine Parallel_ComputeIntegral(vec, jacobian, res)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: vec
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: jacobian
    real(kind=CUSTOM_REAL), intent(inout) :: res

    double precision, dimension(NGLLX, NGLLY, NGLLZ) :: wgll_cube
    integer :: i, j, k, ispec
    real(kind=8) :: integl, weight, jacobianl
    real(kind=4) :: integl_tmp

    call build_gll_weight(wgll_cube)

    integl = 0.0
    do ispec = 1, NSPEC
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            weight = wgll_cube(i, j, k)
            jacobianl = jacobian(i, j, k, ispec)
            integl = integl + jacobianl * weight * vec(i,j,k,ispec)
          enddo
         enddo
      enddo
    enddo

    integl_tmp = REAL(integl, CUSTOM_REAL)
    res = 0.0
    call sum_all_all_cr(integl_tmp, res)

  end subroutine Parallel_ComputeIntegral

  subroutine exit_mpi(error_msg)
    character(len=*) :: error_msg

    integer :: ier

    ! write error message to screen
    write(*,*) error_msg(1:len(error_msg))
    write(*,*) 'Error detected, aborting MPI... proc ',myrank
    call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
  end subroutine exit_mpi

  subroutine init_mpi()
    integer :: ier
    call MPI_INIT(ier)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ier)
  end subroutine init_mpi

end module global_var
