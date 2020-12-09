module AdiosIO

  use mpi
  use adios_read_mod
  use adios_write_mod
  use adios_helpers_mod

  use global_var, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, nspec, &
    ADIOS_TRANSPORT_METHOD, ADIOS_BUFFER_SIZE_IN_MB, &
    myrank, nprocs

  implicit none

  private

  public :: read_bp_file_int
  public :: read_bp_file_real
  public :: write_bp_file
  public :: calculate_jacobian_matrix

  interface read_bp_file_int
    module procedure read_bp_file_int_1d
    module procedure read_bp_file_int_4d
  end interface read_bp_file_int

  interface read_bp_file_real
    module procedure read_bp_file_real_1d
    module procedure read_bp_file_real_4d
    module procedure read_bp_file_real_5d
  end interface read_bp_file_real

  contains

  subroutine read_bp_file_int_1d(bp_file, varname, data_array)
    ! Used to read one variable(integer) with dimension of 1
    ! For example, the variable idoubling: (NSPEC)
    character(len=*), intent(in) :: bp_file
    character(len=*), intent(in):: varname
    integer, dimension(:), intent(inout):: data_array

    ! local variables
    integer :: ier, local_dim
    integer(kind=8) :: read_handle
    integer(kind=8) :: sel
    integer(kind=8), dimension(1) :: start, count_ad

    call adios_read_open_file(read_handle, trim(bp_file), 0, &
                              MPI_COMM_WORLD, ier)

    call adios_get_scalar(read_handle, trim(varname)//"/local_dim",&
                          local_dim, ier)
    start(1) = local_dim * myrank
    count_ad(1) = local_dim

    call adios_selection_boundingbox(sel, 1, start, count_ad)
    call adios_schedule_read( &
            read_handle, sel, trim(varname)//"/array", &
            0, 1, data_array, ier)

    call adios_perform_reads(read_handle, ier)
    call adios_read_close(read_handle,ier)
  end subroutine read_bp_file_int_1d

  subroutine read_bp_file_int_4d(bp_file, varname, data_array)
    ! Used to read one variable(integer) with dimension of 4.
    ! For example, the variable ibool: (NGLLX, NGLLY, NGLLZ, NSPEC)
    character(len=*), intent(in) :: bp_file
    character(len=*), intent(in):: varname
    integer, dimension(:, :, :, :), intent(inout):: data_array

    ! local variables
    integer :: ier, local_dim
    integer(kind=8) :: read_handle
    integer(kind=8) :: sel
    integer(kind=8), dimension(1) :: start, count_ad

    call adios_read_open_file(read_handle, trim(bp_file), 0, &
                              MPI_COMM_WORLD, ier)

    call adios_get_scalar(read_handle, trim(varname)//"/local_dim",&
                          local_dim, ier)
    start(1) = local_dim * myrank
    count_ad(1) = local_dim

    call adios_selection_boundingbox(sel, 1, start, count_ad)
    call adios_schedule_read( &
            read_handle, sel, trim(varname)//"/array", &
            0, 1, data_array, ier)

    call adios_perform_reads(read_handle, ier)
    call adios_read_close(read_handle,ier)
  end subroutine read_bp_file_int_4d

  subroutine read_bp_file_real_1d(bp_file, varname, data_array)
    ! Used to read one variable(real) with dimension of 1
    ! For example, the coordinate -- x_glob, y_glob and z_glob: (NGLOB)
    character(len=*), intent(in) :: bp_file
    character(len=*), intent(in):: varname
    real(kind=CUSTOM_REAL), dimension(:), intent(inout):: data_array

    ! local variables
    integer :: ier, local_dim
    integer(kind=8) :: read_handle
    integer(kind=8) :: sel
    integer(kind=8), dimension(1) :: start, count_ad

    data_array = 0.0

    call adios_read_open_file(read_handle, trim(bp_file), 0, &
                              MPI_COMM_WORLD, ier)

    call adios_get_scalar(read_handle, trim(varname)//"/local_dim",&
                          local_dim, ier)
    start(1) = local_dim * myrank
    !count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec
    count_ad(1) = local_dim

    call adios_selection_boundingbox(sel, 1, start, count_ad)
    call adios_schedule_read( &
            read_handle, sel, trim(varname)//"/array", &
            0, 1, data_array, ier)

    call adios_perform_reads(read_handle, ier)
    call adios_read_close(read_handle,ier)

  end subroutine read_bp_file_real_1d

  subroutine read_bp_file_real_4d(bp_file, varname, data_array)
    ! Mainly used to read in multiple kernels(real) or models(vpv, vsv, ...) at once
    ! Dimension of data_array: (NGLLX, NGLLY, NGLLZ, NSPEC, nkernels)
    character(len=*), intent(in) :: bp_file
    character(len=*), intent(in):: varname
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(inout):: data_array

    ! local variables
    integer :: ier, local_dim
    integer(kind=8) :: read_handle
    integer(kind=8) :: sel
    integer(kind=8), dimension(1) :: start, count_ad

    data_array = 0.0

    call adios_read_open_file(read_handle, trim(bp_file), 0, &
                              MPI_COMM_WORLD, ier)

    call adios_get_scalar(read_handle, trim(varname)//"/local_dim", &
                          local_dim, ier)
    start(1) = local_dim * myrank
    !count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec
    count_ad(1) = local_dim

    call adios_selection_boundingbox(sel, 1, start, count_ad)
    call adios_schedule_read( &
            read_handle, sel, trim(varname)//"/array", &
            0, 1, data_array, ier)

    call adios_perform_reads(read_handle, ier)
    call adios_read_close(read_handle,ier)

  end subroutine read_bp_file_real_4d

  subroutine read_bp_file_real_5d(bp_file, varnames, data_array)
    ! Mainly used to read in multiple kernels(real) or models(vpv, vsv, ...) at once
    ! Dimension of data_array: (NGLLX, NGLLY, NGLLZ, NSPEC, nkernels)
    character(len=*), intent(in) :: bp_file
    character(len=*), dimension(:), intent(in):: varnames
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout):: data_array

    ! local variables
    integer :: iker, ier, local_dim
    integer :: nvars
    integer(kind=8) :: read_handle
    integer(kind=8) :: sel
    integer(kind=8), dimension(1) :: start, count_ad

    data_array = 0.0

    call adios_read_open_file(read_handle, trim(bp_file), 0, &
                              MPI_COMM_WORLD, ier)

    nvars = size(varnames)
    do iker = 1, nvars
      call adios_get_scalar(read_handle, trim(varnames(iker))//"/local_dim",&
                            local_dim, ier)
      start(1) = local_dim * myrank
      !count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec
      count_ad(1) = local_dim

      call adios_selection_boundingbox(sel, 1, start, count_ad)
      call adios_schedule_read( &
              read_handle, sel, trim(varnames(iker))//"/array", &
              0, 1, data_array(:, :, :, :, iker), ier)
      !call adios_perform_reads(read_handle, ier)
      !call adios_selection_delete(sel)
    enddo

    call adios_perform_reads(read_handle, ier)
    call adios_read_close(read_handle,ier)

  end subroutine read_bp_file_real_5d

  subroutine write_bp_file(data_array, varnames, group_name, filename)
    implicit none
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: data_array
    character(len=*), dimension(:), intent(in) :: varnames
    character(len=*), intent(in) :: group_name
    character(len=*), intent(in) :: filename

    integer(kind=8) :: group, handle, groupsize, totalsize
    integer :: local_dim
    integer :: ier, iker
    integer :: nvars

    nvars = size(varnames)
    call MPI_Barrier(MPI_COMM_WORLD, ier)

    groupsize = 0

    call adios_init_noxml (MPI_COMM_WORLD, ier);
    call adios_set_max_buffer_size(ADIOS_BUFFER_SIZE_IN_MB)
    call adios_declare_group(group, group_name, "", 1, ier)
    call adios_select_method(group, ADIOS_TRANSPORT_METHOD, "", "", ier)

    ! Define Variables
    call define_adios_scalar(group, groupsize, "", "nspec", nspec)
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    do iker=1, nvars
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                       "", trim(varnames(iker)), &
                                       data_array(:,:,:,:,iker))
    enddo

    call MPI_Barrier(MPI_COMM_WORLD, ier)

    ! Write Data
    call adios_open(handle, group_name, filename, "w", &
                    MPI_COMM_WORLD, ier)
    call adios_group_size(handle, groupsize, totalsize, ier)

    call adios_write(handle, "nspec", nspec, ier)
    do iker=1, nvars
      call write_adios_global_1d_array(handle, myrank, nprocs, local_dim, &
                                       trim(varnames(iker)), &
                                       data_array(:,:,:,:,iker))
    enddo

    call adios_close(handle, ier)
    call MPI_Barrier(MPI_COMM_WORLD, ier)
    call adios_finalize(myrank, ier)

  end subroutine write_bp_file

  subroutine calculate_jacobian_matrix(solver_file, jacobian)
    use global_var, only : CUSTOM_REAL, build_jacobian
    character(len=*), intent(in) :: solver_file
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(out) :: jacobian

    ! Don't change the order of varlist, since the order matters when
    ! calling the build_jacobian() function
    character(len=100), parameter :: regionTag = "reg1"
    integer, parameter :: nvars = 9
    character(len=100), dimension(9) :: varlist = &
      (/character(len=100) :: "xixstore", "xiystore", "xizstore", &
                              "etaxstore", "etaystore", "etazstore", &
                              "gammaxstore", "gammaystore", "gammazstore"/)

    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), allocatable :: model
    integer, dimension(:, :, :, :), allocatable :: ibool
    integer :: i

    do i = 1, nvars
     varlist(i) = trim(regionTag)//"/"//trim(varlist(i))
    enddo

    allocate(model(NGLLX, NGLLY, NGLLZ, NSPEC, nvars))
    call read_bp_file_real(solver_file, varlist, model)

    ! read ibool
    allocate(ibool(NGLLX, NGLLY, NGLLZ, NSPEC))
    call read_bp_file_int(solver_file, "reg1/ibool", ibool)

    call build_jacobian( &
           ibool(:,:,:,:), &
           model(:,:,:,:,1), model(:,:,:,:,2), model(:,:,:,:,3), &
           model(:,:,:,:,4), model(:,:,:,:,5), model(:,:,:,:,6), &
           model(:,:,:,:,7), model(:,:,:,:,8), model(:,:,:,:,9), &
           jacobian)

  end subroutine calculate_jacobian_matrix

end module AdiosIO
