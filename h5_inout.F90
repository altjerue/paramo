module h5_inout
   use data_types
   use hdf5
   implicit none

contains

   !
   !  ###### # #      ######  ####
   !  #      # #      #      #
   !  #####  # #      #####   ####
   !  #      # #      #           #
   !  #      # #      #      #    #
   !  #      # ###### ######  ####
   !
   
   !   ----->   create file
   subroutine h5io_createf(hfname, file_id, error)
      implicit none
      character, intent(in) :: hfname*(*)
      integer, intent(inout) :: error
      integer(HID_T), intent(inout) :: file_id

      call h5fcreate_f(trim(hfname), H5F_ACC_TRUNC_F, file_id, error)

   end subroutine h5io_createf
   
   !   ----->   open file
   subroutine h5io_openf(hfname, file_id, error)
      implicit none
      
      character, intent(in) :: hfname*(*)
      integer, intent(inout) :: error
      integer(HID_T), intent(inout) :: file_id

      call h5fopen_f(trim(hfname), H5F_ACC_RDWR_F, file_id, error)

   end subroutine h5io_openf
   
   !   ----->   close file
   subroutine h5io_closef(file_id, error)
      implicit none
      integer, intent(inout) :: error
      integer(HID_T), intent(inout) :: file_id

      call h5fclose_f(file_id, error)

   end subroutine h5io_closef


   ! 
   !   ####  #####   ####  #    # #####   ####
   !  #    # #    # #    # #    # #    # #
   !  #      #    # #    # #    # #    #  ####
   !  #  ### #####  #    # #    # #####       #
   !  #    # #   #  #    # #    # #      #    #
   !   ####  #    #  ####   ####  #       ####
   ! 

   
   !   ----->   create group
   subroutine h5io_createg(file_id, hgname, group_id, error)
      implicit none
      character(len=*), intent(in) :: hgname
      integer, intent(inout) :: error
      integer(HID_T), intent(inout) :: file_id, group_id

      call h5gcreate_f(file_id, trim(hgname), group_id, error)

   end subroutine h5io_createg
   
   !   ----->   open group
   subroutine h5io_openg(file_id, hgname, group_id, error)
      implicit none
      character(len=*), intent(in) :: hgname
      integer, intent(inout) :: error
      integer(HID_T), intent(inout) :: file_id, group_id

      call h5gopen_f(file_id, trim(hgname), group_id, error)

   end subroutine h5io_openg
   
   !   ----->   close group
   subroutine h5io_closeg(group_id, error)
      implicit none
      integer, intent(inout) :: error
      integer(HID_T), intent(inout) :: group_id

      call h5gclose_f(group_id, error)

   end subroutine h5io_closeg


   ! 
   !  # #    # ##### ######  ####  ###### #####   ####
   !  # ##   #   #   #      #    # #      #    # #
   !  # # #  #   #   #####  #      #####  #    #  ####
   !  # #  # #   #   #      #  ### #      #####       #
   !  # #   ##   #   #      #    # #      #   #  #    #
   !  # #    #   #   ######  ####  ###### #    #  ####
   ! 

   !   ----->   write integer scalar
   subroutine h5io_wint0(num_id, iname, buf, error)
      implicit none
      integer, intent(in) :: buf
      integer(HID_T), intent(in) :: num_id
      character(len=*), intent(in) :: iname
      integer, intent(inout) :: error
      integer, dimension(1) :: buf_tmp, dims
      dims(1) = 1
      buf_tmp(1) = buf
      call h5io_write_int(num_id, iname, 1, dims, buf_tmp, error)
   end subroutine h5io_wint0


   !   ----->   read integer scalar
   subroutine h5io_rint0(num_id, iname, buf, error)
      implicit none
      integer(HID_T), intent(in) :: num_id
      character(len=*), intent(in) :: iname
      integer, intent(inout) :: error
      integer, intent(out) :: buf
      integer, dimension(1) :: buf_tmp, dims
      dims(1) = 1
      call h5io_read_int(num_id, iname, 1, dims, buf_tmp, error)
      buf = buf_tmp(1)
   end subroutine h5io_rint0


   !   ----->   write integer vector
   subroutine h5io_wint1(num_id, iname, buf, error)
      implicit none
      integer, intent(in), dimension(:) :: buf
      integer(HID_T), intent(in) :: num_id
      character(len=*), intent(in) :: iname
      integer, intent(inout) :: error
      integer, dimension(1) :: dims
      dims(1) = size(buf)
      call h5io_write_int(num_id, iname, 1, dims, buf, error)
   end subroutine h5io_wint1


   !   ----->   read integer 1D array
   subroutine h5io_rint1(num_id, iname, buf, error)
      implicit none
      integer(HID_T), intent(in) :: num_id
      character(len=*), intent(in) :: iname
      integer, intent(inout) :: error
      integer, intent(out), dimension(:) :: buf
      integer, dimension(1) :: dims
      dims(1) = size(buf)
      call h5io_read_int(num_id, iname, 1, dims, buf, error)
   end subroutine h5io_rint1


   !   ----->   write integer 2D array
   subroutine h5io_wint2(num_id, iname, buf, error)
      implicit none
      integer, intent(in), dimension(:, :) :: buf
      integer(HID_T), intent(in) :: num_id
      character(len=*), intent(in) :: iname
      integer, intent(inout) :: error
      integer, dimension(2) :: dims
      dims(1) = size(buf, dim=1)
      dims(2) = size(buf, dim=2)
      call h5io_write_int(num_id, iname, 2, dims, buf, error)
   end subroutine h5io_wint2


   !   ----->   read integer scalar
   subroutine h5io_rint2(num_id, iname, buf, error)
      implicit none
      integer(HID_T), intent(in) :: num_id
      character(len=*), intent(in) :: iname
      integer, intent(inout) :: error
      integer, intent(out), dimension(:, :) :: buf
      integer, dimension(2) :: dims
      dims(1) = size(buf, dim=1)
      dims(2) = size(buf, dim=2)
      call h5io_read_int(num_id, iname, 2, dims, buf, error)
   end subroutine h5io_rint2


   !   ----->   write integer in general
   subroutine h5io_write_int(num_id, iname, num_dims, dims, buf, error)
      implicit none
      integer, intent(in) :: num_dims
      integer(HID_T), intent(in) :: num_id
      integer, intent(in), dimension(num_dims) :: dims
      integer, intent(in), dimension(*) :: buf
      character(len=*), intent(in) :: iname
      integer, intent(inout) :: error
      integer(HID_T) :: ispace, iset
      integer(HSIZE_T), dimension(num_dims) :: idims
      idims(1:num_dims) = dims(1:num_dims)
      call h5screate_simple_f(num_dims, idims, ispace, error)
      call h5dcreate_f(num_id, trim(iname), H5T_NATIVE_INTEGER, ispace, iset, error)
      call h5dwrite_f(iset, H5T_NATIVE_INTEGER, buf, idims, error)
      call h5dclose_f(iset, error)
      call h5sclose_f(ispace, error)
   end subroutine h5io_write_int


   !   ----->   read integer in general
   subroutine h5io_read_int(num_id, iname, num_dims, dims, buf, error)
      implicit none
      integer, intent(in) :: num_dims
      integer(HID_T), intent(in) :: num_id
      integer, intent(in), dimension(:) :: dims
      character(len=*), intent(in) :: iname
      integer, intent(inout) :: error
      integer, intent(out), dimension(*) :: buf
      integer(HID_T) :: iset
      integer(HSIZE_T), dimension(num_dims) :: idims
      idims(1:num_dims) = dims(1:num_dims)
      call h5dopen_f(num_id, trim(iname), iset, error)
      call h5dread_f(iset, H5T_NATIVE_INTEGER, buf, idims, error)
      call h5dclose_f(iset, error)
   end subroutine h5io_read_int


   !
   !  #####   ####  #    # #####  #      ######  ####
   !  #    # #    # #    # #    # #      #      #
   !  #    # #    # #    # #####  #      #####   ####
   !  #    # #    # #    # #    # #      #           #
   !  #    # #    # #    # #    # #      #      #    #
   !  #####   ####   ####  #####  ###### ######  ####
   !

   !   ----->   write double precision scalar
   subroutine h5io_wdble0(num_id, dname, buf, error)
      implicit none
      integer(HID_T), intent(in) :: num_id
      real(dp), intent(in) :: buf
      character(len=*), intent(in) :: dname
      integer, intent(inout) :: error
      integer, dimension(1) :: dims
      real(dp), dimension(1) :: buf_tmp
      dims(1) = 1
      buf_tmp(1) = buf
      call h5io_write_double(num_id, dname, 1, dims, buf_tmp, error)
   end subroutine h5io_wdble0


   !   ----->   read double precision scalar
   subroutine h5io_rdble0(num_id, dname, buf, error)
      implicit none
      integer(HID_T), intent(in) :: num_id
      character(len=*), intent(in) :: dname
      integer, intent(inout) :: error
      real(dp), intent(out) :: buf
      integer, dimension(1) :: dims
      real(dp), dimension(1) :: buf_tmp
      dims(1) = 1
      call h5io_read_double(num_id, dname, 1, dims, buf_tmp, error)
      buf = buf_tmp(1)
   end subroutine h5io_rdble0


   !   ----->   write double precision 1D array
   subroutine h5io_wdble1(num_id, dname, buf, error)
      implicit none
      integer(HID_T), intent(in) :: num_id
      real(dp), intent(in), dimension(:) :: buf
      character(len=*), intent(in) :: dname
      integer, intent(inout) :: error
      integer, dimension(1) :: dims
      dims(1) = size(buf)
      call h5io_write_double(num_id, dname, 1, dims, buf, error)
   end subroutine h5io_wdble1


   !   ----->   read double precision 1D array
   subroutine h5io_rdble1(num_id, dname, buf, error)
      implicit none
      integer(HID_T), intent(in) :: num_id
      character(len=*), intent(in) :: dname
      integer, intent(inout) :: error
      real(dp), intent(out), dimension(:) :: buf
      integer, dimension(1) :: dims
      dims(1) = size(buf)
      call h5io_read_double(num_id, dname, 1, dims, buf, error)
   end subroutine h5io_rdble1


   !   ----->   write double precision 2D array
   subroutine h5io_wdble2(num_id, dname, buf, error)
      implicit none
      integer(HID_T), intent(in) :: num_id
      real(dp), intent(in), dimension(:, :) :: buf
      character(len=*), intent(in) :: dname
      integer, intent(inout) :: error
      integer, dimension(2) :: dims
      dims(1) = size(buf, dim=1)
      dims(2) = size(buf, dim=2)
      call h5io_write_double(num_id, dname, 2, dims, buf, error)
   end subroutine h5io_wdble2


   !   ----->   read double precision 2D array
   subroutine h5io_rdble2(num_id, dname, buf, error)
      implicit none
      integer(HID_T), intent(in) :: num_id
      character(len=*), intent(in) :: dname
      integer, intent(inout) :: error
      real(dp), intent(out), dimension(:, :) :: buf
      integer, dimension(2) :: dims
      dims(1) = size(buf, dim=1)
      dims(2) = size(buf, dim=2)
      call h5io_read_double(num_id, dname, 2, dims, buf, error)
   end subroutine h5io_rdble2


   !   ----->   write double precision in general
   subroutine h5io_write_double(num_id, dname, num_dims, dims, buf, error)
      implicit none
      integer, intent(in) :: num_dims
      integer(HID_T), intent(in) :: num_id
      integer, intent(in), dimension(num_dims) :: dims
      real(dp), intent(in), dimension(*) :: buf
      character(len=*), intent(in) :: dname
      integer, intent(inout) :: error
      integer(HID_T) :: dspace, dset
      integer(HSIZE_T), dimension(num_dims) :: ddims
      ddims(1:num_dims) = dims(1:num_dims)
      call h5screate_simple_f(num_dims, ddims, dspace, error)
      call h5dcreate_f(num_id, trim(dname), H5T_NATIVE_DOUBLE, dspace, dset, error)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, buf, ddims, error)
      call h5dclose_f(dset, error)
      call h5sclose_f(dspace, error)
   end subroutine h5io_write_double


   !   ----->   read double precision in general
   subroutine h5io_read_double(num_id, dname, num_dims, dims, buf, error)
      implicit none
      integer, intent(in) :: num_dims
      integer(HID_T), intent(in) :: num_id
      integer, intent(in), dimension(num_dims) :: dims
      character(len=*), intent(in) :: dname
      integer, intent(inout) :: error
      real(dp), intent(out), dimension(*) :: buf
      integer(HID_T) :: dset
      integer(HSIZE_T), dimension(num_dims) :: ddims
      ddims(1:num_dims) = dims(1:num_dims)
      call h5dopen_f(num_id, trim(dname), dset, error)
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, buf, ddims, error)
      call h5dclose_f(dset, error)
   end subroutine h5io_read_double

end module h5_inout
