      subroutine inimov_hdf_xcut

      use param
      use mpih
      use hdf5

      IMPLICIT none
      integer j,k

      integer hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: cordin_dset_rid
      integer(HID_T) :: cordin_dset_zid
      integer(HID_T) :: cordin_dspace_id
      integer(HSIZE_T) :: dims(2)

      real ymov(m2m,m3m),zmov(m2m,m3m)
      character(70) namfile

      call h5open_f(hdf_error)

      !RO   Fill up the arrays for the constant z cut

      do k=1,n3m
        do j=1,n2m
          ymov(j,k) = ym(j)
          zmov(j,k) = zm(k)
        end do
      end do

      !RO   Initialize hdf file

      namfile='flowmov/cordin_xcut_info.h5'

      dims(1)=n2m
      dims(2)=n3m


      if (myid.eq.0) then
      call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error)

      !RO   Create dataspace
      call h5screate_simple_f(2, dims, cordin_dspace_id, hdf_error)


      !RO   Create datasets to write to and write
      call h5dcreate_f(file_id, 'y', H5T_NATIVE_DOUBLE,cordin_dspace_id,cordin_dset_rid, hdf_error)
      call h5dwrite_f(cordin_dset_rid, H5T_NATIVE_DOUBLE,ymov, dims, hdf_error)
      call h5dcreate_f(file_id, 'z', H5T_NATIVE_DOUBLE,cordin_dspace_id,cordin_dset_zid, hdf_error)
      call h5dwrite_f(cordin_dset_zid, H5T_NATIVE_DOUBLE, zmov, dims, hdf_error)

      !RO   Close datasets and dataspace

      call h5dclose_f(cordin_dset_rid, hdf_error)
      call h5dclose_f(cordin_dset_zid, hdf_error)
      call h5sclose_f(cordin_dspace_id, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      endif

      call h5close_f(hdf_error)

      return                                                            
      end                                                               


