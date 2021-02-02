!================================================
      subroutine mpi_write_continua
      use param
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      use local_arrays, only: dens,q2,q3,q1,dsal,isbody1,isbody2,co2
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_q1
      integer(HID_T) :: dset_q2
      integer(HID_T) :: dset_q3
      integer(HID_T) :: dset_dens
      integer(HID_T) :: dset_dsal
      integer(HID_T) :: dset_body
      integer(HID_T) :: dset_co2

      integer(HID_T) :: plist_id

      integer(HSIZE_T), dimension(3) :: dims
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character filnam1*30
      character filnam2*30
      character filnam3*30
      character filnam4*30
      character filnam5*30
      character filnam6*30
      character filnam7*30
      character filnam8*30
      character filnam9*30
      character filnam10*30
      character filnam*50

      call h5open_f(hdf_error)

      ! Sort out MPI definitions
      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      ! Form the name of the file

      filnam1 = 'continua_dens.h5'
      filnam2 = 'continua_q1.h5'
      filnam3 = 'continua_q2.h5'
      filnam4 = 'continua_q3.h5'
      filnam5 = 'continua_dsal.h5'
      filnam6 = 'continua_pr.h5'
      filnam8 = 'continua_co2.h5'
      filnam9 = 'continua_body1.h5'
      filnam10 = 'continua_body2.h5'

      ! Set offsets and element counts
   
      ndims = 3

      dims(1)=n1
      dims(2)=n2
      dims(3)=n3m

      data_count(1) = n1
      data_count(2) = n2
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

      ! dens

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'dens', H5T_NATIVE_DOUBLE,filespace, dset_dens, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_dens, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_dens, H5T_NATIVE_DOUBLE,dens(1:n1,1:n2,kstart:kend), dims, hdf_error, file_space_id = slabspace,mem_space_id = memspace,xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dens, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      ! q1
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vth', H5T_NATIVE_DOUBLE,filespace, dset_q1, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_q1, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_q1, H5T_NATIVE_DOUBLE,q1(1:n1,1:n2,kstart:kend), dims, hdf_error, file_space_id = slabspace,mem_space_id = memspace, xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q1, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      ! q2
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vr', H5T_NATIVE_DOUBLE,filespace, dset_q2, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_q2, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_q2, H5T_NATIVE_DOUBLE,q2(1:n1,1:n2,kstart:kend), dims,hdf_error, file_space_id = slabspace,mem_space_id = memspace, xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q2, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      ! q2
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE,filespace, dset_q3, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_q3, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_q3, H5T_NATIVE_DOUBLE,q3(1:n1,1:n2,kstart:kend), dims, hdf_error, file_space_id = slabspace,mem_space_id = memspace, xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q3, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      ! isbody1
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam9, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'isbody1', H5T_NATIVE_DOUBLE,filespace, dset_body, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_body, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_body, H5T_NATIVE_DOUBLE,isbody1(1:n1,1:n2,kstart:kend), dims, hdf_error, file_space_id = slabspace,mem_space_id = memspace, xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_body, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      ! isbody2
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam10, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'isbody2', H5T_NATIVE_DOUBLE,filespace, dset_body, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_body, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_body, H5T_NATIVE_DOUBLE,isbody2(1:n1,1:n2,kstart:kend), dims, hdf_error, file_space_id = slabspace,mem_space_id = memspace, xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_body, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      ! co2
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam8, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'co2', H5T_NATIVE_DOUBLE,filespace, dset_co2, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_co2, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_co2, H5T_NATIVE_DOUBLE,co2(1:n1,1:n2,kstart:kend), dims, hdf_error, file_space_id = slabspace,mem_space_id = memspace, xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_co2, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
      ! Set offsets and element counts

      ndims = 3
      dims(1)=n1r
      dims(2)=n2r
      dims(3)=n3mr

      data_count(1) = n1r
      data_count(2) = n2r
      data_count(3) = kendr-kstartr+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstartr-1

      ! dsal

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam5, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'dsal', H5T_NATIVE_DOUBLE,filespace, dset_dsal, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_dsal, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_dsal, H5T_NATIVE_DOUBLE,dsal(1:n1r,1:n2r,kstartr:kendr), dims,hdf_error, file_space_id = slabspace,mem_space_id = memspace,xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dsal, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
   
      call h5close_f(hdf_error)


      ! continua_grid
      if (myid .eq. 0) then
        filnam   =trim('./continua_grid.h5')
        call HdfStart
        call HdfCreateBlankFile(filnam)
        call HdfSerialWriteIntScalar('n1',filnam,n1)
        call HdfSerialWriteIntScalar('n2',filnam,n2)
        call HdfSerialWriteIntScalar('n3',filnam,n3)
        call HdfSerialWriteRealScalar('rext1',filnam,rext1)
        call HdfSerialWriteRealScalar('rext2',filnam,rext2)
        call HdfSerialWriteRealScalar('time',filnam,time)
        call HdfSerialWriteIntScalar('istr3',filnam,istr3)
        call HdfSerialWriteIntScalar('str3',filnam,str3)
        call HdfSerialWriteIntScalar('mref1',filnam,mref1)
        call HdfSerialWriteIntScalar('mref2',filnam,mref2)
        call HdfSerialWriteIntScalar('mref3',filnam,mref3)
        call HdfClose
      endif

      return
      end subroutine mpi_write_continua

!================================================      
      subroutine mpi_read_continua(n1o,n2o,n3o,ks,ke,intvar,qua)
      use mpih
      use param
      use hdf5
      implicit none
      integer, intent(in) :: ks,ke,n2o,n1o,n3o
      real, dimension(1:n1o,1:n2o,ks-lvlhalo:ke+lvlhalo),intent(out)::qua
      integer k,j,i

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_qua

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      integer, intent(in) :: intvar
      character(70) :: filnam1
      character(10) :: dsetname

      call h5open_f(hdf_error)

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      ! Select file and dataset based on intvar

      select case (intvar)
        case (1)
          dsetname = trim('Vth')
          filnam1 = trim('continua_q1.h5')
        case (2)
          dsetname = trim('Vr')
          filnam1 = trim('continua_q2.h5')
        case (3)
          dsetname = trim('Vz')
          filnam1 = trim('continua_q3.h5')
        case (4)
          dsetname = trim('dens')
          filnam1 = trim('continua_dens.h5')
        case (5)
          dsetname = trim('dsal')
          filnam1 = trim('continua_dsal.h5')
        case (6)
          dsetname = trim('dens_in')
          filnam1 = trim('continua_dens_in.h5')
        case (7)
          dsetname = trim('dsal_in')
          filnam1 = trim('continua_dsal_in.h5')
        case (8)
          dsetname = trim('co2')
          filnam1 = trim('continua_co2.h5')
        case (9)
          dsetname = trim('isbody1')
          filnam1 = trim('continua_body1.h5')
        case (10)
          dsetname = trim('isbody2')
          filnam1 = trim('continua_body2.h5')
      end select

      do k=ks,ke
        do j=1,n2o
          do i=1,n1o
            qua(i,j,k)=0.d0
          enddo
        enddo
      enddo

      ! Set offsets and element counts
   
      ndims = 3

      dims(1)=n1o
      dims(2)=n2o
      dims(3)=n3o-1


      data_count(1) = n1o
      data_count(2) = n2o
      data_count(3) = ke-ks+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = ks-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id,hdf_error)

      call h5dopen_f(file_id, dsetname, dset_qua, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_qua, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE,qua(1:n1o,1:n2o,ks:ke), dims, hdf_error, file_space_id = slabspace,mem_space_id = memspace, xfer_prp = plist_id)

      call h5dclose_f(dset_qua, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5pclose_f(plist_id, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      call h5close_f(hdf_error)

      if(myid.eq.0)write(*,'(5x,a)')'reading complete: '//filnam1

      return
      end subroutine mpi_read_continua
