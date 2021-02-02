      subroutine mkmov_hdf_ycut

      use param
      use hdf5
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only: q2,q3,q1,dsal,co2,dens !,pr

      IMPLICIT none

      integer im,jm,km
      integer ic,jc,kc
      integer ip,jp,kp

      integer hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_q1v
      integer(HID_T) :: dset_q2v
      integer(HID_T) :: dset_q3v
      !integer(HID_T) :: dset_prv
      integer(HID_T) :: dset_sal
      integer(HID_T) :: dset_den
      integer(HID_T) :: dset_co2

      integer(HSIZE_T) :: dims(2)


      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 

      integer :: comm, info
      integer :: ndims

      real tprfi
      real, dimension(:,:), allocatable :: prx,v1,v2,v3,densx,dsalx,co2xx
      integer itime

      character(70) namfile,xdmnam
      character(5) ipfi

      !allocate(prx(n1m,n3m))
      allocate(v1(n1m,n3m))
      allocate(v2(n1m,n3m))
      allocate(v3(n1m,n3m))
      allocate(densx(n1m,n3m))
      allocate(dsalx(n1m,n3m))
      allocate(co2xx(n1m,n3m))

      call h5open_f(hdf_error)
      !RO   first of all, fill up velocity arrays
      !RO   so that we have node-centred velocities

      ndims=2
      jc=FLOOR(0.55*dx2) !n2m/2
      jp=jc+1
      
      do kc=kstart,kend
       kp=kc+1

       do ic=1,n1m                                                   
        ip=ic+1

        v1(ic,kc) = (q1(ic,jc,kc)+q1(ip,jc,kc))*0.5d0
        v2(ic,kc) = (q2(ic,jc,kc)+q2(ic,jp,kc))*0.5d0
        v3(ic,kc) = (q3(ic,jc,kc)+q3(ic,jc,kp))*0.5d0
        !prx(ic,kc)= pr(ic,jc,kc)
        dsalx(ic,kc) = dsal(ic,jc,kc)
        densx(ic,kc) = dens(ic,jc,kc)
        co2xx(ic,kc) = co2(ic,jc,kc)

       end do
      end do

      !RO   File writing part

      !RO   Form the name of the file

      tprfi = 1.d0/tframe
      itime=nint(time*tprfi)
      write(ipfi,'(i5.5)')itime

      namfile='flowmov/frame_'//ipfi//'_ycut.h5'
      xdmnam='flowmov/frame_'//ipfi//'_ycut.xmf'

      !RO   Sort out MPI definitions and open file

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error,&
     &                 access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      !RO   Create dataspace

      dims(1)=n1m
      dims(2)=n3m
      call h5screate_simple_f(ndims, dims, &
     &                        filespace, hdf_error)
         

      !RO   Create the dataset with default properties.
      call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, filespace,dset_q1v, hdf_error)
      call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE,filespace,dset_q2v, hdf_error)
      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE,filespace,dset_q3v, hdf_error)
      !call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE,filespace,dset_prv, hdf_error)
      call h5dcreate_f(file_id, 'Sal', H5T_NATIVE_DOUBLE,filespace,dset_sal, hdf_error)
      call h5dcreate_f(file_id, 'Den', H5T_NATIVE_DOUBLE,filespace,dset_den, hdf_error)
      call h5dcreate_f(file_id, 'CO2', H5T_NATIVE_DOUBLE,filespace,dset_co2, hdf_error)

      !RO   Set offsets and element counts
      data_count(1)=n1m
      data_count(2)=kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = kstart-1

      !RO   Create dataspace in memory
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      !RO   Select hyperslab  and then write it
      call h5dget_space_f(dset_q1v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
       call h5dwrite_f(dset_q1v, H5T_NATIVE_DOUBLE,&
     &   v1(1:n1m,kstart:kend), dims, &
     &   hdf_error, file_space_id = filespace, mem_space_id = memspace, &
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dget_space_f(dset_q2v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
       call h5dwrite_f(dset_q2v, H5T_NATIVE_DOUBLE,&
     &   v2(1:n1m,kstart:kend), dims, &
     &   hdf_error, file_space_id = filespace, mem_space_id = memspace, &
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dget_space_f(dset_q3v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
       call h5dwrite_f(dset_q3v, H5T_NATIVE_DOUBLE, &
     &   v3(1:n1m,kstart:kend), dims, &
     &   hdf_error, file_space_id = filespace, mem_space_id = memspace, &
     &   xfer_prp = plist_id)

!      call h5dget_space_f(dset_prv, filespace, hdf_error)
!      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
!      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
!       call h5dwrite_f(dset_prv, H5T_NATIVE_DOUBLE, &
!     &   prx(1:n1m,kstart:kend), dims, &
!     &   hdf_error, file_space_id = filespace, mem_space_id = memspace, &
!     &   xfer_prp = plist_id)

      call h5dget_space_f(dset_sal, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
       call h5dwrite_f(dset_sal, H5T_NATIVE_DOUBLE, &
     &   dsalx(1:n1m,kstart:kend), dims, &
     &   hdf_error, file_space_id = filespace, mem_space_id = memspace, &
     &   xfer_prp = plist_id)

      !Dens
      call h5dget_space_f(dset_den, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
      call h5dwrite_f(dset_den, H5T_NATIVE_DOUBLE, densx(1:n1m,kstart:kend), dims,  hdf_error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      !CO2
      call h5dget_space_f(dset_co2, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
      call h5dwrite_f(dset_co2, H5T_NATIVE_DOUBLE, co2xx(1:n1m,kstart:kend), dims,  hdf_error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      !RO   Close properties and file
      call h5dclose_f(dset_q1v, hdf_error)
      call h5dclose_f(dset_q2v, hdf_error)
      call h5dclose_f(dset_q3v, hdf_error)
      !call h5dclose_f(dset_prv, hdf_error)
      call h5dclose_f(dset_sal, hdf_error)
      call h5dclose_f(dset_den, hdf_error)
      call h5dclose_f(dset_co2, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5pclose_f(plist_id, hdf_error)
      call h5fclose_f(file_id, hdf_error)
       
      if (myid.eq.0) then

      open(45,file=xdmnam,status='unknown')
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""ycut"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""2DSMesh"" NumberOfElements=&
     &""",i4," ",i4,"""/>")') n3m,n1m
      write(45,'("<Geometry GeometryType=""X_Y"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("cordin_ycut_info.h5:/x")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("cordin_ycut_info.h5:/z")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""X-Velocity""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_",i5.5,"_ycut.h5:/Vx")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Y-Velocity""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_",i5.5,"_ycut.h5:/Vy")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Z-Velocity""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_",i5.5,"_ycut.h5:/Vz")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
!      write(45,'("<Attribute Name=""Pressure""&
!     & AttributeType=""Scalar"" Center=""Node"">")')
!      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
!     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
!      write(45,'("frame_",i5.5,"_ycut.h5:/Pr")') itime
!      write(45,'("</DataItem>")')
!      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Salinity""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_",i5.5,"_ycut.h5:/Sal")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Temperature""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_",i5.5,"_ycut.h5:/Den")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""CO2""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_",i5.5,"_ycut.h5:/CO2")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')

      write(45,'("<Time Value=""",e12.5""" />")')time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45)

      end if


      call h5close_f(hdf_error)

      !if (allocated(prx)) deallocate(prx)
      if (allocated(v1)) deallocate(v1)
      if (allocated(v2)) deallocate(v2)
      if (allocated(v3)) deallocate(v3)
      if (allocated(densx)) deallocate(densx)
      if (allocated(dsalx)) deallocate(dsalx)
      if (allocated(co2xx)) deallocate(co2xx)

      return
      end
