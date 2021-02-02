!===========================================================
! solve wave equations
!***********************************************************
      subroutine calcOutflowBCs
      use param
      use local_arrays, only: q1,q2,q3,dsal,dens,co2
      use mpi_param, only: kstart,kend,kstartr,kendr
      use mpih
      use outflow_vars
      use inflow_vars
      use walls_vars
      implicit none
      integer :: i,j,k
      real :: d1,d2,d3
      real :: ds,dde,dco2
      real :: cou
      real :: qout,qin,area,darea,cor

      cou=0.3d0

      !-- Inflow
      do k=kstart,kend
        if(zm(k)<Inzlevel+InOutzlen/2.0 .and. zm(k)>Inzlevel-InOutzlen/2.0) then
          do i=1,n1
            dq1bw(i,k)=0.d0-q1bw(i,k)
            dq2bw(i,k)=InOutvel-q2bw(i,k)
            dq3bw(i,k)=0.d0-q3bw(i,k)
            ddensbw(i,k)=0.d0-densbw(i,k)  !Fixed inflow value
            dco2bw(i,k)=0.d0-co2bw(i,k)    !Fixed inflow value
          enddo
        else
          dq1bw(1:n1,k)=0.d0-q1bw(1:n1,k)
          dq2bw(1:n1,k)=0.d0-q2bw(1:n1,k)
          dq3bw(1:n1,k)=0.d0-q3bw(1:n1,k)
          ddensbw(1:n1,k)=dens(1:n1,1,k)-densbw(1:n1,k)  !Adiabatic BC
          dco2bw(1:n1,k)=co2(1:n1,1,k)-co2bw(1:n1,k)     !Adiabatic BC
        endif
      enddo
      do k=kstartr,kendr
        if(zmr(k)<Inzlevel+InOutzlen/2.0 .and. zmr(k)>Inzlevel-InOutzlen/2.0) then
          do i=1,n1r
             ddsalbw(i,k)=0.d0-dsalbw(i,k) !Fixed inflow value
          enddo
        else
          ddsalbw(1:n1r,k)=dsal(1:n1r,1,k)-dsalbw(1:n1r,k)   !Adiabatic BC
        endif
      enddo

      ! Radiative B.C. at outflow (j=n2m)
      ! Courant number gives speed of waves.
      do k=kstart,kend
        if(zm(k)<Outzlevel+InOutzlen/2.0 .and. zm(k)>Outzlevel-InOutzlen/2.0) then
          do i=1,n1
            !--  radiative b.c. for x2 velocity at the outflow  
            d1=(q1be(i,k)-q1(i,n2m,k))*dx2
            dq1be(i,k)=-dt*(ga*d1+ro*dq1beo(i,k))*cou

            !--  radiative b.c. for x2 velocity at the outflow  
            d2=(q2be(i,k)-q2(i,n2m,k))*dx2
            dq2be(i,k)=-dt*(ga*d2+ro*dq2beo(i,k))*cou

            !--  radiative b.c. for x2 velocity at the outflow  
            d3=(q3be(i,k)-q3(i,n2m,k))*dx2
            dq3be(i,k)=-dt*(ga*d3+ro*dq3beo(i,k))*cou

            !--  radiative b.c. for dens
            dde=(densbe(i,k)-dens(i,n2m,k))*2.*dx2
            ddensbe(i,k)=-dt*(ga*dde+ro*ddensbeo(i,k))*cou

            !--  radiative b.c. for co2
            dco2=(co2be(i,k)-co2(i,n2m,k))*2.*dx2
            dco2be(i,k)=-dt*(ga*dco2+ro*dco2beo(i,k))*cou

            !--  save data for next time step
            dq1beo(i,k)=d1
            dq2beo(i,k)=d2
            dq3beo(i,k)=d3
            ddensbeo(i,k)=dde
            dco2beo(i,k)=dco2
          enddo
        else
          dq1be(1:n1,k)=0.d0-q1be(1:n1,k)
          dq2be(1:n1,k)=0.d0-q2be(1:n1,k)
          dq3be(1:n1,k)=0.d0-q3be(1:n1,k)
          ddensbe(1:n1,k)=dens(1:n1,n2m,k)-densbe(1:n1,k)   !Adiabatic BC
          dco2be(1:n1,k)=co2(1:n1,n2m,k)-co2be(1:n1,k)      !Adiabatic BC
        endif
      enddo

      do k=kstartr,kendr
        if(zmr(k)<Outzlevel+InOutzlen/2.0 .and. zmr(k)>Outzlevel-InOutzlen/2.0) then
          do i=1,n1r
            !--  radiative b.c. for dsal at the outflow  
            ds=(dsalbe(i,k)-dsal(i,n2mr,k))*2.*dx2
            ddsalbe(i,k)=-dt*(ga*ds+ro*ddsalbeo(i,k))*cou
            
            !--  save data for next time step
            ddsalbeo(i,k)=ds
          enddo
        else
          ddsalbe(1:n1r,k)=dsal(1:n1r,n2mr,k)-dsalbe(1:n1r,k)   !Adiabatic BC
        endif
      enddo 

      !-- Calculate qout
      area=0.
      qout=0.
      qin =0.
      do k=kstart,kend
        if(zm(k)<Inzlevel+InOutzlen/2.0 .and. zm(k)>Inzlevel-InOutzlen/2.0) then
        do i=1,n1m
          darea=1.0d0/dx3/dx1
          qin =qin +dq2bw(i,k)*darea
        end do
        endif
        if(zm(k)<Outzlevel+InOutzlen/2.0 .and. zm(k)>Outzlevel-InOutzlen/2.0) then
        do i=1,n1m
          darea=1.0d0/dx3/dx1
          area = area+darea
          qout=qout+dq2be(i,k)*darea
        end do
        endif
      enddo
      call mpi_globalsum_double_arr(area,1)
      call mpi_globalsum_double_arr(qout,1)
      call mpi_globalsum_double_arr(qin,1)

      !-- Correct qout
      cor=(qin-qout)/area
      do k=kstart,kend
        if(zm(k)<Outzlevel+InOutzlen/2.0 .and. zm(k)>Outzlevel-InOutzlen/2.0) then
        do i=1,n1m
          dq2be(i,k)=dq2be(i,k)+cor
        end do
        endif
      end do

      return
      end

!===========================================================
! update outflow BCs
!***********************************************************
      subroutine updateOutflowBCs
      use param
      use outflow_vars
      use mpi_param, only: kstart,kend,kstartr,kendr
      use walls_vars, only: q1ybn,q2ybn,q3ybn,dsalybn,densybn,co2ybn
      implicit none
      real :: darea
      integer :: i,k

      !-- outflowBC
      do k=kstart,kend
        do i=1,n1
          ! vel
          q1be(i,k) = q1be(i,k)+dq1be(i,k)
          q2be(i,k) = q2be(i,k)+dq2be(i,k)
          q3be(i,k) = q3be(i,k)+dq3be(i,k)
          q1ybn(i,k) = q1be(i,k)
          q2ybn(i,k) = q2be(i,k)
          q3ybn(i,k) = q3be(i,k)

          ! scalar
          densbe(i,k) = densbe(i,k)+ddensbe(i,k)
          co2be(i,k) = co2be(i,k)+dco2be(i,k)
          densybn(i,k) = densbe(i,k)
          co2ybn(i,k) = co2be(i,k)
        enddo
      enddo

      do k=kstartr,kendr
        do i=1,n1r
          dsalbe(i,k) = dsalbe(i,k)+ddsalbe(i,k)
          dsalybn(i,k) = dsalbe(i,k)
        enddo
      enddo

      return
      end

!===========================================================
! update inflow BCs
!***********************************************************
      subroutine updateInflowBCs
      use param
      use mpi_param, only: kstart,kend,kstartr,kendr
      use inflow_vars
      use walls_vars, only: q1ybs,q2ybs,q3ybs,dsalybs,densybs,co2ybs
      implicit none
      real :: darea
      integer :: i,k

      !-- inflowBC
      do k=kstart,kend
        do i=1,n1
          ! vel
          q1bw(i,k) = q1bw(i,k)+dq1bw(i,k)
          q2bw(i,k) = q2bw(i,k)+dq2bw(i,k)
          q3bw(i,k) = q3bw(i,k)+dq3bw(i,k)
          q1ybs(i,k) = q1bw(i,k)
          q2ybs(i,k) = q2bw(i,k)
          q3ybs(i,k) = q3bw(i,k)

          ! scalar
          densbw(i,k) = densbw(i,k)+ddensbw(i,k)
          co2bw(i,k) = co2bw(i,k)+dco2bw(i,k)
          densybs(i,k) = densbw(i,k)
          co2ybs(i,k) = co2bw(i,k)
        enddo
      enddo

      do k=kstartr,kendr
        do i=1,n1r
          dsalbw(i,k) = dsalbw(i,k)+ddsalbw(i,k)
          dsalybs(i,k) = dsalbw(i,k)
        enddo
      enddo

      return
      end

!===========================================================
! memory allocations
!***********************************************************
      subroutine outflowBCs_mem_alloc
      use mpih
      use param
      use mpi_param, only: kstart,kend,kstartr,kendr
      use outflow_vars
      implicit none

      allocate( q1be(1:n1,kstart:kend), dq1be(1:n1,kstart:kend) )
      allocate( q2be(1:n1,kstart:kend), dq2be(1:n1,kstart:kend) )
      allocate( q3be(1:n1,kstart:kend), dq3be(1:n1,kstart:kend) )
      allocate( dq1beo(1:n1,kstart:kend) )
      allocate( dq2beo(1:n1,kstart:kend) )
      allocate( dq3beo(1:n1,kstart:kend) )
      allocate( dsalbe(1:n1r,kstartr:kendr), ddsalbe(1:n1r,kstartr:kendr) )
      allocate( ddsalbeo(1:n1r,kstartr:kendr) )
      allocate( densbe(1:n1r,kstartr:kendr), ddensbe(1:n1r,kstartr:kendr) )
      allocate( ddensbeo(1:n1r,kstartr:kendr) )
      allocate( co2be(1:n1r,kstartr:kendr), dco2be(1:n1r,kstartr:kendr) )
      allocate( dco2beo(1:n1r,kstartr:kendr) )
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      q1be(:,:)=0.0d0
      dq1be(:,:)=0.0d0
      dq1beo(:,:)=0.0d0
      q2be(:,:)=0.0d0
      dq2be(:,:)=0.0d0
      dq2beo(:,:)=0.0d0
      q3be(:,:)=0.0d0
      dq3be(:,:)=0.0d0
      dq3beo(:,:)=0.0d0
      dsalbe(:,:)=0.0d0
      ddsalbe(:,:)=0.0d0
      ddsalbeo(:,:)=0.0d0
      densbe(:,:)=0.0d0
      ddensbe(:,:)=0.0d0
      ddensbeo(:,:)=0.0d0
      co2be(:,:)=0.0d0
      dco2be(:,:)=0.0d0
      dco2beo(:,:)=0.0d0

      return
      end

!===========================================================
      subroutine inflowBCs_mem_alloc
      use mpih
      use param
      use mpi_param, only: kstart,kend,kstartr,kendr
      use inflow_vars
      implicit none

      allocate( q1bw(1:n1,kstart:kend), dq1bw(1:n1,kstart:kend)  )
      allocate( q2bw(1:n1,kstart:kend), dq2bw(1:n1,kstart:kend)  )
      allocate( q3bw(1:n1,kstart:kend), dq3bw(1:n1,kstart:kend)  )
      allocate( dsalbw(1:n1r,kstartr:kendr), ddsalbw(1:n1r,kstartr:kendr) )
      allocate( densbw(1:n1r,kstartr:kendr), ddensbw(1:n1r,kstartr:kendr) )
      allocate( co2bw(1:n1r,kstartr:kendr), dco2bw(1:n1r,kstartr:kendr) )
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      q1bw(:,:)=0.0d0
      dq1bw(:,:)=0.0d0
      q2bw(:,:)=0.0d0
      dq2bw(:,:)=0.0d0
      q3bw(:,:)=0.0d0
      dq3bw(:,:)=0.0d0
      dsalbw(:,:)=0.0d0
      ddsalbw(:,:)=0.0d0
      densbw(:,:)=0.0d0
      ddensbw(:,:)=0.0d0
      co2bw(:,:)=0.0d0
      dco2bw(:,:)=0.0d0

      return
      end

!===========================================================
! continuation write BC arrays to files
!***********************************************************
      subroutine inoutflowBCs_write_continua
      use mpih
      use hdf5
      use param
      use mpi_param, only: kstart,kend,kstartr,kendr
      use outflow_vars
      use inflow_vars
      implicit none
      integer :: comm, info
      integer :: ndims
      integer hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: memspace
      integer(HID_T) :: dset_q1v
      integer(HID_T) :: dset_q2v
      integer(HID_T) :: dset_q3v
      integer(HID_T) :: dset_sal,dset_dens,dset_co2
      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 
      integer(HSIZE_T) :: dims(2)
      character(70) namfile

      namfile='continua_inoutflowBCs.h5'

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      ndims = 2
      dims(1)=n1
      dims(2)=n3m


      data_count(1) = n1
      data_count(2) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = kstart-1

      !-- open
      call h5open_f(hdf_error)
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id,hdf_error)

      !-- header
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'dq1beo', H5T_NATIVE_DOUBLE,filespace, dset_q1v, hdf_error)
      call h5dcreate_f(file_id, 'dq2beo', H5T_NATIVE_DOUBLE,filespace, dset_q2v, hdf_error)
      call h5dcreate_f(file_id, 'dq3beo', H5T_NATIVE_DOUBLE,filespace, dset_q3v, hdf_error)
      call h5dcreate_f(file_id, 'ddsalbeo', H5T_NATIVE_DOUBLE,filespace, dset_sal, hdf_error)
      call h5dcreate_f(file_id, 'ddensbeo', H5T_NATIVE_DOUBLE,filespace, dset_dens, hdf_error)
      call h5dcreate_f(file_id, 'dco2beo', H5T_NATIVE_DOUBLE,filespace, dset_co2, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      !-- memspace
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      !-- write
      call h5dget_space_f(dset_q1v, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_q1v, H5T_NATIVE_DOUBLE,dq1beo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)

      call h5dget_space_f(dset_q2v, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_q2v, H5T_NATIVE_DOUBLE,dq2beo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)

      call h5dget_space_f(dset_q3v, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_q3v, H5T_NATIVE_DOUBLE,dq3beo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)

      call h5dget_space_f(dset_sal, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_sal, H5T_NATIVE_DOUBLE,ddsalbeo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)

      call h5dget_space_f(dset_dens, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_dens, H5T_NATIVE_DOUBLE,ddensbeo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)

      call h5dget_space_f(dset_co2, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_co2, H5T_NATIVE_DOUBLE,dco2beo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)

      !-- close header
      call h5dclose_f(dset_q1v, hdf_error)
      call h5dclose_f(dset_q2v, hdf_error)
      call h5dclose_f(dset_q3v, hdf_error)
      call h5dclose_f(dset_sal, hdf_error)
      call h5dclose_f(dset_dens, hdf_error)
      call h5dclose_f(dset_co2, hdf_error)

      !-- close
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
      call h5close_f(hdf_error)

      if(myid.eq.0)write(*,'(5x,a)')'write complete: '//namfile

      return
      end

!===========================================================
! continuation read BC arrays from files
!***********************************************************
      subroutine inoutflowBCs_read_continua
      use mpih
      use hdf5
      use param
      use mpi_param, only: kstart,kend,kstartr,kendr
      use outflow_vars
      use inflow_vars
      implicit none
      integer :: comm, info
      integer :: ndims
      integer hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: memspace
      integer(HID_T) :: dset_q1v
      integer(HID_T) :: dset_q2v
      integer(HID_T) :: dset_q3v
      integer(HID_T) :: dset_sal,dset_dens,dset_co2
      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 
      integer(HSIZE_T) :: dims(2)
      character(70) namfile

      namfile='continua_inoutflowBCs.h5'

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      ndims = 2
      dims(1)=n1
      dims(2)=n3m


      data_count(1) = n1
      data_count(2) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = kstart-1

      !-- open
      call h5open_f(hdf_error)
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fopen_f(namfile, H5F_ACC_RDONLY_F, file_id, hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id,hdf_error)

      !-- header
      call h5dopen_f(file_id, 'dq1beo', dset_q1v, hdf_error)
      call h5dopen_f(file_id, 'dq2beo', dset_q2v, hdf_error)
      call h5dopen_f(file_id, 'dq3beo', dset_q3v, hdf_error)
      call h5dopen_f(file_id, 'ddsalbeo', dset_sal, hdf_error)
      call h5dopen_f(file_id, 'ddensbeo', dset_dens, hdf_error)
      call h5dopen_f(file_id, 'dco2beo', dset_co2, hdf_error)

      !-- memspace
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      !-- write
      call h5dget_space_f(dset_q1v, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dread_f(dset_q1v, H5T_NATIVE_DOUBLE,dq1beo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)

      call h5dget_space_f(dset_q2v, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dread_f(dset_q2v, H5T_NATIVE_DOUBLE,dq2beo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)

      call h5dget_space_f(dset_q3v, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dread_f(dset_q3v, H5T_NATIVE_DOUBLE,dq3beo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)

      call h5dget_space_f(dset_sal, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dread_f(dset_sal, H5T_NATIVE_DOUBLE,ddsalbeo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)

      call h5dget_space_f(dset_dens, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dread_f(dset_dens, H5T_NATIVE_DOUBLE,ddensbeo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)

      call h5dget_space_f(dset_co2, filespace, hdf_error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dread_f(dset_co2, H5T_NATIVE_DOUBLE,dco2beo(1:n1,kstart:kend), dims,hdf_error, file_space_id = filespace,mem_space_id = memspace,xfer_prp = plist_id)


      !-- close header
      call h5dclose_f(dset_q1v, hdf_error)
      call h5dclose_f(dset_q2v, hdf_error)
      call h5dclose_f(dset_q3v, hdf_error)
      call h5dclose_f(dset_sal, hdf_error)
      call h5dclose_f(dset_dens, hdf_error)
      call h5dclose_f(dset_co2, hdf_error)

      !-- close
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
      call h5close_f(hdf_error)

      if(myid.eq.0)write(*,'(5x,a)')'reading complete: '//namfile

      return
      end
