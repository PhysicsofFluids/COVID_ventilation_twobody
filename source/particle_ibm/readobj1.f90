!-------------------------------------------------------
!     reading in obj file and storing geometry data
!-------------------------------------------------------
      subroutine readobj1
      USE param
      USE mpih
      USE local_arrays, only: isbody1
      USE mpi_param, only: kstartr,kendr
      IMPLICIT NONE
      character*50  :: input_file_name
      integer :: jc,kc,ic
      logical :: inside
      real :: p(3),chksum,buf,val_max,val_min,val_max_all,val_min_all
      real :: xmax,ymax,zmax, xmin,ymin,zmin, xdiff,ydiff,zdiff, alldiff

      integer ( kind = 4 ), allocatable, dimension(:,:) :: face_node
      integer ( kind = 4 ) face_num
      integer ( kind = 4 ), allocatable, dimension(:) :: face_order
      integer ( kind = 4 ) ierror
      integer ( kind = 4 ) node_num
      real ( kind = 8 ), allocatable, dimension(:,:) :: node_xyz
      real ( kind = 8 ), allocatable, dimension(:,:) :: normal_vector
      integer ( kind = 4 ) normal_num
      integer ( kind = 4 ) order_max
      integer ( kind = 4 ), allocatable, dimension(:,:) :: vertex_normal

      input_file_name = objfx1

      call obj_size ( input_file_name, node_num, face_num, normal_num, &
       order_max )
   
      allocate ( face_node(order_max,face_num) )
      allocate ( face_order(face_num) )
      allocate ( node_xyz(3,node_num) )
      allocate ( normal_vector(3,normal_num) )
      allocate ( vertex_normal(order_max,face_num) )
   
      call obj_read ( input_file_name, node_num, face_num, normal_num, &
       order_max, node_xyz, face_order, face_node, normal_vector, vertex_normal )
   
      !-- rescale obj
      xmax=maxval(node_xyz(1,:))
      xmin=minval(node_xyz(1,:))
      ymax=maxval(node_xyz(2,:))
      ymin=minval(node_xyz(2,:))
      zmax=maxval(node_xyz(3,:))
      zmin=minval(node_xyz(3,:))
      xdiff=xmax-xmin
      ydiff=ymax-ymin
      zdiff=zmax-zmin
      if(xdiff.ge.ydiff .and. xdiff.ge.zdiff) then 
        alldiff=xdiff
      elseif(ydiff.ge.xdiff .and. ydiff.ge.zdiff) then
        alldiff=ydiff
      else
        alldiff=zdiff
      endif
      
      node_xyz(1,:)=(node_xyz(1,:)-(xmax+xmin)/2.d0)/alldiff
      node_xyz(2,:)=(node_xyz(2,:)-(ymax+ymin)/2.d0)/alldiff
      node_xyz(3,:)=(node_xyz(3,:)-(zmax+zmin)/2.d0)/alldiff

      node_xyz(1,:)=node_xyz(1,:)*sclf1+rext1*xpos1
      node_xyz(2,:)=node_xyz(2,:)*sclf1+rext2/2.d0
      node_xyz(3,:)=node_xyz(3,:)*sclf1+sclf1/2.d0 !alx3/4.d0
      isbody1(:,:,:) = 0.d0
      do kc=kstartr,kendr
        do jc=1,n2r
          do ic=1,n1r
            p = (/xmr(ic),ymr(jc),zmr(kc)/)
            call polyhedron_contains_point_3d ( node_num, face_num, order_max, node_xyz, face_order, face_node, p, inside )
            if(inside) then
              isbody1(ic,jc,kc) = 1.d0
            else
              isbody1(ic,jc,kc) = 0.d0
            endif
          enddo
        enddo
      enddo

      chksum=sum(isbody1(:,:,:))
      call MPI_REDUCE(chksum,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) write(*,*)'  Sum isbody1: ',buf

      deallocate ( face_node )
      deallocate ( face_order )
      deallocate ( node_xyz )
      deallocate ( normal_vector )
      deallocate ( vertex_normal )


      return
      end subroutine readobj1
