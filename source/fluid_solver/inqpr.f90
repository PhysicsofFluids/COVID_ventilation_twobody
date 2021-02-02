!***********************************************************************
!                       INITIAL CONDITION                              *
!***********************************************************************
      subroutine inqpr
      use param
      use local_arrays, only: q2,q3,dens,q1,dsal,co2
      use mpi_param, only: kstart,kend,kstartr,kendr
      use mpih
      implicit none
      integer :: j,k,i
      real(8) :: eps, varptb

      ! prescribed initial distributions of T and S
      eps = 1.d-2
      call random_seed()

      do k=kstart,kend
       do j=1,n2m
        do i=1,n1m
          dens(i,j,k) = 0.0d0
          co2(i,j,k) = 0.0d0
        enddo
       enddo
      end do


      do k=kstartr,kendr
       do j=1,n2mr
        do i=1,n1mr
         !call random_number(varptb)
         dsal(i,j,k) = 0.d0 !dsalbot+(dsaltop-dsalbot)*zmr(k)/alx3+eps*(2.d0*varptb-1.d0)!+eps*dsin((i-1)*1.0d0/(n1mr-1)*3.141596)!+eps*(2.d0*varptb-1.d0)!+eps*dsin((i-1)*1.0d0/(n1mr-1)*3.141596)
        enddo
       enddo
      end do

      !  velocity field
      do k=kstart-1,kend+1
        do j=1,n2
          do i=1,n1
            q1(i,j,k)=0.d0
            q2(i,j,k)=0.d0
            q3(i,j,k)=0.d0
          enddo
        enddo
      enddo

      return
      end
