!  this subroutine perform the calculation of dph , periodic direction
!  along x3 and x1to use the real fourier transform
      subroutine phcalc
      use param
      use local_arrays, only: dph
      use mpi_param
      use mpih
      implicit none
!      include "mkl_service.fi"
      integer :: i,j,k,info
      real :: coefnorm
      real :: xr(m2m,m1m)
      real :: xa(m2m,m1m)

      integer :: n2mh,jmh
      real,allocatable,dimension(:,:,:) :: dpht,dpho

      real, dimension(n3m) :: b,r,u
      real, dimension(n3m,3) :: opdz
      logical losing

      allocate(dpht(1:n3m,1:n1m,jstart:jend))
      allocate(dpho(1:n1m,1:n2m,kstart:kend))

      coefnorm = 1.d0/(4*dble(n1m-1)*dble(n2m-1))

      do k=kstart,kend
        do j=1,n2m
          do i=1,n1m
           xr(j,i)=dph(i,j,k)
          enddo
        enddo
        
        call dfftw_execute_r2r(fwd_plan,xr,xa)

        do j=1,n2m
         do i=1,n1m
         dpho(i,j,k)=xa(j,i)*coefnorm
        enddo
        enddo
      end do

      !=====Transpose
      call PackZ_UnpackRP(dpho,dpht)

      !=====Tridiagonal matrix solver
      opdz(:,:) = 0.0d0
      opdz(2:n3m,1) = amphk(2:n3m)
      opdz(1:n3m,2) = acphk(1:n3m)
      opdz(1:n3m-1,3) = apphk(1:n3m-1)

      losing=.false.

      do j=jstart,jend
        do i=1,n1m
          do k = 1,n3m
            b(k) = opdz(k,2)-ak2(j)-ak1(i)
            r(k) = dpht(k,i,j)
          enddo
          if(i.eq.1 .and. j.eq.1) losing=.true.
          call tridag(opdz(:,1),b,opdz(:,3),r,u,n3m,losing)
          losing=.false.
          dpht(:,i,j) = 0.0d0
          do k = 1,n3m
            dpht(k,i,j) = u(k)
          enddo
        enddo
      enddo


      !=====Transpose
      call PackR_UnpackZP(dpht,dpho)

      do k=kstart,kend
       do j=1,n2m
        do i=1,n1m
          xa(j,i)=dpho(i,j,k)
        enddo
       end do

      call dfftw_execute_r2r(bck_plan,xa,xr)

       do j=1,n2m
         do i=1,n1m
           dph(i,j,k)=xr(j,i)
         enddo
       end do
      end do

      if(allocated(dpht)) deallocate(dpht)
      if(allocated(dpho)) deallocate(dpho)

      return
      end

!========================================================
 subroutine tridag(a,b,c,r,u,ns,losing)
!------------------------------------------------
!     Solve tridiagonal matrix (based on Numerical Recipes)
!     a.u(i-1)+bu(i)+cu(i+1)=r
!------------------------------------------------
      use param
      use local_arrays, only: dph
      use mpi_param
      use mpih
      implicit none
!      include "mkl_service.fi"

      !--- local variables
      integer, intent(in) :: ns
      real, dimension(ns), intent(in ) :: a,b,c,r
      real, dimension(ns), intent(out) :: u
      integer, parameter  :: nsmax=4096
      integer  :: ntr
      real :: bet
      real, dimension(nsmax) :: gamtmp
      logical  :: losing


      if(b(1) == 0.0d0) then
         write(*,*) 'tridag: rewrite equations'
      endif
      bet=b(1)
      u(1)=r(1)/bet
    
      do ntr=2,ns-1
           gamtmp(ntr)=c(ntr-1)/bet
           bet=b(ntr)-a(ntr)*gamtmp(ntr)
           u(ntr)=(r(ntr)-a(ntr)*u(ntr-1))/bet
      enddo

      !--- Check for singular matrix
      ntr=ns
      gamtmp(ntr)=c(ntr-1)/bet
      bet=b(ntr)-a(ntr)*gamtmp(ntr)
      if(losing) then
         u(ntr) = 0.0d0
      else
         u(ntr)=(r(ntr)-a(ntr)*u(ntr-1))/bet
      endif


      do ntr=ns-1,1,-1
         u(ntr)=u(ntr)-gamtmp(ntr+1)*u(ntr+1)
      enddo

      return
end subroutine tridag
