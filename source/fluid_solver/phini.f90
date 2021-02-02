      subroutine fftqua
      use param
      integer :: n2mh,n2mp,j,i,n1mh,n1mp

      !  Eigenvalue
      do i=1,n1m
        ak1(i)=(2.d0*sin(pi*(i-1)/2/(n1m-1))*(dble(n1m)/rext1))**2
      enddo

      do j=1,n2m
        ak2(j)=(2.d0*sin(pi*(j-1)/2/(n2m-1))*(dble(n2m)/rext2))**2
      enddo

      return
      end

!=====================================================
      subroutine phini
      use param
      use mpi_param
      use mpih
      implicit none
      integer,parameter ::  fftw_es =64

      real, dimension(m2m,m1m) :: xr
      real, dimension(m2m,m1m) :: xa
    
      !RO   Initialize tridiag matrices
      call tridiag_matrices   

      !m    Initialize FFTW
      call fftqua

      call dfftw_plan_r2r_2d(fwd_plan,m2m,m1m,xr,xa,3,3,fftw_es)
      call dfftw_plan_r2r_2d(bck_plan,m2m,m1m,xa,xr,3,3,fftw_es)

      return
      end
      
!=======================================================================
      subroutine tridiag_matrices
      use param
      implicit none
      integer  :: kc,km,kp
      real :: ugmmm,a33icc,a33icp

      !   tridiagonal matrix coefficients at each k and i
      !   x1 and x3 cartesian coordinates
      do kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        a33icc=kmc(kc)*dx3q/g3rc(kc)
        a33icp=kpc(kc)*dx3q/g3rc(kp)
        ugmmm=1.0d0/g3rm(kc)
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        if(kc.eq.1)then
        acphk(kc)=-(apphk(kc))
        elseif(kc.eq.n3m)then
        acphk(kc)=-(amphk(kc))
        else
        acphk(kc)=-(amphk(kc)+apphk(kc))
        endif
      enddo

      end subroutine tridiag_matrices

!==================================================
      subroutine phend
      use param
      implicit none

      call dfftw_destroy_plan(fwd_plan)
      call dfftw_destroy_plan(bck_plan)

      call dfftw_cleanup_threads()

      return
      end subroutine phend
