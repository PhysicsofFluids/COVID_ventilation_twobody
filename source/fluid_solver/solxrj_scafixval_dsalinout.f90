      subroutine solxrj_scafixval_dsalinout(betadx)
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstartr,kendr
      use outflow_vars, only: ddsalbe
      use inflow_vars, only: ddsalbw
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amjl(0:m2mr),apjl(0:m2mr),acjl(0:m2r)
      real :: fjl(0:m2r,m1mr)

      acjl(0)     = 1.0d0
      acjl(1:n2mr) = (1.d0+2.d0*betadx)
      acjl(n2r)    = 1.0d0

      amjl(0:n2mr-1) = -betadx
      amjl(n2mr)     = 0.0d0

      apjl(0)     = 0.0d0
      apjl(1:n2mr) = -betadx

      call ddttrfb(n2r+1,amjl(0),acjl(0),apjl(0),info)

      do kc=kstartr,kendr
        do ic=1,n1mr
          fjl(0,ic)=ddsalbw(ic,kc)
          do jc=1,n2mr
            fjl(jc,ic)=rhs(ic,jc,kc)
          enddo
          fjl(n2r,ic)=ddsalbe(ic,kc)
        end do

        call ddttrsb('N',n2r+1,n1mr,amjl(0),acjl(0),apjl(0),fjl(0,1),n2r+1,info)

        do ic=1,n1mr
          do jc=1,n2mr
            rhs(ic,jc,kc) = fjl(jc,ic)
          enddo
        end do
      end do 

      return
      end
