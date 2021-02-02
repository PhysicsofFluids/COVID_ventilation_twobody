      subroutine solxj_scafixval_co2inout(betadx)
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      use outflow_vars, only: dco2be
      use inflow_vars, only: dco2bw
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amjl(0:m2m),apjl(0:m2m),acjl(0:m2)
      real :: fjl(0:m2,m1m)

      acjl(0)     = 1.0d0
      acjl(1:n2m) = (1.d0+2.d0*betadx)
      acjl(n2)    = 1.0d0

      amjl(0:n2m-1) = -betadx
      amjl(n2m)     = 0.0d0

      apjl(0)     = 0.0d0
      apjl(1:n2m) = -betadx

      call ddttrfb(n2+1,amjl(0),acjl(0),apjl(0),info)

      do kc=kstart,kend
        do ic=1,n1m
          fjl(0,ic)=dco2bw(ic,kc)
          do jc=1,n2m
            fjl(jc,ic)=rhs(ic,jc,kc)
          enddo
          fjl(n2,ic)=dco2be(ic,kc)
        end do

        call ddttrsb('N',n2+1,n1m,amjl(0),acjl(0),apjl(0),fjl(0,1),n2+1,info)

        do ic=1,n1m
          do jc=1,n2m
            rhs(ic,jc,kc) = fjl(jc,ic)
          enddo
        end do
      end do 

      return
      end
