      subroutine solxj_ns_q1out(betadx)
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      use outflow_vars, only: dq1be
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amjl(m2m),apjl(m2m),acjl(m2)
      real :: fjl(m2,m1m)

      amjl(1:n2m-1)=-betadx
      amjl(n2m)=0.d0
      apjl(1:n2m)=-betadx

      acjl(1) = (1.d0+3.d0*betadx)
      acjl(2:n2m) = (1.d0+2.d0*betadx)
      acjl(n2) = (1.d0)

      call ddttrfb(n2,amjl,acjl,apjl,info)

      do kc=kstart,kend
        do ic=1,n1m
          do jc=1,n2m
            fjl(jc,ic)=rhs(ic,jc,kc)
          enddo
          fjl(n2,ic)=dq1be(ic,kc)
        end do

        call ddttrsb('N',n2,n1m,amjl,acjl,apjl,fjl,n2,info)

        do ic=1,n1m
          do jc=1,n2m
            rhs(ic,jc,kc) = fjl(jc,ic)
          enddo
        end do
      end do 

      return
      end
