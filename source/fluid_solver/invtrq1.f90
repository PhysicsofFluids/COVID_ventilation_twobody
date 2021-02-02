      subroutine invtrq1
      use param
      use local_arrays, only: pr,rhs,ru1,q1,dq
      use mpi_param, only: kstart,kend
      use walls_vars

      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      real    :: udx1,amm,acc,app
      real    :: dcq1,dpx11
      real    :: d22q1,d33q1,d11q1
      real    :: alre,udx1q,udx2q

      alre=al*nu

      udx1=dx1*al
      udx1q=dx1q
      udx2q=dx2q

      do kc=kstart,kend
        km=kmv(kc)
        kp=kpv(kc)
        amm=am3sk(kc)
        acc=ac3sk(kc)
        app=ap3sk(kc)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=2,n1m
            im=imv(ic)
            ip=ipv(ic)

            if(ic.eq.2) then
               d11q1=(q1(ip,jc,kc)-3.d0*q1(ic,jc,kc)+0.0d0)*udx1q
            elseif(ic.eq.n1m) then
               d11q1=(0.0d0-3.d0*q1(ic,jc,kc)+q1(im,jc,kc))*udx1q
            else
               d11q1=(q1(ip,jc,kc)-2.d0*q1(ic,jc,kc)+q1(im,jc,kc))*udx1q
            endif
            
            if(jc.eq.1) then
               d22q1=(q1(ic,jp,kc) -2.d0*q1(ic,jc,kc)+q1ybs(ic,kc))*udx2q
            elseif(jc.eq.n2m) then
               d22q1=(q1ybn(ic,kc) -2.d0*q1(ic,jc,kc)+q1(ic,jm,kc))*udx2q
            else
               d22q1=(q1(ic,jp,kc) -2.d0*q1(ic,jc,kc)+q1(ic,jm,kc))*udx2q
            endif
            
            if(kc.eq.1) then
               d33q1=q1(ic,jc,kp)*app +q1(ic,jc,kc)*acc+0.0d0*amm
            elseif(kc.eq.n3m) then
               d33q1=0.0d0*app +q1(ic,jc,kc)*acc+q1(ic,jc,km)*amm
            else
               d33q1=q1(ic,jc,kp)*app +q1(ic,jc,kc)*acc+q1(ic,jc,km)*amm
            endif

            dcq1=d11q1+d22q1+d33q1
            
            dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*udx1

            rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc)+alre*dcq1-dpx11)*dt
            ru1(ic,jc,kc)=dq(ic,jc,kc)
          enddo
        enddo
      enddo

      !-- implicit matrix solver
      call solxi_endwall(beta*al*dx1q)
      call solxj_ns_q1out(beta*al*dx2q)
      call solq1k

      return
      end
