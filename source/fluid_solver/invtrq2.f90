      subroutine invtrq2
      use param
      use local_arrays, only: q2,pr,rhs,dph,ru2
      use mpi_param, only: kstart,kend
      use walls_vars

      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      real    :: udx2,amm,app,acc
      real    :: dcq2,dpx22
      real    :: d22q2,d33q2,d11q2
      real    :: alre,udx1q,udx2q
      real    :: a33,a33p,a33m

      alre=al*nu

      udx2=dx2*al
      udx1q=dx1q
      udx2q=dx2q
      do kc=kstart,kend
          km=kmv(kc)
          kp=kpv(kc)
          amm=am3sk(kc)
          acc=ac3sk(kc)
          app=ap3sk(kc)
          do jc=2,n2m
           jm=jmv(jc)
           jp=jpv(jc)
            do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)

            if(ic.eq.1) then
               d11q2=(q2(ip,jc,kc)-3.d0*q2(ic,jc,kc)+0.0d0)*udx1q
               !d11q2=(q2(ip,jc,kc)-1.d0*q2(ic,jc,kc)+0.0d0)*udx1q
            elseif(ic.eq.n1m) then
               d11q2=(0.0d0-3.d0*q2(ic,jc,kc)+q2(im,jc,kc))*udx1q
               !d11q2=(0.0d0-1.d0*q2(ic,jc,kc)+q2(im,jc,kc))*udx1q
            else
               d11q2=(q2(ip,jc,kc)-2.d0*q2(ic,jc,kc)+q2(im,jc,kc))*udx1q
            endif

            if(jc.eq.2) then
               d22q2=(q2(ic,jp,kc)-2.d0*q2(ic,jc,kc)+q2ybs(ic,kc))*udx2q
            elseif(jc.eq.n2m) then
               d22q2=(q2ybn(ic,kc)-2.d0*q2(ic,jc,kc)+q2(ic,jm,kc))*udx2q
            else
               d22q2=(q2(ic,jp,kc)-2.d0*q2(ic,jc,kc)+q2(ic,jm,kc))*udx2q
            endif

            if(kc.eq.1) then
               d33q2=q2(ic,jc,kp)*app +q2(ic,jc,kc)*acc+0.0d0*amm
            elseif(kc.eq.n3m) then
               d33q2=0.0d0*app +q2(ic,jc,kc)*acc+q2(ic,jc,km)*amm
            else
               d33q2=q2(ic,jc,kp)*app +q2(ic,jc,kc)*acc+q2(ic,jc,km)*amm
            endif

            dcq2=d22q2+d33q2+d11q2

            dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2

            rhs(ic,jc,kc)=(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc) +alre*dcq2 -dpx22)*dt
            ru2(ic,jc,kc)=dph(ic,jc,kc)      
         enddo
       enddo
      enddo

      !-- implicit matrix solver
      !call solxi_fs(beta*al*dx1q)
      call solxi_ns(beta*al*dx1q)
      call solxj_endwall(beta*al*dx2q)   
      call solq2k

      return
      end
