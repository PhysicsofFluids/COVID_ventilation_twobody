      subroutine invtrsa
      use param
      use local_arrays, only:q2,dsal,hsal,rusal,rhsr
      use mpi_param, only: kstartr,kendr
      use walls_vars
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real    :: dq32,dq33,dcq3,dq31
      real    :: app,acc,amm
      real    :: alpecl,udx1qr,udx2qr
      real    :: del1, del2, fcder
      real    :: dsaltop_outflow

      alpecl=al *Dw
      udx1qr=dx1qr
      udx2qr=dx2qr

      do kc=kstartr,kendr

      !-- inner
      if( (kc.ge.2) .and. (kc.le.n3mr-1) ) then
        km=kmvr(kc)
        kp=kpvr(kc)
        app=ap3sskr(kc)
        acc=ac3sskr(kc)
        amm=am3sskr(kc)
        do jc=1,n2mr
          jm=jmvr(jc)
          jp=jpvr(jc)
          do ic=1,n1mr
            im=imvr(ic)
            ip=ipvr(ic)
            dq31=(dsal(ip,jc,kc)-2.d0*dsal(ic,jc,kc) +dsal(im,jc,kc))*udx1qr
            if(jc.eq.1) then
               dq32=(dsal(ic,jp,kc)-2.d0*dsal(ic,jc,kc)+dsalybs(ic,kc))*udx2qr
            elseif(jc.eq.n2m) then
               dq32=(dsalybn(ic,kc)-2.d0*dsal(ic,jc,kc)+dsal(ic,jm,kc))*udx2qr
            else
               dq32=(dsal(ic,jp,kc)-2.d0*dsal(ic,jc,kc)+dsal(ic,jm,kc))*udx2qr
            endif
            dq32=(dsal(ic,jp,kc)-2.d0*dsal(ic,jc,kc)+dsal(ic,jm,kc))*udx2qr
            dq33= dsal(ic,jc,kp)*app+dsal(ic,jc,kc)*acc+dsal(ic,jc,km)*amm
            dcq3=dq32+dq33+dq31

            rhsr(ic,jc,kc)=(ga*hsal(ic,jc,kc)+ro*rusal(ic,jc,kc)+alpecl*dcq3)*dts
            rusal(ic,jc,kc)=hsal(ic,jc,kc)
            enddo
        enddo
      endif

      !-- bottom
      if(kc.eq.1) then
        app=ap3sskr(kc)
        acc=ac3sskr(kc)
        amm=am3sskr(kc)
        kp = kc + 1
        do jc=1,n2mr
          jm=jmvr(jc)
          jp=jpvr(jc)
          do ic=1,n1mr
            im=imvr(ic)
            ip=ipvr(ic)
            dq31=(dsal(ip,jc,kc)-2.d0*dsal(ic,jc,kc)+dsal(im,jc,kc))*udx1qr
            dq32=(dsal(ic,jp,kc)-2.d0*dsal(ic,jc,kc)+dsal(ic,jm,kc))*udx2qr
            dq33=dsal(ic,jc,kp)*app+dsal(ic,jc,kc)*acc
            !dq33=dsal(ic,jc,kp)*app+dsal(ic,jc,kc)*acc+dsalbe(ic,jc)*amm
            dcq3=dq32+dq33+dq31
            rhsr(ic,jc,kc)=(ga*hsal(ic,jc,kc)+ro*rusal(ic,jc,kc) +alpecl*dcq3)*dts
            rusal(ic,jc,kc)=hsal(ic,jc,kc)
          enddo
        enddo
      endif

      !-- top
      if(kc.eq.n3mr) then
        app=ap3sskr(kc)
        acc=ac3sskr(kc)
        amm=am3sskr(kc)
        km = kc - 1
        do jc=1,n2mr
          jm=jmvr(jc)
          jp=jpvr(jc)
          do ic=1,n1mr
            im=imvr(ic)
            ip=ipvr(ic)
            dq31=( dsal(ip,jc,kc)-2.d0*dsal(ic,jc,kc)+dsal(im,jc,kc))*udx1qr
            dq32=(dsal(ic,jp,kc)-2.d0*dsal(ic,jc,kc)+dsal(ic,jm,kc))*udx2qr
            !dq33=dsaltop*app+dsal(ic,jc,kc)*acc+dsal(ic,jc,km)*amm
            dq33=dsal(ic,jc,kc)*acc+dsal(ic,jc,km)*amm
            !dq33=dsalbw(ic,jc)*app+dsal(ic,jc,kc)*acc+dsal(ic,jc,km)*amm
            dcq3=dq32+dq33+dq31

            rhsr(ic,jc,kc)=(ga*hsal(ic,jc,kc)+ro*rusal(ic,jc,kc)+alpecl*dcq3)*dts
            rusal(ic,jc,kc)=hsal(ic,jc,kc)
          enddo
        enddo
      endif
      enddo

      !-- implicit matrix solver
      call solxri( alpecl*dts*0.5d0*dx1qr )
      call solxrj_scafixval_dsalinout( alpecl*dts*0.5d0*dx2qr )
      call solsak

      return
      end
