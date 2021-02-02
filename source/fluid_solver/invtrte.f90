      subroutine invtrte
      use param
      use local_arrays, only: dens,hro,ruro,rhs
      use mpi_param, only: kstart,kend
      use walls_vars
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real    :: dq32,dq33,dcq3,dq31
      real    :: app,acc,amm
      real    :: alpecl,udx1q,udx2q
      real    :: del1, del2, fcder
      !!real    :: dsaltop_outflow

      alpecl=al *Dw
      udx1q=dx1q
      udx2q=dx2q

      do kc=kstart,kend

      !-- inner
      if( (kc.ge.2) .and. (kc.le.n3m-1) ) then
        km=kmv(kc)
        kp=kpv(kc)
        app=ap3ssk(kc)
        acc=ac3ssk(kc)
        amm=am3ssk(kc)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
            dq31=(dens(ip,jc,kc)-2.d0*dens(ic,jc,kc) +dens(im,jc,kc))*udx1q
            if(jc.eq.1) then
               dq32=(dens(ic,jp,kc)-2.d0*dens(ic,jc,kc)+densybs(ic,kc))*udx2q
            elseif(jc.eq.n2m) then
               dq32=(densybn(ic,kc)-2.d0*dens(ic,jc,kc)+dens(ic,jm,kc))*udx2q
            else
               dq32=(dens(ic,jp,kc)-2.d0*dens(ic,jc,kc)+dens(ic,jm,kc))*udx2q
            endif
            dq33= dens(ic,jc,kp)*app+dens(ic,jc,kc)*acc+dens(ic,jc,km)*amm
            dcq3=dq32+dq33+dq31

            rhs(ic,jc,kc)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc)+alpecl*dcq3)*dts
            ruro(ic,jc,kc)=hro(ic,jc,kc)
            enddo
        enddo
      endif

      !-- bottom
      if(kc.eq.1) then
        app=ap3ssk(kc)
        acc=ac3ssk(kc)
        amm=am3ssk(kc)
        kp = kc + 1
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
            dq31=(dens(ip,jc,kc)-2.d0*dens(ic,jc,kc)+dens(im,jc,kc))*udx1q
            dq32=(dens(ic,jp,kc)-2.d0*dens(ic,jc,kc)+dens(ic,jm,kc))*udx2q
            dq33=dens(ic,jc,kp)*app+dens(ic,jc,kc)*acc
            dcq3=dq32+dq33+dq31
            rhs(ic,jc,kc)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc) +alpecl*dcq3)*dts
            ruro(ic,jc,kc)=hro(ic,jc,kc)
          enddo
        enddo
      endif

      !-- top
      if(kc.eq.n3m) then
        app=ap3ssk(kc)
        acc=ac3ssk(kc)
        amm=am3ssk(kc)
        km = kc - 1
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
            dq31=( dens(ip,jc,kc)-2.d0*dens(ic,jc,kc)+dens(im,jc,kc))*udx1q
            dq32=(dens(ic,jp,kc)-2.d0*dens(ic,jc,kc)+dens(ic,jm,kc))*udx2q
            !dq33=dsaltop*app+dsal(ic,jc,kc)*acc+dsal(ic,jc,km)*amm
            dq33=dens(ic,jc,kc)*acc+dens(ic,jc,km)*amm
            dcq3=dq32+dq33+dq31

            rhs(ic,jc,kc)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc)+alpecl*dcq3)*dts
            ruro(ic,jc,kc)=hro(ic,jc,kc)
          enddo
        enddo
      endif
      enddo

      !-- implicit matrix solver
      call solxi( alpecl*dts*0.5d0*dx1q )
      call solxj_scafixval_densinout( alpecl*dts*0.5d0*dx2q )
      call soltek

      return
      end
