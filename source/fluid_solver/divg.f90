!
!     this subroutine calculates divg(q).
!     q are the "fluxes"
!
      subroutine divg
      use param
      use local_arrays, only: q1,q2,q3,dph
      use mpi_param, only: kstart,kend
      use walls_vars
      implicit none
      integer :: jc,jp,kc,kp,ic,ip
      real    :: usdtal,dqcap
      real    :: q1up,q1low,q2up,q2low,q3up,q3low

      usdtal = 1.d0/(dt*al)

      do kc=kstart,kend

        q1up=1.0;q1low=1.0;q2up=1.0;q2low=1.0;q3up=1.0;q3low=1.0
        kp=kc+1
        do jc=1,n2m
          jp=jpv(jc)
          do ic=1,n1m
            ip=ipv(ic)

            if(ic.eq.n1m) then
               q1up=0.0d0
            else
               q1up=1.0d0
            endif
            if(ic.eq.1) then
               q1low=0.0d0
            else
               q1low=1.0d0
            endif

            if(kc.eq.n3m) then
               q3up=0.0d0
            else
               q3up=1.0d0
            endif

            if(kc.eq.1) then
               q3low=0.0d0
            else
               q3low=1.0d0
            endif

            if(jc.eq.1) then
               dqcap= (q1(ip,jc,kc)*q1up-q1(ic,jc,kc)*q1low)*dx1&
        &            +(q2(ic,jp,kc)-q2ybs(ic,kc))*dx2&
        &            +(q3(ic,jc,kp)*q3up-q3(ic,jc,kc)*q3low)*udx3m(kc)
            elseif(jc.eq.n2m) then
               dqcap= (q1(ip,jc,kc)*q1up-q1(ic,jc,kc)*q1low)*dx1&
        &            +(q2ybn(ic,kc)-q2(ic,jc,kc))*dx2&
        &            +(q3(ic,jc,kp)*q3up-q3(ic,jc,kc)*q3low)*udx3m(kc)
            else
               dqcap= (q1(ip,jc,kc)*q1up-q1(ic,jc,kc)*q1low)*dx1&
        &            +(q2(ic,jp,kc)-q2(ic,jc,kc))*dx2&
        &            +(q3(ic,jc,kp)*q3up-q3(ic,jc,kc)*q3low)*udx3m(kc)
            endif

            dph(ic,jc,kc)=dqcap*usdtal
          enddo
        enddo
      enddo

      return
      end
