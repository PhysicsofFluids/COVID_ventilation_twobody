      subroutine hdnlq3
      use param
      use local_arrays, only: q1,q2,q3,qcap,dens,co2
      use local_arrays, only:dsal
      use mgrd_arrays, only:dsalc
      use mpi_param, only: kstart,kend
      use walls_vars
      implicit none
      integer :: ic,jc,kc,js
      integer :: km,kp,jm,jp,im,ip
      real    :: h32,h33,h31
      real    :: densit,dsalit,co2it,udx1,udx2
      integer :: kstartp

      if(kstart.eq.1) then
        kstartp=2
      else
        kstartp=kstart
      endif

      udx1=dx1*0.25d0
      udx2=dx2*0.25d0

      do kc=kstartp,kend
        km=kc-1
        kp=kc+1
        do jc=1,n2m
          jm=jc-1
          jp=jc+1
          do ic=1,n1m
            im=ic-1
            ip=ic+1
            if(ic.eq.1) then
               h31=( (q1(ip,jc,kc)+q1(ip,jc,km))*(q3(ip,jc,kc)+q3(ic,jc,kc)))*udx1
            elseif(ic.eq.n1m) then
               h31=(-(q1(ic,jc,kc)+q1(ic,jc,km))*(q3(ic,jc,kc)+q3(im,jc,kc)))*udx1
            else
               h31=( (q1(ip,jc,kc)+q1(ip,jc,km))*(q3(ip,jc,kc)+q3(ic,jc,kc))&
        &           -(q1(ic,jc,kc)+q1(ic,jc,km))*(q3(ic,jc,kc)+q3(im,jc,kc)))*udx1
            endif

            if(jc.eq.1) then
               h32=( (q2(ic,jp,kc)+q2(ic,jp,km))*(q3(ic,jp,kc)+q3(ic,jc,kc))&
        &           -(q2ybs(ic,kc)+q2ybs(ic,km))*(q3(ic,jc,kc)+q3ybs(ic,kc)))*udx2
            elseif(jc.eq.n2m) then
               h32=( (q2ybn(ic,kc)+q2ybn(ic,km))*(q3ybn(ic,kc)+q3(ic,jc,kc))&
        &           -(q2(ic,jc,kc)+q2(ic,jc,km))*(q3(ic,jc,kc)+q3(ic,jm,kc)))*udx2
            else
               h32=( (q2(ic,jp,kc)+q2(ic,jp,km))*(q3(ic,jp,kc)+q3(ic,jc,kc))&
        &           -(q2(ic,jc,kc)+q2(ic,jc,km))*(q3(ic,jc,kc)+q3(ic,jm,kc)))*udx2
            endif

            if(kc.eq.2) then
               h33=( (q3(ic,jc,kp)+q3(ic,jc,kc))*(q3(ic,jc,kp)+q3(ic,jc,kc))&
        &           -(q3(ic,jc,kc)+       0.0d0)*(q3(ic,jc,kc)+      0.0d0))*udx3c(kc)*0.25d0
            elseif(kc.eq.n3m) then
               h33=( (0.0d0       +q3(ic,jc,kc))*(0.0d0       +q3(ic,jc,kc))&
        &           -(q3(ic,jc,kc)+q3(ic,jc,km))*(q3(ic,jc,kc)+q3(ic,jc,km)))*udx3c(kc)*0.25d0
            else
               h33=( (q3(ic,jc,kp)+q3(ic,jc,kc))*(q3(ic,jc,kp)+q3(ic,jc,kc))&
        &           -(q3(ic,jc,kc)+q3(ic,jc,km))*(q3(ic,jc,kc)+q3(ic,jc,km)))*udx3c(kc)*0.25d0
            endif

            densit = (dens(ic,jc,kc)+dens(ic,jc,kc-1))*0.5d0
            dsalit = (dsalc(ic,jc,kc)+dsalc(ic,jc,kc-1))*0.5d0
            co2it  = (co2(ic,jc,kc)+co2(ic,jc,kc-1))*0.5d0

            qcap(ic,jc,kc) = -(h31+h32+h33)+densit*byct+dsalit*bycs+co2it*bycco2
          enddo
        enddo
      enddo
      return
      end
