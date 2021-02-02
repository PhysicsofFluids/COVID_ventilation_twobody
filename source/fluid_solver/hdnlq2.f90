      subroutine hdnlq2
      use param
      use local_arrays, only: q1,q2,q3,dph
      use mgrd_arrays, only: dsalc
      use mpi_param, only: kstart,kend
      use walls_vars
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip
      integer :: kmm,kpp
      real    :: h22,h23,udx1,udx2,h21
      real    :: dsalit

      udx1=dx1*0.25d0
      udx2=dx2*0.25d0

    do kc=kstart,kend
        kmm=kc-1
        kpp=kc+1
        kp=kc+1
        do jc=2,n2m
          jm=jc-1
          jp=jc+1
          do ic=1,n1m
            im=ic-1
            ip=ic+1
            if(ic.eq.1) then
               h21=( (q2(ip,jc,kc)+q2(ic,jc,kc))*(q1(ip,jc,kc)+q1(ip,jm,kc)))*udx1
            elseif(ic.eq.n1m) then
               h21=(-(q2(ic,jc,kc)+q2(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jm,kc)))*udx1
            else
               h21=( (q2(ip,jc,kc)+q2(ic,jc,kc))*(q1(ip,jc,kc)+q1(ip,jm,kc))&
        &           -(q2(ic,jc,kc)+q2(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jm,kc)))*udx1
            endif

            if(jc.eq.2) then
               h22=( (q2(ic,jp,kc)+q2(ic,jc,kc))*(q2(ic,jp,kc)+q2(ic,jc,kc))&
        &           -(q2ybs(ic,kc)+q2(ic,jc,kc))*(q2ybs(ic,kc)+q2(ic,jc,kc)))*udx2
            elseif(jc.eq.n2m) then
               h22=( (q2ybn(ic,kc)+q2(ic,jc,kc))*(q2ybn(ic,kc)+q2(ic,jc,kc))&
        &           -(q2(ic,jm,kc)+q2(ic,jc,kc))*(q2(ic,jm,kc)+q2(ic,jc,kc)))*udx2
            else
               h22=( (q2(ic,jp,kc)+q2(ic,jc,kc))*(q2(ic,jp,kc)+q2(ic,jc,kc))&
        &           -(q2(ic,jm,kc)+q2(ic,jc,kc))*(q2(ic,jm,kc)+q2(ic,jc,kc)))*udx2
            endif

            if(kc.eq.1) then
               h23=( (q3(ic,jc,kp)+q3(ic,jm,kp))*(q2(ic,jc,kpp)+q2(ic,jc,kc)))*udx3m(kc)*0.25d0
            elseif(kc.eq.n3m) then
               h23=(-(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,kmm)))*udx3m(kc)*0.25d0
            else
               h23=( (q3(ic,jc,kp)+q3(ic,jm,kp))*(q2(ic,jc,kpp)+q2(ic,jc,kc))&
        &           -(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,kmm)))*udx3m(kc)*0.25d0
            endif

            dph(ic,jc,kc)=-(h21+h22+h23)

          enddo
        enddo
      enddo
      
      return
      end

