      subroutine hdnlq1
      use param
      use local_arrays, only: q1,q2,q3,dq
      use mpi_param, only: kstart,kend
      use walls_vars
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip,js
      integer :: kmm,kpp
      real    :: h11,h12,h13,udx1,udx2

      udx1=dx1*0.25d0
      udx2=dx2*0.25d0

      do kc=kstart,kend
        kmm=kc-1
        kpp=kc+1
        kp=kc+1
        do jc=1,n2m
          jm=jc-1
          jp=jc+1
          do ic=2,n1m
            im=ic-1
            ip=ic+1

            if(ic.eq.2) then
               h11=( (q1(ip,jc,kc)+q1(ic,jc,kc))*(q1(ip,jc,kc)+q1(ic,jc,kc))&
        &           -(0.0d0+q1(ic,jc,kc))*(0.0d0+q1(ic,jc,kc)))*udx1
            elseif(ic.eq.n1m) then
               h11=( (0.0d0+q1(ic,jc,kc))*(0.0d0+q1(ic,jc,kc))&
        &           -(q1(im,jc,kc)+q1(ic,jc,kc))*(q1(im,jc,kc)+q1(ic,jc,kc)))*udx1
            else
               h11=( (q1(ip,jc,kc)+q1(ic,jc,kc))*(q1(ip,jc,kc)+q1(ic,jc,kc))&
        &           -(q1(im,jc,kc)+q1(ic,jc,kc))*(q1(im,jc,kc)+q1(ic,jc,kc)))*udx1
            endif

            if(jc.eq.1) then
               h12=( (q2(ic,jp,kc)+q2(im,jp,kc))*(q1(ic,jp,kc)+q1(ic,jc,kc))&
        &           -(q2ybs(ic,kc)+q2ybs(im,kc))*(q1(ic,jc,kc)+q1ybs(ic,kc)))*udx2
            elseif(jc.eq.n2m) then
               h12=( (q2ybn(ic,kc)+q2ybn(im,kc))*(q1ybn(ic,kc)+q1(ic,jc,kc))&
        &           -(q2(ic,jc,kc)+q2(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jm,kc)))*udx2
            else
               h12=( (q2(ic,jp,kc)+q2(im,jp,kc))*(q1(ic,jp,kc)+q1(ic,jc,kc))&
        &           -(q2(ic,jc,kc)+q2(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jm,kc)))*udx2
            endif

            if(kc.eq.1) then
               h13=( (q3(ic,jc,kp)+q3(im,jc,kp))*(q1(ic,jc,kpp)+q1(ic,jc,kc)))*udx3m(kc)*0.25d0
            elseif(kc.eq.n3m) then
               h13=( -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,kmm)))*udx3m(kc)*0.25d0
            else
               h13=( (q3(ic,jc,kp)+q3(im,jc,kp))*(q1(ic,jc,kpp)+q1(ic,jc,kc))&
        &           -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,kmm)))*udx3m(kc)*0.25d0
            endif

            dq(ic,jc,kc)=-(h11+h12+h13)

          enddo
        enddo
      enddo

      return
      end

