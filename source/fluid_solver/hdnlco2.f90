      subroutine hdnlco2
      use param
      use local_arrays, only: q1,q2,q3,co2,hco2
      use mpi_param, only: kstart,kend
      use walls_vars
      implicit none
      integer :: jc,kc,ic,js
      integer :: kpp,kmm,kp,jp,jmm,jpp,ip,imm,ipp
      real    :: h31,h32,h33,udx2,udx1

      udx1=dx1*0.5d0
      udx2=dx2*0.5d0
      do kc=kstart,kend
        kmm=kc-1
        kpp=kc+1
        kp=kc+1
        do jc=1,n2m
          jp=jc+1
          jmm=jc-1
          jpp=jc+1
          do ic=1,n1m
            ip=ic+1
            ipp=ic+1
            imm=ic-1
            if(ic.eq.1) then
               h31=( q1(ip,jc,kc)*(co2(ipp,jc,kc)+co2(ic,jc,kc)))*udx1
            elseif(ic.eq.n1m) then
               h31=(-q1(ic,jc,kc)*(co2(ic,jc,kc)+co2(imm,jc,kc)))*udx1
            else
               h31=( q1(ip,jc,kc)*(co2(ipp,jc,kc)+co2(ic,jc,kc)) &
        &           -q1(ic,jc,kc)*(co2(ic,jc,kc)+co2(imm,jc,kc)))*udx1
            endif

            if(jc.eq.1) then
               h32=( q2(ic,jp,kc)*(co2(ic,jpp,kc)+co2(ic,jc,kc)) &
        &           -q2ybs(ic,kc)*(co2(ic,jc,kc)+co2ybs(ic,kc)))*udx2
            elseif(jc.eq.n2m) then
               h32=( q2ybn(ic,kc)*(co2(ic,jc,kc)+co2ybn(ic,kc)) &
        &           -q2(ic,jc,kc)*(co2(ic,jc,kc)+co2(ic,jmm,kc)))*udx2
            else
               h32=( q2(ic,jp,kc)*(co2(ic,jpp,kc)+co2(ic,jc,kc)) &
        &           -q2(ic,jc,kc)*(co2(ic,jc,kc)+co2(ic,jmm,kc)))*udx2
            endif

            if(kc.eq.1) then
               h33=( (co2(ic,jc,kpp)*g3rm(kc)+co2(ic,jc,kc)*g3rm(kpp))&
        &              /(g3rm(kc)+g3rm(kpp))*q3(ic,jc,kp))*udx3m(kc)
            elseif(kc.eq.n3m) then
               h33=(-(co2(ic,jc,kc)*g3rm(kmm)+co2(ic,jc,kmm)*g3rm(kc))&
        &              /(g3rm(kc)+g3rm(kmm))*q3(ic,jc,kc)&
        &          )*udx3m(kc)
            else
               h33=( (co2(ic,jc,kpp)*g3rm(kc)+co2(ic,jc,kc)*g3rm(kpp))&
        &              /(g3rm(kc)+g3rm(kpp))*q3(ic,jc,kp)&
        &           -(co2(ic,jc,kc)*g3rm(kmm)+co2(ic,jc,kmm)*g3rm(kc))&
        &              /(g3rm(kc)+g3rm(kmm))*q3(ic,jc,kc)&
        &          )*udx3m(kc)
            endif

            hco2(ic,jc,kc)= -(h31+h32+h33)
     
          enddo
        enddo
      enddo

      return
      end

