      subroutine addbodyobj
      use param
      use mpih
      use local_arrays, only: q2,q3,q1,dsal,dens,co2,isbody1,isbody2,exhale_signal
      use mpi_param, only: kstart,kend
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic,im,jm,km
      real :: time_shift,breath_interval,time_signal,vel_peak,time_signal_exp,tprefactor,sprefactor,prefactor

      real :: injectedvol,injectmeanq1,injectmeanq2,injectmeanq3,injectmeandsal,injectmeandens,injectmeanco2

      !-- vel barrier
      do kc=kstart,kend
         do jc=1,n2
            do ic=1,n1
               km = kc-1
               if(kc.eq.1) km=kc
               q3(ic,jc,kc) = q3(ic,jc,kc)*((1.0-isbody1(ic,jc,km)-isbody2(ic,jc,km)) + (1.0-isbody1(ic,jc,km)-isbody2(ic,jc,km)) )*0.5
            enddo
         enddo
      enddo

      do kc=kstart,kend
         do jc=1,n2
            do ic=1,n1
               im = ic-1
               if(ic.eq.1) im=ic
               q1(ic,jc,kc) = q1(ic,jc,kc)*((1.0-isbody1(ic,jc,km)-isbody2(ic,jc,km)) + (1.0-isbody1(ic,jc,km)-isbody2(ic,jc,km)) )*0.5
            enddo
         enddo
      enddo

      do kc=kstart,kend
         do jc=1,n2
            do ic=1,n1
               jm = jc-1
               if(jc.eq.1) jm=jc
               q2(ic,jc,kc) = q2(ic,jc,kc)*((1.0-isbody1(ic,jc,km)-isbody2(ic,jc,km)) + (1.0-isbody1(ic,jc,km)-isbody2(ic,jc,km)) )*0.5
            enddo
         enddo
      enddo

      !-- dens and co2 boj
      do kc=kstart,kend
         do jc=1,n2
            do ic=1,n1
              if(isbody1(ic,jc,kc).eq.1.0 .or. isbody2(ic,jc,kc).eq.1.0) then
                     dens(ic,jc,kc) = 1.0d0
                     co2(ic,jc,kc) = 0.0d0
              endif
            enddo
         enddo
      enddo

      !-- sal boj
      do kc=kstartr,kendr
         do jc=1,n2r
            do ic=1,n1r
              if(isbody1(ic,jc,kc).eq.1.0 .or. isbody2(ic,jc,kc).eq.1.0) then
                     dsal(ic,jc,kc) = 0.0d0
              endif
            enddo
         enddo
      enddo

      ! compute injection quantities
      injectedvol  = 5e-4  /3.0/3.0/3.0  ! normalized 0.5L (by length scale 3m)
      injectmeanq1 = 0.0
      injectmeanq2 = -0.5*dcos(pi/3.0)/0.71    ! normalized 0.5m/s with angle (by free fall vel 0.71m/s)
      injectmeanq3 = -0.5*dsin(pi/3.0)/0.71    ! normalized 0.5m/s with angle (by free fall vel 0.71m/s)
      injectmeandsal = 1.0
      injectmeandens = 1.0
      injectmeanco2 = 1.0

      ! set temporal Gaussian func
      time_shift        =2.0/4.25      ! normalized 2s (by free fall time 4.25s)
      breath_interval   =4.25/4.25     ! normalized 4.25s (by free fall time 4.25s)
      
      time_signal_exp=exp(-0.5*( 2.0*(modulo(time,breath_interval)-time_shift)/kernel_width_time)**2)
      if(time_signal_exp.le.1d-8) time_signal_exp=0.d0
      tprefactor = (2.0/(2.0*pi)**0.5)/kernel_width_time
      time_signal = tprefactor*time_signal_exp

      !-- set breath
      do kc=kstart,kend
        do jc=1,n2
          do ic=1,n1
            !-- set region with breath only
            if (exhale_signal(ic,jc,kc) .ge. 1d-5) then

              !-- prefactor
              sprefactor = (2.0/(2.0*pi)**0.5)**3.0/kernel_width_space/kernel_width_space/kernel_width_space
              prefactor  = (sprefactor*exhale_signal(ic,jc,kc)*injectedvol)*(time_signal*ga*dt)

              !-- add to the quantities
              q1(ic,jc,kc) = q1(ic,jc,kc) + injectmeanq1*prefactor
              q2(ic,jc,kc) = q2(ic,jc,kc) + injectmeanq2*prefactor
              q3(ic,jc,kc) = q3(ic,jc,kc) + injectmeanq3*prefactor

              dsal(ic,jc,kc) = dsal(ic,jc,kc) + injectmeandsal*prefactor
              dens(ic,jc,kc) = dens(ic,jc,kc) + injectmeandens*prefactor
              co2(ic,jc,kc) = co2(ic,jc,kc) + injectmeanco2*prefactor



            endif
          enddo
        enddo
      enddo

      return
      end subroutine addbodyobj
!=====================================

      subroutine InitBreath
      use param
      use mpih
      use local_arrays, only: exhale_signal
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic
      real :: xcen1,xcen2,ycen,zcen,temp_signal

      exhale_signal(:,:,:)=0.d0

!      xcen=0.5d0*rext1
      xcen1=xpos1*rext1
      xcen2=xpos2*rext1
      ycen=0.472
      zcen=0.312
      do kc=kstart,kend
         do jc=1,n2
            do ic=1,n1
              exhale_signal(ic,jc,kc)=exp(-0.5*(  (2.0*(xc(ic)-xcen1)/kernel_width_space)**2  + (2.0*(yc(jc)-ycen)/kernel_width_space)**2+   (2.0*(zc(kc)-zcen)/kernel_width_space)**2))
              exhale_signal(ic,jc,kc)=exhale_signal(ic,jc,kc)+exp(-0.5*(  (2.0*(xc(ic)-xcen2)/kernel_width_space)**2  + (2.0*(yc(jc)-ycen)/kernel_width_space)**2+   (2.0*(zc(kc)-zcen)/kernel_width_space)**2))
            enddo
         enddo
      enddo

      return
      end subroutine InitBreath
