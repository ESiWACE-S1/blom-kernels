program advect_eitvel
  use wallclock_mod
  use random_mod
  implicit none
  integer, parameter :: nbdy = 3
#ifdef BLOM_CHANNEL_SMALL
  integer, parameter :: idm = 64
  integer, parameter :: jdm = 128
#elif (defined BLOM_CHANNEL_MEDIUM)
  integer, parameter :: idm = 128
  integer, parameter :: jdm = 256
#elif (defined BLOM_CHANNEL_LARGE)
  integer, parameter :: idm = 208
  integer, parameter :: jdm = 512
#endif
  integer, parameter :: kdm = 53

  integer :: l,i,j,k,km,kn,mm,nn, m, n, nstep
  integer :: ii,jj,kk, i0,j0
  integer, parameter :: ms = 10
  integer :: isu(1-nbdy:jdm+nbdy-1)
  integer :: ifu(-1:jdm+2, ms) ! first section index per section
  integer :: ilu(-1:jdm+2, ms) ! last section index per section
  integer :: isv(1-nbdy:jdm+nbdy-1)
  integer :: ifv(-1:jdm+2, ms) ! first section index per section
  integer :: ilv(-1:jdm+2, ms) ! last section index per section
  ! mod_state
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: utotm, scuy, umax, vtotm, scvx, vmax
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy, 2*kdm) :: u, umfltd, dpu, v, vmfltd, dpv
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy, 2) :: ubflxs_p, vbflxs_p, pbu, pbv
  real :: dlt, onemm, delt1
  real :: utotm_, vtotm_

  integer, parameter :: MAX_ITERATIONS=1000
  integer :: n_it
  integer, parameter :: ONE = 1
  !https://stackoverflow.com/a/6880672
  real(8)::t1,delta, delta_orig

  i0 = 0
  ii = idm
  j0 = 0
  jj = jdm
  kk = kdm
! --- letter 'm' refers to mid-time level (example: dp(i,j,km) )
! --- letter 'n' refers to old and new time level
  nstep = 1
  m=mod(nstep  ,2)+1
  n=mod(nstep+1,2)+1
  mm=(m-1)*kk
  nn=(n-1)*kk

  utotm=1.0; scuy=1.0; umax=1.0; vtotm=1.0; scvx=1.0; vmax=1.0
  u=1.0; umfltd=1.0; dpu=1.0; v=1.0; vmfltd=1.0; dpv=1.0
  ubflxs_p=1.0; vbflxs_p=1.0; pbu=1.0; pbv=1.0
  dlt=1.0; onemm=1.0; delt1=1.0
  write(*,*) "COMPUTE: eitvel"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "DEFAULT: iterating over",MAX_ITERATIONS, " iterations..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do j=-1,jj+2
     isu(j) = 1
     ifu(j,1) = -1
     ilu(j,1) = ii+2
     isv(j) = 1
     ifv(j,1) = -1
     ilv(j,1) = ii+2
  end do

  delta=0.0
  do n_it=1, MAX_ITERATIONS
     t1=wallclock()
     do k=1,kk
        km=k+mm
        kn=k+nn
        !
        ! --- --- advective and diffusive velocity at mid time level
        !
        do j=-1,jj+2
           do l=1,isu(j)
              do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
                 utotm(i,j)=u(i,j,km) &
                      +(ubflxs_p(i,j,m)*dlt/pbu(i,j,m)&
                      +umfltd(i,j,km)/max(onemm,dpu(i,j,kn)))&
                      /(delt1*scuy(i,j))
                 utotm(i,j)=max(-umax(i,j),min(umax(i,j),utotm(i,j)))
              enddo
           enddo
        enddo
        do j=0,jj+2
           do l=1,isv(j)
              do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
                 vtotm(i,j)=v(i,j,km)&
                      +(vbflxs_p(i,j,m)*dlt/pbv(i,j,m)&
                      +vmfltd(i,j,km)/max(onemm,dpv(i,j,kn)))&
                      /(delt1*scvx(i,j))
                 vtotm(i,j)=max(-vmax(i,j),min(vmax(i,j),vtotm(i,j)))
              enddo
           enddo
        enddo
        !
        !!     call remap_eitvel(scuy,scvx,scp2i,scp2,pbmin,   &
        !!                            pbu(1-nbdy,1-nbdy,n),pbv(1-nbdy,1-nbdy,n),  &
        !!                            p(1-nbdy,1-nbdy,k+1),utotm,vtotm,delt1,1, &
        !!                            dp(1-nbdy,1-nbdy,kn), &
        !!                            temp(1-nbdy,1-nbdy,kn), &
        !!                            saln(1-nbdy,1-nbdy,kn), &
        !!                            uflx(1-nbdy,1-nbdy,km), &
        !!                            vflx(1-nbdy,1-nbdy,km), &
        !!                            utflx(1-nbdy,1-nbdy,km), &
        !!                            vtflx(1-nbdy,1-nbdy,km), &
        !!                            usflx(1-nbdy,1-nbdy,km), &
        !!                            vsflx(1-nbdy,1-nbdy,km) &
        !!#ifdef TRC
        !!                           ,kn,trc &
        !!#endif
        !!                           )
        !
     enddo
     delta=delta+(wallclock()-t1)
  end do
  delta_orig=delta
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

  i = random_uniform(ONE,idm)
  j = random_uniform(ONE,jdm)
  write(*,*) 'random', i,j, utotm(i,j), vtotm(i,j)

#ifdef ADVECT_EITVEL_OPT1
  write(*,*) "OPT1: eitvel"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "full loops: iterating over",MAX_ITERATIONS, " iterations..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  delta=0.0
  do n_it=1, MAX_ITERATIONS
     t1=wallclock()
     do k=1,kk
        km=k+mm
        kn=k+nn
        !
        ! --- --- advective and diffusive velocity at mid time level
        !
        do j=-1,jj+2
           do i=0,ii+2
              utotm(i,j)=u(i,j,km) &
                   +(ubflxs_p(i,j,m)*dlt/pbu(i,j,m)&
                   +umfltd(i,j,km)/max(onemm,dpu(i,j,kn)))&
                   /(delt1*scuy(i,j))
              utotm(i,j)=max(-umax(i,j),min(umax(i,j),utotm(i,j)))
           enddo
        enddo
        do j=0,jj+2
           do i=-1,ii+2
              vtotm(i,j)=v(i,j,km)&
                   +(vbflxs_p(i,j,m)*dlt/pbv(i,j,m)&
                   +vmfltd(i,j,km)/max(onemm,dpv(i,j,kn)))&
                   /(delt1*scvx(i,j))
              vtotm(i,j)=max(-vmax(i,j),min(vmax(i,j),vtotm(i,j)))
           enddo
        enddo
        !
        !!     call remap_eitvel(scuy,scvx,scp2i,scp2,pbmin,   &
        !!                            pbu(1-nbdy,1-nbdy,n),pbv(1-nbdy,1-nbdy,n),  &
        !!                            p(1-nbdy,1-nbdy,k+1),utotm,vtotm,delt1,1, &
        !!                            dp(1-nbdy,1-nbdy,kn), &
        !!                            temp(1-nbdy,1-nbdy,kn), &
        !!                            saln(1-nbdy,1-nbdy,kn), &
        !!                            uflx(1-nbdy,1-nbdy,km), &
        !!                            vflx(1-nbdy,1-nbdy,km), &
        !!                            utflx(1-nbdy,1-nbdy,km), &
        !!                            vtflx(1-nbdy,1-nbdy,km), &
        !!                            usflx(1-nbdy,1-nbdy,km), &
        !!                            vsflx(1-nbdy,1-nbdy,km) &
        !!#ifdef TRC
        !!                           ,kn,trc &
        !!#endif
        !!                           )
        !
     enddo
     delta=delta+(wallclock()-t1)
  end do
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

  i = random_uniform(ONE,idm)
  j = random_uniform(ONE,jdm)
  write(*,*) 'random', i,j, utotm(i,j), vtotm(i,j)
#endif

#ifdef ADVECT_EITVEL_OPT3
  write(*,*) "OPT3: eitvel"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "full loops+scalars: iterating over",MAX_ITERATIONS, " iterations..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  delta=0.0
  do n_it=1, MAX_ITERATIONS
     t1=wallclock()
     do k=1,kk
        do j=-1,jj+2
           do i=0,ii+2
              utotm_=u(i,j,k+mm) &
                   +(ubflxs_p(i,j,m)*dlt/pbu(i,j,m)&
                   +umfltd(i,j,k+mm)/max(onemm,dpu(i,j,k+nn)))&
                   /(delt1*scuy(i,j))
              utotm_=max(-umax(i,j),min(umax(i,j),utotm_))
           enddo
           utotm(i,j) = utotm_
        enddo
        do j=0,jj+2
           do i=-1,ii+2
              vtotm_=v(i,j,k+mm)&
                   +(vbflxs_p(i,j,m)*dlt/pbv(i,j,m)&
                   +vmfltd(i,j,k+mm)/max(onemm,dpv(i,j,k+nn)))&
                   /(delt1*scvx(i,j))
              vtotm_=max(-vmax(i,j),min(vmax(i,j),vtotm_))
           enddo
           vtotm(i,j) = vtotm_
        enddo
     enddo
     delta=delta+(wallclock()-t1)
  end do
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

  i = random_uniform(ONE,idm)
  j = random_uniform(ONE,jdm)
  write(*,*) 'random', i,j, utotm(i,j), vtotm(i,j)
#endif

end program advect_eitvel
