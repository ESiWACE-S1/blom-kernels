program remap_vel_u
  use wallclock_mod
  use random_mod
  implicit none
  integer, parameter :: nbdy = 3
  integer, parameter :: idm = 64
  integer, parameter :: jdm = 128
  integer, parameter :: kdm = 53
  
  integer :: l,i,j,iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne
  integer :: ii,jj,kk, i0,j0
  integer, parameter :: ms = 10
  integer :: isu(1-nbdy:jdm+nbdy)
  integer :: ifu(1-nbdy:jdm+nbdy, ms) ! first section index per section
  integer :: ilu(1-nbdy:jdm+nbdy, ms) ! last section index per section
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: dp, pup, plo
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fdu, fdv, ftu, ftv, fsu, fsv, cu, cv
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: dx, dy, xd, yd, tx, ty, temp, td, sx, sy, saln, sd
  real :: dxi, dyi, dpsw, dps, dpc, dpe, dpnw, dpn, dpne, dpse, dpw, dgmx, dfmx, dfmn, q, q1, q2, q3, q4, tgmx, tgmn, tfmx, tfmn
  real :: sgmx, sgmn, sfmx, sfmn
  integer :: mrg
  ! bigrid
  logical :: lperiodi, lperiodj, larctic
  integer :: nreg
  
  ! masks
  integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ip, iu, iv, iq
  integer :: ip_
  real :: p_
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: util1, util2, util3
  real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: depth
  
  integer, parameter :: MAX_ITERATIONS=10000
  integer :: n_it
  !https://stackoverflow.com/a/6880672
  real(8)::t1,delta, delta_orig
  
  real :: xd_, yd_, sx_, sy_, sd_, tx_, ty_, td_, dx_, dy_, dp_
  real :: dp_min, dp_max, t_min, t_max, s_min, s_max
  real :: te, tw, tn, ts, se, sw, sn, ss

  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: cuc, cvc, scp2, scp2i, pbu, uflx, utflx, usflx
  real :: xm, ym, xc0, xc1, x4, y4, dl, fd, qx, qy, x2, y2
  real :: a,ax,ay,axx,ayy,axy
  real :: cu_, cuc_, cvc_, scp2_, scp2i_, pbu_, uflx_, utflx_, usflx_, fdu_, ftu_, fsu_
  
  i0 = 0
  ii = idm
  j0 = 0
  jj = jdm
  kk = kdm
  mrg = 0
  ! bigrid init
  nreg = 1 ! periodic domain in i-index
  larctic=.false.
  lperiodi=.true.
  lperiodj=.false.
  depth=0.0
  depth(:,2-nbdy:jdm+nbdy)=1.0 ! dummy fixed depth (/=0)
  depth(:,1)=0.0 ! dummy fixed depth (/=0)
  depth(:,jdm)=0.0 ! dummy fixed depth (/=0)
  !call xctilr(depth,1,1, nbdy,nbdy, halo_ps)
!! --- is the domain periodic in i-index?
!      depmax=0.0
!      if     (i0+ii.eq.idm) then
!        do j= 1,jj
!          depmax=max(depmax,depth(ii,j))
!        enddo
!      endif
!      lperiodi=depmax.gt.0.0
!! --- is the domain periodic in j-index?
!      depmax=0.0
!      if     (j0+jj.eq.jdm) then
!        do i= 1,ii
!          depmax=max(depmax,depth(i,jj))
!        enddo
!      endif
!      larctic=depmax.gt.0.0 .and. nreg.eq.2
!      lperiodj=depmax.gt.0.0 .and. nreg.ne.2
  ! --- allow for non-periodic and non-arctic boundaries (part I).
  if     (.not.lperiodj .and. j0.eq.0) then
     ! ---   south boundary is all land.
     do j=1-nbdy,0
        do i=1-nbdy,ii+nbdy
           depth(i,j) = 0.0
        enddo
     enddo
  endif
  !
  if     (.not.lperiodj .and. .not.larctic .and. j0+jj.eq.jdm) then
     ! ---   north boundary is all land.
     do j=jj+1,jj+nbdy
        do i=1-nbdy,ii+nbdy
           depth(i,j) = 0.0
        enddo
     enddo
  endif
  !
  if     (.not.lperiodi .and. i0.eq.0) then
     ! ---   west boundary is all land.
     do j=1-nbdy,jj+nbdy
        do i=1-nbdy,0
           depth(i,j) = 0.0
        enddo
     enddo
  endif
  !
  if     (.not.lperiodi .and. i0+ii.eq.idm) then
     ! ---   east boundary is all land.
     do j=1-nbdy,jj+nbdy
        do i=ii+1,ii+nbdy
           depth(i,j) = 0.0
        enddo
     enddo
  endif
  ! --- start out with masks as land everywhere
  !$OMP PARALLEL DO PRIVATE(i)
  do j=1-nbdy,jdm+nbdy
     do i=1-nbdy,idm+nbdy
        ip(i,j)=0
        iq(i,j)=0
        iu(i,j)=0
        iv(i,j)=0
        util1(i,j)=0.
        util2(i,j)=0.
        util3(i,j)=0.
     enddo
  enddo
  !$OMP END PARALLEL DO
  !
  ! --- mass points are defined where water depth is greater than zero
  !$OMP PARALLEL DO PRIVATE(i)
  do j=1-nbdy,jj+nbdy
     do i=1-nbdy,ii+nbdy
        if (depth(i,j).gt.0.) then
           ip(i,j)=1
        endif
     enddo
  enddo
  !$OMP END PARALLEL DO
  !
  ! --- u,v points are located halfway between any 2 adjoining mass points
  ! --- 'interior' q points require water on all 4 sides.
  ! --- 'promontory' q points require water on 3 (or at least 2
  ! --- diametrically opposed) sides
  !$OMP PARALLEL DO PRIVATE(i)
  do j=1,jj
     do i=1,ii
        if (ip(i-1,j).gt.0.and.ip(i,j).gt.0) then
           iu(i,j)=1
        endif
        if (ip(i,j-1).gt.0.and.ip(i,j).gt.0) then
           iv(i,j)=1
        endif
        if (min(ip(i,j),ip(i-1,j),ip(i,j-1),ip(i-1,j-1)).gt.0) then
           iq(i,j)=1
        elseif ((ip(i  ,j).gt.0.and.ip(i-1,j-1).gt.0).or. &
           &            (ip(i-1,j).gt.0.and.ip(i  ,j-1).gt.0)    ) then
           iq(i,j)=1
        endif
        util1(i,j)=iu(i,j)
        util2(i,j)=iv(i,j)
        util3(i,j)=iq(i,j)
     enddo
  enddo
  !$OMP END PARALLEL DO
  if(lperiodi) then
     do j= 1-nbdy,jj+nbdy
        do i= 1,nbdy
           util1(i-nbdy,j) = util1(ii-nbdy+i,j)
           util2(i-nbdy,j) = util2(ii-nbdy+i,j)
           util3(i-nbdy,j) = util3(ii-nbdy+i,j)
           util1(ii+i,j) = util1(i,j)
           util2(ii+i,j) = util2(i,j)
           util3(ii+i,j) = util3(i,j)
        enddo
     end do
  end if
  !!call xctilr(util1,1,1, nbdy,nbdy, halo_us)
  !!call xctilr(util2,1,1, nbdy,nbdy, halo_vs)
  !!call xctilr(util3,1,1, nbdy,nbdy, halo_qs)
  !$OMP PARALLEL DO PRIVATE(i)
  do j= 1-nbdy,jj+nbdy
     do i= 1-nbdy,ii+nbdy
        iu(i,j)=nint(util1(i,j))
        iv(i,j)=nint(util2(i,j))
        iq(i,j)=nint(util3(i,j))
     enddo
  enddo
  !$OMP END PARALLEL DO
  !
  ! --- allow for non-periodic and non-arctic boundaries (part II).
  if     (.not.lperiodj .and. j0.eq.0) then
     ! ---   south boundary is all land.
     do j=1-nbdy,0
        do i=1-nbdy,ii+nbdy
           iq(i,j) = 0
           iu(i,j) = 0
           iv(i,j) = 0
        enddo
     enddo
  endif
  !
  if     (.not.lperiodj .and. .not.larctic .and. j0+jj.eq.jdm) then
     ! ---   north boundary is all land.
     do j=jj+1,jj+nbdy
        do i=1-nbdy,ii+nbdy
           iq(i,j) = 0
           iu(i,j) = 0
           iv(i,j) = 0
        enddo
     enddo
  endif
  !
  if     (.not.lperiodi .and. i0.eq.0) then
     ! ---   west boundary is all land.
     do j=1-nbdy,jj+nbdy
        do i=1-nbdy,0
           iq(i,j) = 0
           iu(i,j) = 0
           iv(i,j) = 0
        enddo
     enddo
  endif
  !
  if     (.not.lperiodi .and. i0+ii.eq.idm) then
     ! ---   east boundary is all land.
     do j=1-nbdy,jj+nbdy
        do i=ii+1,ii+nbdy
           iq(i,j) = 0
           iu(i,j) = 0
           iv(i,j) = 0
        enddo
     enddo
  endif

  write(*,*) 'nreg',nreg
  write(*,*) 'lperiodi',lperiodi
  write(*,*) 'lperiodj',lperiodj
  write(*,*) 'larctic',larctic
  ! --- - u-components of fluxes.
  !
  write(*,*) "COMPUTE: remap_limited_grad"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "DEFAULT: iterating over",MAX_ITERATIONS, " iterations..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call indxi(iu,ifu,ilu,isu)

  call init_seed
  do j=1-nbdy,jdm+nbdy
     do i=1-nbdy,idm+nbdy
       dp(i,j) = random_uniform(-1.0, 1.0)
       pup(i,j) = random_uniform(-1.0, 1.0)
       plo(i,j) = random_uniform(-1.0, 1.0)
       fdu(i,j) = random_uniform(-1.0, 1.0)
       fdv(i,j) = random_uniform(-1.0, 1.0)
       ftu(i,j) = random_uniform(-1.0, 1.0)
       ftv(i,j) = random_uniform(-1.0, 1.0)
       fsu(i,j) = random_uniform(-1.0, 1.0)
       fsv(i,j) = random_uniform(-1.0, 1.0)
       cu(i,j) = random_uniform(-1.0, 1.0)
       cv(i,j) = random_uniform(-1.0, 1.0)
       dx(i,j) = random_uniform(-1.0, 1.0)
       dy(i,j) = random_uniform(-1.0, 1.0)
       xd(i,j) = random_uniform(-1.0, 1.0)
       yd(i,j) = random_uniform(-1.0, 1.0)
       tx(i,j) = random_uniform(-1.0, 1.0)
       ty(i,j) = random_uniform(-1.0, 1.0)
       temp(i,j) = random_uniform(-1.0, 1.0)
       td(i,j) = random_uniform(-1.0, 1.0)
       sx(i,j) = random_uniform(-1.0, 1.0)
       sy(i,j) = random_uniform(-1.0, 1.0)
       saln(i,j) = random_uniform(-1.0, 1.0)
       sd(i,j) = random_uniform(-1.0, 1.0)
       cuc(i,j) = random_uniform(-1.0, 1.0)
       cvc(i,j) = random_uniform(-1.0, 1.0)
       scp2(i,j) = random_uniform(-1.0, 1.0)
       scp2i(i,j) = random_uniform(-1.0, 1.0)
       pbu(i,j) = random_uniform(-1.0, 1.0)
       uflx(i,j) = random_uniform(-1.0, 1.0)
       utflx(i,j) = random_uniform(-1.0, 1.0)
       usflx(i,j) = random_uniform(-1.0, 1.0)
    end do
 end do
  
  delta=0.0
  do n_it=1, MAX_ITERATIONS
     t1=wallclock()
     do j=1-mrg,jj+mrg
        !
        do l=1,isu(j)
           do i=max(1-mrg,ifu(j,l)),min(ii+mrg+1,ilu(j,l))
              !
              ! --- --- Assuming coordinate [0,0] at the u-point, the non-dimensional
              ! --- --- fluxing area is defined as the area of a polygon with vertices
              ! --- --- [0,1/2], [-cuc(i,j+1),-cvc(i,j+1)+1/2], [xm,ym],
              ! --- --- [-cuc(i,j),-cvc(i,j)-1/2], and [0,-1/2]. The vertex [xm,ym] is
              ! --- --- defined so that the polygon area is equal to cu(i,j).
              !
              ym=-.5*(cvc(i,j)+cvc(i,j+1))
              xm=((ym+.5)*cuc(i,j)-(ym-.5)*cuc(i,j+1)-2.*cu(i,j)) &
                     /(1.+cvc(i,j)-cvc(i,j+1))
              !
              if (cu(i,j).gt.0.) then
                 !
                 if (cvc(i,j).gt.0.) then
                    !
                    ! --- ------- Add contributions from grid cell (i-1,j-1). Assuming
                    ! --- ------- coordinate [0,0] at the cell center, the contributions are
                    ! --- ------- flux integrals over the triangle with vertices
                    ! --- ------- [xc1+1/2,1/2], [-cuc(i,j)+1/2,-cvc(i,j)+1/2], and
                    ! --- ------- [1/2,1/2].
                    !
                    xc0=(xm*cvc(i,j)-cuc(i,j)*(ym+.5))/(cvc(i,j)+ym+.5)
                    xc1=xc0*scp2(i-1,j)*scp2i(i-1,j-1)
                    x4=xc0+.5
                    y4=-.5
                    call triint(scp2(i-1,j-1), &
                              xc1+.5,.5,-cuc(i,j)+.5,-cvc(i,j)+.5,.5,.5, &
                                        a,ax,ay,axx,ayy,axy &
#if defined(TRC) && defined(ATRC)
                                       ,axxx,ayyy,axxy,axyy &
#endif
                                       )
                    dl=min(dp(i-1,j-1),max(0.,pbu(i,j)-pup(i-1,j-1)))
                    fd=a*dl+ax*dx(i-1,j-1)+ay*dy(i-1,j-1)
                    fdu(i,j)=fdu(i,j)+fd
                    qx=ax*dl+axx*dx(i-1,j-1)+axy*dy(i-1,j-1)
                    qy=ay*dl+axy*dx(i-1,j-1)+ayy*dy(i-1,j-1)
                    ftu(i,j)=ftu(i,j)+fd*td(i-1,j-1) &
                                    +qx*tx(i-1,j-1)+qy*ty(i-1,j-1)
                    fsu(i,j)=fsu(i,j)+fd*sd(i-1,j-1) &
                                    +qx*sx(i-1,j-1)+qy*sy(i-1,j-1)
#ifdef TRC
#  ifdef ATRC
                    qxx=axx*dl+axxx*dx(i-1,j-1)+axxy*dy(i-1,j-1)
                    qyy=ayy*dl+axyy*dx(i-1,j-1)+ayyy*dy(i-1,j-1)
                    qxy=axy*dl+axxy*dx(i-1,j-1)+axyy*dy(i-1,j-1)
                    do nt=1,natr
                       fdt=fd*trd(nt,i-1,j-1) &
                                    +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                       ftru(nt,i,j)=ftru(nt,i,j)+fdt
                       fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i-1,j-1) &
                                             +(qx *trd(nt,i-1,j-1)  &
                                              +qxx*trx(nt,i-1,j-1)  &
                                              +qxy*try(nt,i-1,j-1))*agx(nt,i-1,j-1) &
                                             +(qy *trd(nt,i-1,j-1)  &
                                              +qxy*trx(nt,i-1,j-1)  &
                                              +qyy*try(nt,i-1,j-1))*agy(nt,i-1,j-1)
                    enddo
                    do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                       if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                       ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j-1) &
                                             +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                    enddo
#  else
                    do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                       if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                       ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j-1) &
                                             +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                    enddo
#  endif
#endif
                 else
                    x4=-cuc(i,j)+.5
                    y4=-cvc(i,j)-.5
                 endif
                 !
                 if (cvc(i,j+1).lt.0.) then
                    !
                    ! --- ------- Add contributions from grid cell (i-1,j+1). Assuming
                    ! --- ------- coordinate [0,0] at the cell center, the contributions are
                    ! --- ------- flux integrals over the triangle with vertices
                    ! --- ------- [xc1+1/2,-1/2], [1/2,-1/2], and
                    ! --- ------- [-cuc(i,j+1)+1/2,-cvc(i,j+1)-1/2].
                    !
                    xc0=(xm*cvc(i,j+1)-cuc(i,j+1)*(ym-.5))/(cvc(i,j+1)+ym-.5)
                    xc1=xc0*scp2(i-1,j)*scp2i(i-1,j+1)
                    x2=xc0+.5
                    y2=.5
                    call triint(scp2(i-1,j+1), &
                              xc1+.5,-.5,.5,-.5,-cuc(i,j+1)+.5,-cvc(i,j+1)-.5,&
                                        a,ax,ay,axx,ayy,axy &
#if defined(TRC) && defined(ATRC)
                                       ,axxx,ayyy,axxy,axyy &
#endif
                                       )
                    dl=min(dp(i-1,j+1),max(0.,pbu(i,j)-pup(i-1,j+1)))
                    fd=a*dl+ax*dx(i-1,j+1)+ay*dy(i-1,j+1)
                    fdu(i,j)=fdu(i,j)+fd
                    qx=ax*dl+axx*dx(i-1,j+1)+axy*dy(i-1,j+1)
                    qy=ay*dl+axy*dx(i-1,j+1)+ayy*dy(i-1,j+1)
                    ftu(i,j)=ftu(i,j)+fd*td(i-1,j+1) &
                                    +qx*tx(i-1,j+1)+qy*ty(i-1,j+1)
                    fsu(i,j)=fsu(i,j)+fd*sd(i-1,j+1) &
                                    +qx*sx(i-1,j+1)+qy*sy(i-1,j+1)
#ifdef TRC
#  ifdef ATRC
                    qxx=axx*dl+axxx*dx(i-1,j+1)+axxy*dy(i-1,j+1)
                    qyy=ayy*dl+axyy*dx(i-1,j+1)+ayyy*dy(i-1,j+1)
                    qxy=axy*dl+axxy*dx(i-1,j+1)+axyy*dy(i-1,j+1)
                    do nt=1,natr
                       fdt=fd*trd(nt,i-1,j+1) &
                                    +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
                       ftru(nt,i,j)=ftru(nt,i,j)+fdt
                       fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i-1,j+1)& 
                                             +(qx *trd(nt,i-1,j+1) &
                                              +qxx*trx(nt,i-1,j+1) &
                                              +qxy*try(nt,i-1,j+1))*agx(nt,i-1,j+1) &
                                             +(qy *trd(nt,i-1,j+1) &
                                              +qxy*trx(nt,i-1,j+1) &
                                              +qyy*try(nt,i-1,j+1))*agy(nt,i-1,j+1)
                    enddo
                    do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                       if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                       ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j+1) &
                                             +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
                    enddo
#  else
                    do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                       if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                       ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j+1) &
                                             +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
                    enddo
#  endif
#endif
                 else
                    x2=-cuc(i,j+1)+.5
                    y2=-cvc(i,j+1)+.5
                 endif
                 !
                 ! --- ----- Add contributions from grid cell (i-1,j). Assuming
                 ! --- ----- coordinate [0,0] at the cell center, the contributions are
                 ! --- ----- flux integrals over the pentagon with vertices [1/2,1/2],
                 ! --- ----- [x2,y2], [xm+1/2,ym], [x4,y4], and [1/2,-1/2].
                 !
                 call penint(scp2(i-1,j), &
                         .5,.5,x2,y2,xm+.5,ym,x4,y4,.5,-.5, &
                                   a,ax,ay,axx,ayy,axy &
#if defined(TRC) && defined(ATRC)
                                  ,axxx,ayyy,axxy,axyy &
#endif
                                  )
                 dl=min(dp(i-1,j),max(0.,pbu(i,j)-pup(i-1,j)))
                 fd=a*dl+ax*dx(i-1,j)+ay*dy(i-1,j)
                 fdu(i,j)=fdu(i,j)+fd
                 qx=ax*dl+axx*dx(i-1,j)+axy*dy(i-1,j)
                 qy=ay*dl+axy*dx(i-1,j)+ayy*dy(i-1,j)
                 ftu(i,j)=ftu(i,j)+fd*td(i-1,j) &
                               +qx*tx(i-1,j)+qy*ty(i-1,j)
                 fsu(i,j)=fsu(i,j)+fd*sd(i-1,j) &
                               +qx*sx(i-1,j)+qy*sy(i-1,j)
#ifdef TRC
#  ifdef ATRC
                 qxx=axx*dl+axxx*dx(i-1,j)+axxy*dy(i-1,j)
                 qyy=ayy*dl+axyy*dx(i-1,j)+ayyy*dy(i-1,j)
                 qxy=axy*dl+axxy*dx(i-1,j)+axyy*dy(i-1,j)
                 do nt=1,natr
                    fdt=fd*trd(nt,i-1,j) &
                               +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                    ftru(nt,i,j)=ftru(nt,i,j)+fdt
                    fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i-1,j) &
                                        +(qx *trd(nt,i-1,j) &
                                         +qxx*trx(nt,i-1,j) &
                                         +qxy*try(nt,i-1,j))*agx(nt,i-1,j) &
                                        +(qy *trd(nt,i-1,j) &
                                         +qxy*trx(nt,i-1,j) &
                                         +qyy*try(nt,i-1,j))*agy(nt,i-1,j)
                 enddo
                 do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                    if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                    ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j) &
                                        +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                 enddo
#  else
                 do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                    if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                    ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i-1,j) &
                                        +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                 enddo
#  endif
#endif
                 !
              else
                 !
                 if (cvc(i,j).gt.0.) then
                    !
                    ! --- ------- Add contributions from grid cell (i,j-1). Assuming
                    ! --- ------- coordinate [0,0] at the cell center, the contributions are
                    ! --- ------- flux integrals over the triangle with vertices
                    ! --- ------- [xc1-1/2,1/2], [-cuc(i,j)-1/2,-cvc(i,j)+1/2], and
                    ! --- ------- [-1/2,1/2].
                    !
                    xc0=(xm*cvc(i,j)-cuc(i,j)*(ym+.5))/(cvc(i,j)+ym+.5)
                    xc1=xc0*scp2(i,j)*scp2i(i,j-1)
                    x4=xc0-.5
                    y4=-.5
                    call triint(scp2(i,j-1), &
                              xc1-.5,.5,-cuc(i,j)-.5,-cvc(i,j)+.5,-.5,.5, &
                                        a,ax,ay,axx,ayy,axy &
#if defined(TRC) && defined(ATRC)
                                       ,axxx,ayyy,axxy,axyy &
#endif
                                       )
                    dl=min(dp(i,j-1),max(0.,pbu(i,j)-pup(i,j-1)))
                    fd=a*dl+ax*dx(i,j-1)+ay*dy(i,j-1)
                    fdu(i,j)=fdu(i,j)+fd
                    qx=ax*dl+axx*dx(i,j-1)+axy*dy(i,j-1)
                    qy=ay*dl+axy*dx(i,j-1)+ayy*dy(i,j-1)
                    ftu(i,j)=ftu(i,j)+fd*td(i,j-1) &
                                    +qx*tx(i,j-1)+qy*ty(i,j-1)
                    fsu(i,j)=fsu(i,j)+fd*sd(i,j-1) &
                                    +qx*sx(i,j-1)+qy*sy(i,j-1)
#ifdef TRC
#  ifdef ATRC
                    qxx=axx*dl+axxx*dx(i,j-1)+axxy*dy(i,j-1)
                    qyy=ayy*dl+axyy*dx(i,j-1)+ayyy*dy(i,j-1)
                    qxy=axy*dl+axxy*dx(i,j-1)+axyy*dy(i,j-1)
                    do nt=1,natr
                       fdt=fd*trd(nt,i,j-1) &
                                    +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                       ftru(nt,i,j)=ftru(nt,i,j)+fdt
                       fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i,j-1) &
                                             +(qx *trd(nt,i,j-1)  &
                                              +qxx*trx(nt,i,j-1)  &
                                              +qxy*try(nt,i,j-1))*agx(nt,i,j-1) &
                                             +(qy *trd(nt,i,j-1)  &
                                              +qxy*trx(nt,i,j-1)  &
                                              +qyy*try(nt,i,j-1))*agy(nt,i,j-1)
                    enddo
                    do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                       if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                       ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j-1) &
                                             +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                    enddo
#  else
                    do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                       if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                       ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j-1) &
                                             +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                    enddo
#  endif
#endif
                 else
                    x4=-cuc(i,j)-.5
                    y4=-cvc(i,j)-.5
                 endif
                 !
                 if (cvc(i,j+1).lt.0.) then
                    !
                    ! --- ------- Add contributions from grid cell (i,j+1). Assuming
                    ! --- ------- coordinate [0,0] at the cell center, the contributions are
                    ! --- ------- flux integrals over the triangle with vertices
                    ! --- ------- [xc1-1/2,-1/2], [-1/2,-1/2], and
                    ! --- ------- [-cuc(i,j+1)-1/2,-cvc(i,j+1)-1/2].
                    !
                    xc0=(xm*cvc(i,j+1)-cuc(i,j+1)*(ym-.5))/(cvc(i,j+1)+ym-.5)
                    xc1=xc0*scp2(i,j)*scp2i(i,j+1)
                    x2=xc0-.5
                    y2=.5
                    call triint(scp2(i,j+1), &
                              xc1-.5,-.5,-.5,-.5,-cuc(i,j+1)-.5,-cvc(i,j+1)-.5, &
                                        a,ax,ay,axx,ayy,axy &
#if defined(TRC) && defined(ATRC)
                                       ,axxx,ayyy,axxy,axyy &
#endif
                                       )
                    dl=min(dp(i,j+1),max(0.,pbu(i,j)-pup(i,j+1)))
                    fd=a*dl+ax*dx(i,j+1)+ay*dy(i,j+1)
                    fdu(i,j)=fdu(i,j)+fd
                    qx=ax*dl+axx*dx(i,j+1)+axy*dy(i,j+1)
                    qy=ay*dl+axy*dx(i,j+1)+ayy*dy(i,j+1)
                    ftu(i,j)=ftu(i,j)+fd*td(i,j+1) &
                                    +qx*tx(i,j+1)+qy*ty(i,j+1)
                    fsu(i,j)=fsu(i,j)+fd*sd(i,j+1) &
                                    +qx*sx(i,j+1)+qy*sy(i,j+1)
#ifdef TRC
#  ifdef ATRC
                    qxx=axx*dl+axxx*dx(i,j+1)+axxy*dy(i,j+1)
                    qyy=ayy*dl+axyy*dx(i,j+1)+ayyy*dy(i,j+1)
                    qxy=axy*dl+axxy*dx(i,j+1)+axyy*dy(i,j+1)
                    do nt=1,natr
                       fdt=fd*trd(nt,i,j+1) &
                                    +qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
                       ftru(nt,i,j)=ftru(nt,i,j)+fdt
                       fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i,j+1) &
                                             +(qx *trd(nt,i,j+1)  &
                                              +qxx*trx(nt,i,j+1)  &
                                              +qxy*try(nt,i,j+1))*agx(nt,i,j+1) &
                                             +(qy *trd(nt,i,j+1)  &
                                              +qxy*trx(nt,i,j+1)  &
                                              +qyy*try(nt,i,j+1))*agy(nt,i,j+1)
                    enddo
                    do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                       if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                       ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j+1) &
                                             +qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
                    enddo
#  else
                    do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                       if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                       ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j+1) &
                                             +qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
                    enddo
#  endif
#endif
                 else
                    x2=-cuc(i,j+1)-.5
                    y2=-cvc(i,j+1)+.5
                 endif
                 !
                 ! --- ----- Add contributions from grid cell (i,j). Assuming
                 ! --- ----- coordinate [0,0] at the cell center, the contributions are
                 ! --- ----- flux integrals over the pentagon with vertices [-1/2,1/2],
                 ! --- ----- [x2,y2], [xm-1/2,ym], [x4,y4], and [-1/2,-1/2].
                 !
                 call penint(scp2(i,j), &
                         -.5,.5,x2,y2,xm-.5,ym,x4,y4,-.5,-.5, &
                                   a,ax,ay,axx,ayy,axy &
#if defined(TRC) && defined(ATRC)
                                  ,axxx,ayyy,axxy,axyy &
#endif
                                  )
                 dl=min(dp(i,j),max(0.,pbu(i,j)-pup(i,j)))
                 fd=a*dl+ax*dx(i,j)+ay*dy(i,j)
                 fdu(i,j)=fdu(i,j)+fd
                 qx=ax*dl+axx*dx(i,j)+axy*dy(i,j)
                 qy=ay*dl+axy*dx(i,j)+ayy*dy(i,j)
                 ftu(i,j)=ftu(i,j)+fd*td(i,j) &
                               +qx*tx(i,j)+qy*ty(i,j)
                 fsu(i,j)=fsu(i,j)+fd*sd(i,j) &
                               +qx*sx(i,j)+qy*sy(i,j)
#ifdef TRC
#  ifdef ATRC
                 qxx=axx*dl+axxx*dx(i,j)+axxy*dy(i,j)
                 qyy=ayy*dl+axyy*dx(i,j)+ayyy*dy(i,j)
                 qxy=axy*dl+axxy*dx(i,j)+axyy*dy(i,j)
                 do nt=1,natr
                    fdt=fd*trd(nt,i,j) &
                               +qx*trx(nt,i,j)+qy*try(nt,i,j)
                    ftru(nt,i,j)=ftru(nt,i,j)+fdt
                    fagu(nt,i,j)=fagu(nt,i,j)+fdt*agd(nt,i,j) &
                                        +(qx *trd(nt,i,j) &
                                         +qxx*trx(nt,i,j) &
                                         +qxy*try(nt,i,j))*agx(nt,i,j) &
                                        +(qy *trd(nt,i,j) &
                                         +qxy*trx(nt,i,j) &
                                         +qyy*try(nt,i,j))*agy(nt,i,j)
                 enddo
                 do nt=natr+1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
                    if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                    ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j) &
                                        +qx*trx(nt,i,j)+qy*try(nt,i,j)
                 enddo
#  else
                 do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
                    if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
                    ftru(nt,i,j)=ftru(nt,i,j)+fd*trd(nt,i,j) &
                                        +qx*trx(nt,i,j)+qy*try(nt,i,j)
                 enddo
#  endif
#endif
                 !
              endif
              !
              ! --- --- u-component of mass, heat and salt flux.
              uflx(i,j)=fdu(i,j)
              utflx(i,j)=ftu(i,j)
              usflx(i,j)=fsu(i,j)
              !
           enddo
        enddo
        !
     enddo
     delta=delta+(wallclock()-t1)
  enddo
  delta_orig=delta
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

#define REMAP_VEL_U_OPT1
#ifdef REMAP_VEL_U_OPT1
  write(*,*) "full loops: iterating over",MAX_ITERATIONS, " iterations..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call init_seed
  do j=1-nbdy,jdm+nbdy
     do i=1-nbdy,idm+nbdy
       dp(i,j) = random_uniform(-1.0, 1.0)
       pup(i,j) = random_uniform(-1.0, 1.0)
       plo(i,j) = random_uniform(-1.0, 1.0)
       fdu(i,j) = random_uniform(-1.0, 1.0)
       fdv(i,j) = random_uniform(-1.0, 1.0)
       ftu(i,j) = random_uniform(-1.0, 1.0)
       ftv(i,j) = random_uniform(-1.0, 1.0)
       fsu(i,j) = random_uniform(-1.0, 1.0)
       fsv(i,j) = random_uniform(-1.0, 1.0)
       cu(i,j) = random_uniform(-1.0, 1.0)
       cv(i,j) = random_uniform(-1.0, 1.0)
       dx(i,j) = random_uniform(-1.0, 1.0)
       dy(i,j) = random_uniform(-1.0, 1.0)
       xd(i,j) = random_uniform(-1.0, 1.0)
       yd(i,j) = random_uniform(-1.0, 1.0)
       tx(i,j) = random_uniform(-1.0, 1.0)
       ty(i,j) = random_uniform(-1.0, 1.0)
       temp(i,j) = random_uniform(-1.0, 1.0)
       td(i,j) = random_uniform(-1.0, 1.0)
       sx(i,j) = random_uniform(-1.0, 1.0)
       sy(i,j) = random_uniform(-1.0, 1.0)
       saln(i,j) = random_uniform(-1.0, 1.0)
       sd(i,j) = random_uniform(-1.0, 1.0)
       cuc(i,j) = random_uniform(-1.0, 1.0)
       cvc(i,j) = random_uniform(-1.0, 1.0)
       scp2(i,j) = random_uniform(-1.0, 1.0)
       scp2i(i,j) = random_uniform(-1.0, 1.0)
       pbu(i,j) = random_uniform(-1.0, 1.0)
       uflx(i,j) = random_uniform(-1.0, 1.0)
       utflx(i,j) = random_uniform(-1.0, 1.0)
       usflx(i,j) = random_uniform(-1.0, 1.0)
    end do
 end do
   
  delta=0.0
  do n_it=1, MAX_ITERATIONS
     t1=wallclock()
     do j=1-mrg,jj+mrg
        do i=1-mrg,ii+mrg+1
           if(iu(i,j) /= 0) then
              fdu_ = 0.0
              ftu_ = 0.0
              fsu_ = 0.0
              cu_ = cu(i,j)
              cuc_ = cuc(i,j)
              cvc_ = cvc(i,j)
              ym=-.5*(cvc_+cvc(i,j+1))
              xm=((ym+.5)*cuc_-(ym-.5)*cuc(i,j+1)-2.*cu_) &
                     /(1.+cvc_-cvc(i,j+1))
              if (cu_.gt.0.) then
                 if (cvc_.gt.0.) then ! Add contributions from grid cell (i-1,j-1)
                    ! --- ------- flux integrals over the triangle with vertices
                    ! --- ------- [xc1+1/2,1/2], [-cuc(i,j)+1/2,-cvc(i,j)+1/2], and
                    ! --- ------- [1/2,1/2].
                    xc0=(xm*cvc_-cuc_*(ym+.5))/(cvc_+ym+.5)
                    xc1=xc0*scp2(i-1,j)*scp2i(i-1,j-1)
                    x4=xc0+.5
                    y4=-.5
                    call triint(scp2(i-1,j-1), &
                              xc1+.5,.5,-cuc_+.5,-cvc_+.5,.5,.5, &
                                        a,ax,ay,axx,ayy,axy &
                                       )
                    dl=min(dp(i-1,j-1),max(0.,pbu(i,j)-pup(i-1,j-1)))
                    fd=a*dl+ax*dx(i-1,j-1)+ay*dy(i-1,j-1)
                    fdu_=fdu_+fd
                    qx=ax*dl+axx*dx(i-1,j-1)+axy*dy(i-1,j-1)
                    qy=ay*dl+axy*dx(i-1,j-1)+ayy*dy(i-1,j-1)
                    ftu_=ftu_+fd*td(i-1,j-1) &
                                    +qx*tx(i-1,j-1)+qy*ty(i-1,j-1)
                    fsu_=fsu_+fd*sd(i-1,j-1) &
                                    +qx*sx(i-1,j-1)+qy*sy(i-1,j-1)
                 else
                    x4=-cuc_+.5
                    y4=-cvc_-.5
                 endif
                 !
                 if (cvc(i,j+1).lt.0.) then ! --- ------- Add contributions from grid cell (i-1,j+1). Assuming
                    xc0=(xm*cvc(i,j+1)-cuc(i,j+1)*(ym-.5))/(cvc(i,j+1)+ym-.5)
                    xc1=xc0*scp2(i-1,j)*scp2i(i-1,j+1)
                    x2=xc0+.5
                    y2=.5
                    call triint(scp2(i-1,j+1), &
                              xc1+.5,-.5,.5,-.5,-cuc(i,j+1)+.5,-cvc(i,j+1)-.5,&
                                        a,ax,ay,axx,ayy,axy &
                                       )
                    dl=min(dp(i-1,j+1),max(0.,pbu(i,j)-pup(i-1,j+1)))
                    fd=a*dl+ax*dx(i-1,j+1)+ay*dy(i-1,j+1)
                    fdu_=fdu_+fd
                    qx=ax*dl+axx*dx(i-1,j+1)+axy*dy(i-1,j+1)
                    qy=ay*dl+axy*dx(i-1,j+1)+ayy*dy(i-1,j+1)
                    ftu_=ftu_+fd*td(i-1,j+1) &
                                    +qx*tx(i-1,j+1)+qy*ty(i-1,j+1)
                    fsu_=fsu_+fd*sd(i-1,j+1) &
                                    +qx*sx(i-1,j+1)+qy*sy(i-1,j+1)
                 else
                    x2=-cuc(i,j+1)+.5
                    y2=-cvc(i,j+1)+.5
                 endif
                 ! --- ----- Add contributions from grid cell (i-1,j). Assuming
                 call penint(scp2(i-1,j), &
                         .5,.5,x2,y2,xm+.5,ym,x4,y4,.5,-.5, &
                                   a,ax,ay,axx,ayy,axy &
                                  )
                 dl=min(dp(i-1,j),max(0.,pbu(i,j)-pup(i-1,j)))
                 fd=a*dl+ax*dx(i-1,j)+ay*dy(i-1,j)
                 fdu_=fdu_+fd
                 qx=ax*dl+axx*dx(i-1,j)+axy*dy(i-1,j)
                 qy=ay*dl+axy*dx(i-1,j)+ayy*dy(i-1,j)
                 ftu_=ftu_+fd*td(i-1,j) &
                               +qx*tx(i-1,j)+qy*ty(i-1,j)
                 fsu_=fsu_+fd*sd(i-1,j) &
                               +qx*sx(i-1,j)+qy*sy(i-1,j)
              else ! cu_ < 0
                 if (cvc_.gt.0.) then ! --- ------- Add contributions from grid cell (i,j-1).
                    xc0=(xm*cvc_-cuc_*(ym+.5))/(cvc_+ym+.5)
                    xc1=xc0*scp2(i,j)*scp2i(i,j-1)
                    x4=xc0-.5
                    y4=-.5
                    call triint(scp2(i,j-1), &
                              xc1-.5,.5,-cuc_-.5,-cvc_+.5,-.5,.5, &
                                        a,ax,ay,axx,ayy,axy &
                                       )
                    dl=min(dp(i,j-1),max(0.,pbu(i,j)-pup(i,j-1)))
                    fd=a*dl+ax*dx(i,j-1)+ay*dy(i,j-1)
                    fdu_=fdu_+fd
                    qx=ax*dl+axx*dx(i,j-1)+axy*dy(i,j-1)
                    qy=ay*dl+axy*dx(i,j-1)+ayy*dy(i,j-1)
                    ftu_=ftu_+fd*td(i,j-1) &
                                    +qx*tx(i,j-1)+qy*ty(i,j-1)
                    fsu_=fsu_+fd*sd(i,j-1) &
                                    +qx*sx(i,j-1)+qy*sy(i,j-1)
                 else
                    x4=-cuc_-.5
                    y4=-cvc_-.5
                 endif
                 if (cvc(i,j+1).lt.0.) then ! --- ------- Add contributions from grid cell (i,j+1). Assuming
                    xc0=(xm*cvc(i,j+1)-cuc(i,j+1)*(ym-.5))/(cvc(i,j+1)+ym-.5)
                    xc1=xc0*scp2(i,j)*scp2i(i,j+1)
                    x2=xc0-.5
                    y2=.5
                    call triint(scp2(i,j+1), &
                              xc1-.5,-.5,-.5,-.5,-cuc(i,j+1)-.5,-cvc(i,j+1)-.5, &
                                        a,ax,ay,axx,ayy,axy &
                                       )
                    dl=min(dp(i,j+1),max(0.,pbu(i,j)-pup(i,j+1)))
                    fd=a*dl+ax*dx(i,j+1)+ay*dy(i,j+1)
                    fdu_=fdu_+fd
                    qx=ax*dl+axx*dx(i,j+1)+axy*dy(i,j+1)
                    qy=ay*dl+axy*dx(i,j+1)+ayy*dy(i,j+1)
                    ftu_=ftu_+fd*td(i,j+1) &
                                    +qx*tx(i,j+1)+qy*ty(i,j+1)
                    fsu_=fsu_+fd*sd(i,j+1) &
                                    +qx*sx(i,j+1)+qy*sy(i,j+1)
                 else
                    x2=-cuc(i,j+1)-.5
                    y2=-cvc(i,j+1)+.5
                 endif
                 ! --- ----- Add contributions from grid cell (i,j). Assuming
                 call penint(scp2(i,j), &
                         -.5,.5,x2,y2,xm-.5,ym,x4,y4,-.5,-.5, &
                                   a,ax,ay,axx,ayy,axy &
                                  )
                 dl=min(dp(i,j),max(0.,pbu(i,j)-pup(i,j)))
                 fd=a*dl+ax*dx(i,j)+ay*dy(i,j)
                 fdu_=fdu_+fd
                 qx=ax*dl+axx*dx(i,j)+axy*dy(i,j)
                 qy=ay*dl+axy*dx(i,j)+ayy*dy(i,j)
                 ftu_=ftu_+fd*td(i,j) &
                               +qx*tx(i,j)+qy*ty(i,j)
                 fsu_=fsu_+fd*sd(i,j) &
                               +qx*sx(i,j)+qy*sy(i,j)
              endif
              ! --- --- u-component of mass, heat and salt flux.
              uflx(i,j)=fdu_
              utflx(i,j)=ftu_
              usflx(i,j)=fsu_
           endif
        enddo
        !
     enddo
     delta=delta+(wallclock()-t1)
  enddo
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta
#endif

contains
  subroutine indxi(ipt,if,il,is)
    implicit none

    integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ipt
    integer, dimension (1-nbdy:jdm+nbdy,ms) :: if,il
    integer, dimension (1-nbdy:jdm+nbdy) :: is
    ! --- input array ipt contains 1 at grid point locations, 0 elsewhere
    ! --- output is arrays if, il, is  where
    ! --- if(j,k) gives row index of first point in column j for k-th section
    ! --- il(j,k) gives row index of last point
    ! --- is(j) gives number of sections in column j (maximum: ms)
    integer i,j,k,last
    do j=1-nbdy,jj+nbdy
       is(j) = 0
       do k=1,ms
          if(j,k) = 0
          il(j,k) = 0
       enddo
       k=1
       last = ipt(1-nbdy,j)
       if     (last .eq. 1) then
          if(j,k) = 1-nbdy
       endif
       do i=2-nbdy,ii+nbdy
          if      (last .eq. 1 .and. ipt(i,j) .eq. 0) then
             il(j,k) = i-1
             k = k+1
          elseif (last .eq. 0 .and. ipt(i,j) .eq. 1) then
             if     (k .gt. ms) then
                write(*,'(a,i5)')  'indxi problem on proc ', 0
                write(*,'(a,2i5)') ' error in indxi -- ms too small at i,j =',i0+i,j0+j
                stop '(indxi)'
             endif
             if(j,k) = i
          endif
          last = ipt(i,j)
       enddo
       if     (last .eq. 1) then
          il(j,k) = ii+nbdy
          is(j) = k
       else
          is(j) = k-1
       endif
    enddo
    return
  end subroutine indxi

  subroutine triint(ac,x1,y1,x2,y2,x3,y3,a,ax,ay,axx,ayy,axy &
#if defined(TRC) && defined(ATRC)
                     ,axxx,ayyy,axxy,axyy &
#endif
                     )
    !
    implicit none
    !
    real ac,x1,y1,x2,y2,x3,y3,a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
    real axxx,ayyy,axxy,axyy
#endif
    !
    real r1_3,r1_6,r1_12
    parameter (r1_3=1./3.,r1_6=1./6.,r1_12=1./12.)
#if defined(TRC) && defined(ATRC)
    real r1_10,r1_30
    parameter (r1_10=.1,r1_30=1./30.)
#endif
    !
    real xx,yy,xy1,xy2,xy3,xy
    !
    xx=x1*x2+x2*x3+x1*x3
    yy=y1*y2+y2*y3+y1*y3
    xy1=x1*y1
    xy2=x2*y2
    xy3=x3*y3
    xy=xy1+xy2+xy3
    !
    a=.5*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))*ac
    !
    ax=r1_3*(x1+x2+x3)
    ay=r1_3*(y1+y2+y3)
    axx=r1_6*(9.*ax*ax-xx)
    ayy=r1_6*(9.*ay*ay-yy)
    axy=r1_12*(9.*ax*ay+xy)
#if defined(TRC) && defined(ATRC)
    axxx=r1_10*((18.*axx-3.*xx)*ax+x1*x2*x3)
    ayyy=r1_10*((18.*ayy-3.*yy)*ay+y1*y2*y3)
    axxy=r1_30*(18.*axx*ay+3.*ax*xy+x1*xy1+x2*xy2+x3*xy3)
    axyy=r1_30*(18.*ayy*ax+3.*ay*xy+y1*xy1+y2*xy2+y3*xy3)
#endif
    !
    ax=ax*a
    ay=ay*a
    axx=axx*a
    ayy=ayy*a
    axy=axy*a
#if defined(TRC) && defined(ATRC)
    axxx=axxx*a
    ayyy=ayyy*a
    axxy=axxy*a
    axyy=axyy*a
#endif
    !
  end subroutine triint
  !
  subroutine penint(ac,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5, &
                      a,ax,ay,axx,ayy,axy &
#if defined(TRC) && defined(ATRC) 
                     ,axxx,ayyy,axxy,axyy &
#endif
                     )
    !
    implicit none
    !
    real ac,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,a,ax,ay,axx,ayy,axy
#if defined(TRC) && defined(ATRC)
    real axxx,ayyy,axxy,axyy
#endif
    !
    real r1_3,r1_6,r1_12
    parameter (r1_3=1./3.,r1_6=1./6.,r1_12=1./12.)
#if defined(TRC) && defined(ATRC)
    real r1_10,r1_30
    parameter (r1_10=.1,r1_30=1./30.)
#endif
    !
    real xx123,yy123,xx135,yy135,xx345,yy345,xy1,xy2,xy3,xy4,xy5, &
         xy123,xy135,xy345,a123,a135,a345,ax123,ax135,ax345, &
         ay123,ay135,ay345,axx123,axx135,axx345,ayy123,ayy135,ayy345, &
         axy123,axy135,axy345
#if defined(TRC) && defined(ATRC)
    real axxx123,axxx135,axxx345,ayyy123,ayyy135,ayyy345,&
         axxy123,axxy135,axxy345,axyy123,axyy135,axyy345
#endif
    !
    xx123=x1*x2+x2*x3+x1*x3
    yy123=y1*y2+y2*y3+y1*y3
    xx135=x1*x3+x3*x5+x1*x5
    yy135=y1*y3+y3*y5+y1*y5
    xx345=x3*x4+x4*x5+x3*x5
    yy345=y3*y4+y4*y5+y3*y5
    xy1=x1*y1
    xy2=x2*y2
    xy3=x3*y3
    xy4=x4*y4
    xy5=x5*y5
    xy123=xy1+xy2+xy3
    xy135=xy1+xy3+xy5
    xy345=xy3+xy4+xy5
    !
    a123=.5*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))*ac
    a135=.5*((x3-x1)*(y5-y1)-(y3-y1)*(x5-x1))*ac
    a345=.5*((x4-x3)*(y5-y3)-(y4-y3)*(x5-x3))*ac
    !
    ax123=r1_3*(x1+x2+x3)
    ay123=r1_3*(y1+y2+y3)
    ax135=r1_3*(x1+x3+x5)
    ay135=r1_3*(y1+y3+y5)
    ax345=r1_3*(x3+x4+x5)
    ay345=r1_3*(y3+y4+y5)
    axx123=r1_6*(9.*ax123*ax123-xx123)
    ayy123=r1_6*(9.*ay123*ay123-yy123)
    axy123=r1_12*(9.*ax123*ay123+xy123)
    axx135=r1_6*(9.*ax135*ax135-xx135)
    ayy135=r1_6*(9.*ay135*ay135-yy135)
    axy135=r1_12*(9.*ax135*ay135+xy135)
    axx345=r1_6*(9.*ax345*ax345-xx345)
    ayy345=r1_6*(9.*ay345*ay345-yy345)
    axy345=r1_12*(9.*ax345*ay345+xy345)
#if defined(TRC) && defined(ATRC)
    axxx123=r1_10*((18.*axx123-3.*xx123)*ax123+x1*x2*x3)
    ayyy123=r1_10*((18.*ayy123-3.*yy123)*ay123+y1*y2*y3)
    axxy123=r1_30*(18.*axx123*ay123+3.*ax123*xy123 &
                  +x1*xy1+x2*xy2+x3*xy3)
    axyy123=r1_30*(18.*ayy123*ax123+3.*ay123*xy123 &
                  +y1*xy1+y2*xy2+y3*xy3)
    axxx135=r1_10*((18.*axx135-3.*xx135)*ax135+x1*x3*x5)
    ayyy135=r1_10*((18.*ayy135-3.*yy135)*ay135+y1*y3*y5)
    axxy135=r1_30*(18.*axx135*ay135+3.*ax135*xy135 &
                  +x1*xy1+x3*xy3+x5*xy5)
    axyy135=r1_30*(18.*ayy135*ax135+3.*ay135*xy135 &
                  +y1*xy1+y3*xy3+y5*xy5)
    axxx345=r1_10*((18.*axx345-3.*xx345)*ax345+x3*x4*x5)
    ayyy345=r1_10*((18.*ayy345-3.*yy345)*ay345+y3*y4*y5)
    axxy345=r1_30*(18.*axx345*ay345+3.*ax345*xy345 &
                  +x3*xy3+x4*xy4+x5*xy5)
    axyy345=r1_30*(18.*ayy345*ax345+3.*ay345*xy345 &
                  +y3*xy3+y4*xy4+y5*xy5)
#endif
    !
    a=a123+a135+a345
    ax=ax123*a123+ax135*a135+ax345*a345
    ay=ay123*a123+ay135*a135+ay345*a345
    axx=axx123*a123+axx135*a135+axx345*a345
    ayy=ayy123*a123+ayy135*a135+ayy345*a345
    axy=axy123*a123+axy135*a135+axy345*a345
#if defined(TRC) && defined(ATRC)
    axxx=axxx123*a123+axxx135*a135+axxx345*a345
    ayyy=ayyy123*a123+ayyy135*a135+ayyy345*a345
    axxy=axxy123*a123+axxy135*a135+axxy345*a345
    axyy=axyy123*a123+axyy135*a135+axyy345*a345
#endif
    !
  end subroutine penint
  
end program remap_vel_u
