program remap_zero
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
  integer :: isp(1-nbdy:jdm+nbdy)
  integer :: ifp(1-nbdy:jdm+nbdy, ms) ! first section index per section
  integer :: ilp(1-nbdy:jdm+nbdy, ms) ! last section index per section
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: dp, pup, plo
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fdu, fdv, ftu, ftv, fsu, fsv, cu, cv
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbmin
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: dx, dy, xd, yd, tx, ty, temp, td, sx, sy, saln, sd
  real :: dpeps
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

  real :: pbmin_, xd_, yd_, sx_, sy_, sd_, tx_, ty_, td_, dx_, dy_, dp_
  real :: dp_min, dp_max, t_min, t_max, s_min, s_max
  real :: te, tw, tn, ts, se, sw, sn, ss

  i0 = 0
  ii = idm
  j0 = 0
  jj = jdm
  kk = kdm
  mrg = 0
  dpeps = 1.0
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
  !do j=1-nbdy,jdm+nbdy
  !   write(*,'(2i3,a$)') j, 1, 'depth: '
  !   do i=1-nbdy,idm+nbdy
  !      write(*,'(1(i1,a)$)') nint(depth(i,j)), "|"
  !   enddo
  !   write(*,'(x)')
  !enddo
  !do j=1-nbdy,jdm+nbdy
  !   write(*,'(2i3,a$)') j, 1, ': '
  !   do i=1-nbdy,idm+nbdy
  !      write(*,'(4(i1,a)$)') ip(i,j), ",", iq(i,j), ",", iu(i,j) &
  !           , ",", iv(i,j), "|"
  !   enddo
  !   write(*,'(x)')
  !enddo
  pbmin = 1.0

  write(*,*) "COMPUTE: remap_limited_grad"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "DEFAULT: iterating over",MAX_ITERATIONS, " iterations..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call indxi(ip,ifp,ilp,isp)

  delta=0.0
  do n_it=1, MAX_ITERATIONS
     t1=wallclock()
     ! --- ------------------------------------------------------------------
     ! --- Compute limited gradients, center of mass coordinates, and
     ! --- non-dimensional velocities.
     ! --- ------------------------------------------------------------------
     !
     do j=1-mrg-1,jj+mrg+1
        do l=1,isp(j)
           do i=max(1-mrg-1,ifp(j,l)),min(ii+mrg+1,ilp(j,l))
              ! --- --- Define indices for grid cell neighbors, ensuring that only wet
              ! --- --- points are used.
              iw=i-iu(i  ,j)
              ie=i+iu(i+1,j)
              js=j-iv(i,j  )
              jn=j+iv(i,j+1)
              isw=i*(1-ip(iw,js))+iw*ip(iw,js)
              jsw=j*(1-ip(iw,js))+js*ip(iw,js)
              ise=i*(1-ip(ie,js))+ie*ip(ie,js)
              jse=j*(1-ip(ie,js))+js*ip(ie,js)
              inw=i*(1-ip(iw,jn))+iw*ip(iw,jn)
              jnw=j*(1-ip(iw,jn))+jn*ip(iw,jn)
              ine=i*(1-ip(ie,jn))+ie*ip(ie,jn)
              jne=j*(1-ip(ie,jn))+jn*ip(ie,jn)

              dxi=1./max(1,ie-iw)
              dyi=1./max(1,jn-js)
              ! --- --- Compute limited gradient for layer pressure thickness and
              ! --- --- center of mass coordinate.
              dpsw=max(dpeps,min(pbmin(i,j)-pup(isw,jsw),dp(isw,jsw)))
              dps =max(dpeps,min(pbmin(i,j)-pup(i  ,js ),dp(i  ,js )))
              dpse=max(dpeps,min(pbmin(i,j)-pup(ise,jse),dp(ise,jse)))
              dpw =max(dpeps,min(pbmin(i,j)-pup(iw ,j  ),dp(iw ,j  )))
              dpc =max(dpeps,min(pbmin(i,j)-pup(i  ,j  ),dp(i  ,j  )))
              dpe =max(dpeps,min(pbmin(i,j)-pup(ie ,j  ),dp(ie ,j  )))
              dpnw=max(dpeps,min(pbmin(i,j)-pup(inw,jnw),dp(inw,jnw)))
              dpn =max(dpeps,min(pbmin(i,j)-pup(i  ,jn ),dp(i  ,jn )))
              dpne=max(dpeps,min(pbmin(i,j)-pup(ine,jne),dp(ine,jne)))
              dx(i,j)=(dpe-dpw)*dxi
              dy(i,j)=(dpn-dps)*dyi
              dgmx=.5*(abs(dx(i,j))+abs(dy(i,j)))
              dfmx=max(0.,max(dpsw,dps,dpse,dpw,dpe,dpnw,dpn,dpne)-dpc)
              dfmn=min(0.,min(dpsw,dps,dpse,dpw,dpe,dpnw,dpn,dpne)-dpc)
              if (dfmx.gt.0..and.dfmn.lt.0.) then
                 q=min(dfmx/max(dfmx,dgmx),dfmn/min(dfmn,-dgmx))
                 dx(i,j)=dx(i,j)*q
                 dy(i,j)=dy(i,j)*q
                 xd(i,j)=dx(i,j)/(12.*dp(i,j))
                 yd(i,j)=dy(i,j)/(12.*dp(i,j))
              else
                 dx(i,j)=0.
                 dy(i,j)=0.
                 xd(i,j)=0.
                 yd(i,j)=0.
              endif
              ! --- --- Compute limited gradients for temperature, salinity, and
              ! --- --- density
              tx(i,j)=(temp(ie,j)-temp(iw,j))*dxi
              ty(i,j)=(temp(i,jn)-temp(i,js))*dyi
              q1=tx(i,j)*(-.5-xd(i,j))
              q2=tx(i,j)*( .5-xd(i,j))
              q3=ty(i,j)*(-.5-yd(i,j))
              q4=ty(i,j)*( .5-yd(i,j))
              tgmx=max(q1,q2)+max(q3,q4)
              tgmn=min(q1,q2)+min(q3,q4)
              tfmx=max(0.,max(temp(isw,jsw),temp(i  ,js ),      &
                              temp(ise,jse),temp(iw ,j  ), &
                              temp(ie ,j  ),temp(inw,jnw), &
                              temp(i  ,jn ),temp(ine,jne)) &
                          -temp(i,j))
              tfmn=min(0.,min(temp(isw,jsw),temp(i  ,js ),      &
                              temp(ise,jse),temp(iw ,j  ), &
                              temp(ie ,j  ),temp(inw,jnw), &
                              temp(i  ,jn ),temp(ine,jne)) &
                          -temp(i,j))
              if (tfmx.gt.0..and.tfmn.lt.0.) then
                 q=min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
                 tx(i,j)=tx(i,j)*q
                 ty(i,j)=ty(i,j)*q
                 td(i,j)=temp(i,j)-tx(i,j)*xd(i,j)-ty(i,j)*yd(i,j)
              else
                 tx(i,j)=0.
                 ty(i,j)=0.
                 td(i,j)=temp(i,j)
              endif

              sx(i,j)=(saln(ie,j)-saln(iw,j))*dxi
              sy(i,j)=(saln(i,jn)-saln(i,js))*dyi
              q1=sx(i,j)*(-.5-xd(i,j))
              q2=sx(i,j)*( .5-xd(i,j))
              q3=sy(i,j)*(-.5-yd(i,j))
              q4=sy(i,j)*( .5-yd(i,j))
              sgmx=max(q1,q2)+max(q3,q4)
              sgmn=min(q1,q2)+min(q3,q4)
              sfmx=max(0.,max(saln(isw,jsw),saln(i  ,js ),     &
                                  saln(ise,jse),saln(iw ,j  ), &
                                  saln(ie ,j  ),saln(inw,jnw), &
                                  saln(i  ,jn ),saln(ine,jne)) &
                             -saln(i,j))
              sfmn=min(0.,min(saln(isw,jsw),saln(i  ,js ),     &
                                  saln(ise,jse),saln(iw ,j  ), &
                                  saln(ie ,j  ),saln(inw,jnw), &
                                  saln(i  ,jn ),saln(ine,jne)) &
                             -saln(i,j))
              if (sfmx.gt.0..and.sfmn.lt.0.) then
                 q=min(sfmx/max(sfmx,sgmx),sfmn/min(sfmn,sgmn))
                 sx(i,j)=sx(i,j)*q
                 sy(i,j)=sy(i,j)*q
                 sd(i,j)=saln(i,j)-sx(i,j)*xd(i,j)-sy(i,j)*yd(i,j)
              else
                 sx(i,j)=0.
                 sy(i,j)=0.
                 sd(i,j)=saln(i,j)
              endif
           end do
        end do
     end do
     delta=delta+(wallclock()-t1)
  enddo
  delta_orig=delta
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

  i = random_uniform(1,idm)
  j = random_uniform(1,jdm)

#define REMAP_LIMITED_GRAD_OPT1
#ifdef REMAP_LIMITED_GRAD_OPT1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "opt1: iterating over",MAX_ITERATIONS, " iterations..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  delta=0.0
  do n_it=1, MAX_ITERATIONS
     t1=wallclock()
     ! --- ------------------------------------------------------------------
     ! --- Compute limited gradients, center of mass coordinates, and
     ! --- non-dimensional velocities.
     ! --- ------------------------------------------------------------------
     !
     do j=1-mrg-1,jj+mrg+1
        do i=1-mrg-1,ii+mrg+1
           if(ip(i,j) /= 0) then
              ! --- --- Define indices for grid cell neighbors, ensuring that only wet
              ! --- --- points are used.
              iw=i-iu(i  ,j)
              ie=i+iu(i+1,j)
              dxi=1./max(1,ie-iw)

              js=j-iv(i,j  )
              jn=j+iv(i,j+1)
              dyi=1./max(1,jn-js)

              isw=i*(1-ip(iw,js))+iw*ip(iw,js)
              jsw=j*(1-ip(iw,js))+js*ip(iw,js)
              ise=i*(1-ip(ie,js))+ie*ip(ie,js)
              jse=j*(1-ip(ie,js))+js*ip(ie,js)
              inw=i*(1-ip(iw,jn))+iw*ip(iw,jn)
              jnw=j*(1-ip(iw,jn))+jn*ip(iw,jn)
              ine=i*(1-ip(ie,jn))+ie*ip(ie,jn)
              jne=j*(1-ip(ie,jn))+jn*ip(ie,jn)

              ! --- --- Compute limited gradient for layer pressure thickness and
              ! --- --- center of mass coordinate.
              pbmin_ = pbmin(i,j)
              dpsw=max(dpeps,min(pbmin_-pup(isw,jsw),dp(isw,jsw)))
              dps =max(dpeps,min(pbmin_-pup(i  ,js ),dp(i  ,js )))
              dpse=max(dpeps,min(pbmin_-pup(ise,jse),dp(ise,jse)))
              dpw =max(dpeps,min(pbmin_-pup(iw ,j  ),dp(iw ,j  )))
              dpc =max(dpeps,min(pbmin_-pup(i  ,j  ),dp(i  ,j  )))
              dpe =max(dpeps,min(pbmin_-pup(ie ,j  ),dp(ie ,j  )))
              dpnw=max(dpeps,min(pbmin_-pup(inw,jnw),dp(inw,jnw)))
              dpn =max(dpeps,min(pbmin_-pup(i  ,jn ),dp(i  ,jn )))
              dpne=max(dpeps,min(pbmin_-pup(ine,jne),dp(ine,jne)))
              dx_=(dpe-dpw)*dxi
              dy_=(dpn-dps)*dyi
              ! q1=dx(i,j)*(-.5-xd(i,j))
              ! q2=dx(i,j)*( .5-xd(i,j))
              ! q3=dy(i,j)*(-.5-yd(i,j))
              ! q4=dy(i,j)*( .5-yd(i,j))
              ! dgmx=max(q1,q2)+max(q3,q4)
              ! dgmn=min(q1,q2)+min(q3,q4)
              dfmx=max(0.,max(dpsw,dps,dpse,dpw,dpe,dpnw,dpn,dpne)-dpc)
              dfmn=min(0.,min(dpsw,dps,dpse,dpw,dpe,dpnw,dpn,dpne)-dpc)
              if (dfmx.gt.0..and.dfmn.lt.0.) then
                 dgmx=.5*(abs(dx_)+abs(dy_))
                 q=min(dfmx/max(dfmx,dgmx),dfmn/min(dfmn,-dgmx))
                 dx_=dx_*q
                 dy_=dy_*q
                 xd_=dx_/(12.*dp(i,j))
                 yd_=dy_/(12.*dp(i,j))
              else
                 dx(i,j)=0.
                 dy(i,j)=0.
                 xd_=0.
                 yd_=0.
              endif
              dx(i,j) = dx_
              dy(i,j) = dy_
              ! --- --- Compute limited gradients for temperature, salinity, and
              ! --- --- density
              tx_=(temp(ie,j)-temp(iw,j))*dxi
              ty_=(temp(i,jn)-temp(i,js))*dyi
              td_=temp(i,j)
              tfmx=max(0.,max(temp(isw,jsw),temp(i  ,js ),      &
                              temp(ise,jse),temp(iw ,j  ), &
                              temp(ie ,j  ),temp(inw,jnw), &
                              temp(i  ,jn ),temp(ine,jne)) &
                          -td_)
              tfmn=min(0.,min(temp(isw,jsw),temp(i  ,js ),      &
                              temp(ise,jse),temp(iw ,j  ), &
                              temp(ie ,j  ),temp(inw,jnw), &
                              temp(i  ,jn ),temp(ine,jne)) &
                          -td_)
              if (tfmx.gt.0..and.tfmn.lt.0.) then
                 q1=tx_*(-.5-xd_)
                 q2=tx_*( .5-xd_)
                 q3=ty_*(-.5-yd_)
                 q4=ty_*( .5-yd_)
                 tgmx=max(q1,q2)+max(q3,q4)
                 tgmn=min(q1,q2)+min(q3,q4)
                 q=min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
                 tx_=tx_*q
                 ty_=ty_*q
                 td_=td_-tx_*xd_-ty_*yd_
              else
                 tx_=0.
                 ty_=0.
              endif
              tx(i,j) = tx_
              ty(i,j) = ty_
              td(i,j) = td_

              sx_=(saln(ie,j)-saln(iw,j))*dxi
              sy_=(saln(i,jn)-saln(i,js))*dyi
              sd_=saln(i,j)
              sfmx=max(0.,max(saln(isw,jsw),saln(i  ,js ),     &
                                  saln(ise,jse),saln(iw ,j  ), &
                                  saln(ie ,j  ),saln(inw,jnw), &
                                  saln(i  ,jn ),saln(ine,jne)) &
                             -sd_)
              sfmn=min(0.,min(saln(isw,jsw),saln(i  ,js ),     &
                                  saln(ise,jse),saln(iw ,j  ), &
                                  saln(ie ,j  ),saln(inw,jnw), &
                                  saln(i  ,jn ),saln(ine,jne)) &
                             -sd_)
              if (sfmx.gt.0..and.sfmn.lt.0.) then
                 q1=sx_*(-.5-xd_)
                 q2=sx_*( .5-xd_)
                 q3=sy_*(-.5-yd_)
                 q4=sy_*( .5-yd_)
                 sgmx=max(q1,q2)+max(q3,q4)
                 sgmn=min(q1,q2)+min(q3,q4)
                 q=min(sfmx/max(sfmx,sgmx),sfmn/min(sfmn,sgmn))
                 sx_=sx_*q
                 sy_=sy_*q
                 sd_=sd_-sx_*xd_-sy_*yd_
              else
                 sx_=0.
                 sy_=0.
              endif
              sx(i,j) = sx_
              sy(i,j) = sy_
              sd(i,j) = sd_
           end if
        end do
     end do
     delta=delta+(wallclock()-t1)
  enddo
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta
#endif

#define REMAP_LIMITED_GRAD_OPT2
#ifdef REMAP_LIMITED_GRAD_OPT2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "opt2: iterating over",MAX_ITERATIONS, " iterations..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  delta=0.0
  do n_it=1, MAX_ITERATIONS
     t1=wallclock()
     ! --- ------------------------------------------------------------------
     ! --- Compute limited gradients, center of mass coordinates, and
     ! --- non-dimensional velocities.
     ! --- ------------------------------------------------------------------
     !
     do j=1-mrg-1,jj+mrg+1
        do i=1-mrg-1,ii+mrg+1
           if(ip(i,j) /= 0) then
              ! --- --- Define indices for grid cell neighbors, ensuring that only wet
              ! --- --- points are used.
              iw=i-iu(i  ,j)
              ie=i+iu(i+1,j)
              dxi=1./max(1,ie-iw)

              js=j-iv(i,j  )
              jn=j+iv(i,j+1)
              dyi=1./max(1,jn-js)

              pbmin_ = pbmin(i,j)

              dp_ = max(dpeps,min(pbmin_-pup(i  ,j  ),dp(j  ,j  )))
              dp_max = dp_
              dp_min = dp_
              dpw = dp_
              dpe = dp_
              dps = dp_
              dpn = dp_

              td_ = temp(i  , j  )
              t_max = td_
              t_min = td_
              tw = td_
              te = td_
              ts = td_
              tn = td_

              sd_ = saln(i  , j  )
              s_max = sd_
              s_min = sd_
              sw = sd_
              se = sd_
              ss = sd_
              sn = sd_

              ! iw defined
              if(iw /= 0) then
                 if(js /= 0 .and. ip(i-1,j-1) /= 0) then
                    dp_ = max(dpeps,min(pbmin_-pup(i-1,j-1),dp(i-1,j-1)))
                    dp_max = max(dp_max, dp_)
                    dp_min = min(dp_min, dp_)

                    td_ = temp(i-1, j-1)
                    t_max = max(t_max, td_)
                    t_min = max(t_min, td_)

                    sd_ = saln(i-1, j-1)
                    s_max = max(s_max, sd_)
                    s_min = max(s_min, sd_)
                 end if
                 if(              ip(i-1,j  ) /= 0) then
                    dp_ = max(dpeps,min(pbmin_-pup(i-1,j  ),dp(i-1,j  )))
                    dp_max = max(dp_max, dp_)
                    dp_min = min(dp_min, dp_)

                    td_ = temp(i-1, j  )
                    t_max = max(t_max, td_)
                    t_min = max(t_min, td_)

                    sd_ = saln(i-1, j  )
                    s_max = max(s_max, sd_)
                    s_min = max(s_min, sd_)
                    ! dpw
                    dpw = dp_
                    tw = td_
                    sw = sd_
                 end if
                 if(jn /= 0 .and. ip(i-1,j+1) /= 0) then
                    dp_ = max(dpeps,min(pbmin_-pup(i-1,j+1),dp(i-1,j+1)))
                    dp_max = max(dp_max, dp_)
                    dp_min = min(dp_min, dp_)

                    td_ = temp(i-1, j+1)
                    t_max = max(t_max, td_)
                    t_min = max(t_min, td_)

                    sd_ = saln(i-1, j+1)
                    s_max = max(s_max, sd_)
                    s_min = max(s_min, sd_)
                 end if
              end if

                 if(js /= 0 .and. ip(i  ,j-1) /= 0) then
                    dp_ = max(dpeps,min(pbmin_-pup(i  ,j-1),dp(j  ,j-1)))
                    dp_max = max(dp_max, dp_)
                    dp_min = min(dp_min, dp_)

                    td_ = temp(i-1, j-1)
                    t_max = max(t_max, td_)
                    t_min = max(t_min, td_)

                    sd_ = saln(i-1, j-1)
                    s_max = max(s_max, sd_)
                    s_min = max(s_min, sd_)
                    ! dps
                    dps = dp_
                    ts = td_
                    ss = sd_
                 end if

                 if(jn /= 0 .and. ip(i  ,j+1) /= 0) then
                    dp_ = max(dpeps,min(pbmin_-pup(i  ,j+1),dp(j  ,j+1)))
                    dp_max = max(dp_max, dp_)
                    dp_min = min(dp_min, dp_)

                    td_ = temp(i  , j+1)
                    t_max = max(t_max, td_)
                    t_min = max(t_min, td_)

                    sd_ = saln(i  , j+1)
                    s_max = max(s_max, sd_)
                    s_min = max(s_min, sd_)
                    ! dpn
                    dpn = dp_
                    tn = td_
                    sn = sd_
                 end if

               ! ie defined
              if(ie /= 0) then
                 if(js /= 0 .and. ip(i+1,j-1) /= 0) then
                    dp_ = max(dpeps,min(pbmin_-pup(i+1,j-1),dp(i+1,j-1)))
                    dp_max = max(dp_max, dp_)
                    dp_min = min(dp_min, dp_)

                    td_ = temp(i+1, j-1)
                    t_max = max(t_max, td_)
                    t_min = max(t_min, td_)

                    sd_ = saln(i+1, j-1)
                    s_max = max(s_max, sd_)
                    s_min = max(s_min, sd_)
                 end if
                 if(              ip(i+1,j  ) /= 0) then
                    dp_ = max(dpeps,min(pbmin_-pup(i+1,j  ),dp(i+1,j  )))
                    dp_max = max(dp_max, dp_)
                    dp_min = min(dp_min, dp_)

                    td_ = temp(i+1, j  )
                    t_max = max(t_max, td_)
                    t_min = max(t_min, td_)

                    sd_ = saln(i+1, j  )
                    s_max = max(s_max, sd_)
                    s_min = max(s_min, sd_)
                    ! dpe
                    dpe = dp_
                    te = td_
                    se = sd_
                 end if
                 if(jn /= 0 .and. ip(i+1,j+1) /= 0) then
                    dp_ = max(dpeps,min(pbmin_-pup(i+1,j-1),dp(i+1,j+1)))
                    dp_max = max(dp_max, dp_)
                    dp_min = min(dp_min, dp_)

                    td_ = temp(i+1, j+1)
                    t_max = max(t_max, td_)
                    t_min = max(t_min, td_)

                    sd_ = saln(i+1, j+1)
                    s_max = max(s_max, sd_)
                    s_min = max(s_min, sd_)
                 end if
              end if

              ! --- --- center of mass coordinate.
              dx_=(dpe-dpw)*dxi
              dy_=(dpn-dps)*dyi
              dfmx=max(0.,dp_max-dpc)
              dfmn=min(0.,dp_min-dpc)
              if (dfmx.gt.0..and.dfmn.lt.0.) then
                 ! q1=dx_*(-.5-xd_)
                 ! q2=dx_*( .5-xd_)
                 ! q3=dy_*(-.5-yd_)
                 ! q4=dy_*( .5-yd_)
                 ! dgmx=max(q1,q2)+max(q3,q4)
                 ! dgmn=min(q1,q2)+min(q3,q4)
                 dgmx=.5*(abs(dx_)+abs(dy_))
                 q=min(dfmx/max(dfmx,dgmx),dfmn/min(dfmn,-dgmx))
                 dx_=dx_*q
                 dy_=dy_*q
                 xd_=dx_/(12.*dp(i,j))
                 yd_=dy_/(12.*dp(i,j))
              else
                 dx(i,j)=0.
                 dy(i,j)=0.
                 xd_=0.
                 yd_=0.
              endif
              dx(i,j) = dx_
              dy(i,j) = dy_
              ! --- --- Compute limited gradients for temperature, salinity, and
              ! --- --- density
              tx_=(te-tw)*dxi
              ty_=(tn-ts)*dyi
              td_=temp(i,j)
              tfmx=max(0.,t_max-td_)
              tfmn=min(0.,t_min-td_)
              if (tfmx.gt.0..and.tfmn.lt.0.) then
                 q1=tx_*(-.5-xd_)
                 q2=tx_*( .5-xd_)
                 q3=ty_*(-.5-yd_)
                 q4=ty_*( .5-yd_)
                 tgmx=max(q1,q2)+max(q3,q4)
                 tgmn=min(q1,q2)+min(q3,q4)
                 q=min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
                 tx_=tx_*q
                 ty_=ty_*q
                 td_=td_-tx_*xd_-ty_*yd_
              else
                 tx_=0.
                 ty_=0.
              endif
              tx(i,j) = tx_
              ty(i,j) = ty_
              td(i,j) = td_

              sx_=(se-sw)*dxi
              sy_=(sn-ss)*dyi
              sd_=saln(i,j)
              sfmx=max(0.,s_max-sd_)
              sfmn=min(0.,s_min-sd_)
              if (sfmx.gt.0..and.sfmn.lt.0.) then
                 q1=sx_*(-.5-xd_)
                 q2=sx_*( .5-xd_)
                 q3=sy_*(-.5-yd_)
                 q4=sy_*( .5-yd_)
                 sgmx=max(q1,q2)+max(q3,q4)
                 sgmn=min(q1,q2)+min(q3,q4)
                 q=min(sfmx/max(sfmx,sgmx),sfmn/min(sfmn,sgmn))
                 sx_=sx_*q
                 sy_=sy_*q
                 sd_=sd_-sx_*xd_-sy_*yd_
              else
                 sx_=0.
                 sy_=0.
              endif
              sx(i,j) = sx_
              sy(i,j) = sy_
              sd(i,j) = sd_
           end if
        end do
     end do
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


end program remap_zero
