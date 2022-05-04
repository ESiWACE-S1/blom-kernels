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
  integer :: isp(1-nbdy:jdm+nbdy-1)
  integer :: ifp(-1:jdm+2, ms) ! first section index per section
  integer :: ilp(-1:jdm+2, ms) ! last section index per section
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: dp, pup, plo
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fdu, fdv, ftu, ftv, fsu, fsv, cu, cv
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbmin
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: dx, dy, xd, yd, tx, ty, temp, td, sx, sy, saln, sd
  real :: dpeps
  real :: dxi, dyi, dpsw, dps, dpc, dpe, dpnw, dpn, dpne, dpse, dpw, dgmx, dfmx, dfmn, q, q1, q2, q3, q4, tgmx, tgmn, tfmx, tfmn
  real :: sgmx, sgmn, sfmx, sfmn
  integer :: mrg
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

  i0 = 0
  ii = idm
  j0 = 0
  jj = jdm
  kk = kdm
  mrg = 0
  dpeps = 1.0

  write(*,*) "COMPUTE: remap_zero"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "DEFAULT: iterating over",MAX_ITERATIONS, " iterations..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do j=-1,jj+2
     isp(j) = 1
     ifp(j,1) = -1
     ilp(j,1) = ii+2
  end do

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

end program remap_zero
