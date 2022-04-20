program advect_pbmin
  implicit none
  integer, parameter :: nbdy = 3
  integer, parameter :: idm = 64
  integer, parameter :: jdm = 128
  integer, parameter :: kdm = 53

  integer :: l,i,j,iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne
  integer :: ii,jj,kk, i0,j0
  integer, parameter :: ms = 10
  integer :: isp(ms)
  integer :: ifp(-1:jdm+2, ms) ! first section index per section
  integer :: ilp(-1:jdm+2, ms) ! last section index per section
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbmin
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy, kdm+1) :: p

  ! masks
  integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ip, iu, iv, iq
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: util1, util2, util3
  real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: depth
  real :: depmax
  ! bigrid
  logical :: lperiodi, lperiodj, larctic
  integer :: nreg

  i0 = 0
  ii = idm
  j0 = 0
  jj = jdm
  kk = kdm

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
  do j=1-nbdy,jdm+nbdy
     write(*,'(2i3,a$)') j, 1, 'depth: '
     do i=1-nbdy,idm+nbdy
        write(*,'(1(i1,a)$)') nint(depth(i,j)), "|"
     enddo
     write(*,'(x)')
  enddo
  do j=1-nbdy,jdm+nbdy
     write(*,'(2i3,a$)') j, 1, ': '
     do i=1-nbdy,idm+nbdy
        write(*,'(4(i1,a)$)') ip(i,j), ",", iq(i,j), ",", iu(i,j) &
             , ",", iv(i,j), "|"
     enddo
     write(*,'(x)')
  enddo
  stop '(DEBUG)'
  !$OMP PARALLEL DO PRIVATE( &
  !$OMP&  l,i,iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne)
  do j=-1,jj+2
     do l=1,isp(j)
        do i=max(-1,ifp(j,l)),min(ii+2,ilp(j,l))
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
           pbmin(i,j)= &
                min(p(isw,jsw,kk+1),p(i  ,js ,kk+1),p(ise,jse,kk+1), &
                p(iw ,j  ,kk+1),p(i  ,j  ,kk+1),p(ie ,j  ,kk+1), &
                p(inw,jnw,kk+1),p(i  ,jn ,kk+1),p(ine,jne,kk+1))
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end program advect_pbmin
