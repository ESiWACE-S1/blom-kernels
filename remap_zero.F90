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
  real :: dpeps
  integer :: mrg

  integer, parameter :: MAX_ITERATIONS=100000
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
     do j=1-mrg-2,jj+mrg+2
        do l=1,isp(j)
           do i=max(1-mrg-2,ifp(j,l)),min(ii+mrg+2,ilp(j,l))
              dp(i,j)=max(0.,dp(i,j))+dpeps
              pup(i,j)=plo(i,j)-dp(i,j)
           enddo
        enddo
        do i=1-mrg-1,ii+mrg+1
           fdu(i,j)=0.
           fdv(i,j)=0.
           ftu(i,j)=0.
           ftv(i,j)=0.
           fsu(i,j)=0.
           fsv(i,j)=0.
#ifdef TRC
#  ifdef ATRC
           do nt=1,ntr-natr
#    if defined(TKE) && !defined(TKEADV)
              if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
              ftru(nt,i,j)=0.
              ftrv(nt,i,j)=0.
           enddo
           do nt=1,natr
              fagu(nt,i,j)=0.
              fagv(nt,i,j)=0.
           enddo
#  else
           do nt=1,ntr
#    if defined(TKE) && !defined(TKEADV)
              if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#    endif
              ftru(nt,i,j)=0.
              ftrv(nt,i,j)=0.
           enddo
#  endif
#endif
           cu(i,j)=0.
           cv(i,j)=0.
        enddo
     enddo
     delta=delta+(wallclock()-t1)
  enddo
  delta_orig=delta
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta

  i = random_uniform(1,idm)
  j = random_uniform(1,jdm)
  write(*,'(a,2i8,10e14.6)') 'out', i,j, dp(i,j), pup(i,j),cu(i,j), cv(i,j)&
       , fdu(i,j) &
       , fdv(i,j) &
       , ftu(i,j) &
       , ftv(i,j) &
       , fsu(i,j) &
       , fsv(i,j)

#define REMAP_RETVEIL_OPT1
#ifdef REMAP_RETVEIL_OPT1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "full loops: iterating over",MAX_ITERATIONS, " iterations..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  delta=0.0
  do n_it=1, MAX_ITERATIONS
     t1=wallclock()
     do j=1-mrg-2,jj+mrg+2
        do i=1-mrg-2,ii+mrg+2
           dp(i,j)=max(0.,dp(i,j))+dpeps
           pup(i,j)=plo(i,j)-dp(i,j)
        enddo
     enddo
     do j=1-mrg-2,jj+mrg+2
        do i=1-mrg-1,ii+mrg+1
           fdu(i,j)=0.
           fdv(i,j)=0.
           ftu(i,j)=0.
           ftv(i,j)=0.
           fsu(i,j)=0.
           fsv(i,j)=0.
           cu(i,j)=0.
           cv(i,j)=0.
        enddo
     enddo
     delta=delta+(wallclock()-t1)
  enddo
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta
#endif
     i = random_uniform(1,idm)
     j = random_uniform(1,jdm)
     write(*,'(a,2i8,10e14.6)') 'out', i,j, dp(i,j), pup(i,j),cu(i,j), cv(i,j)&
          , fdu(i,j) &
          , fdv(i,j) &
          , ftu(i,j) &
          , ftv(i,j) &
          , fsu(i,j) &
          , fsv(i,j)

#define REMAP_RETVEIL_OPT2
#ifdef REMAP_RETVEIL_OPT2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "full loops: iterating over",MAX_ITERATIONS, " iterations..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  delta=0.0
  do n_it=1, MAX_ITERATIONS
     t1=wallclock()
     do j=1-mrg-2,jj+mrg+2
        do i=1-mrg-2,ii+mrg+2
           dp(i,j)=max(0.,dp(i,j))+dpeps
           pup(i,j)=plo(i,j)-dp(i,j)
        enddo
     enddo
     do j=1-mrg-2,jj+mrg+2
        do i=1-mrg-1,ii+mrg+1
           fdu(i,j)=0.
        enddo
     enddo
     do j=1-mrg-2,jj+mrg+2
        do i=1-mrg-1,ii+mrg+1
           fdv(i,j)=0.
        enddo
     enddo
     do j=1-mrg-2,jj+mrg+2
        do i=1-mrg-1,ii+mrg+1
           ftu(i,j)=0.
        enddo
     enddo
     do j=1-mrg-2,jj+mrg+2
        do i=1-mrg-1,ii+mrg+1
           ftv(i,j)=0.
        enddo
     enddo
     do j=1-mrg-2,jj+mrg+2
        do i=1-mrg-1,ii+mrg+1
           fsu(i,j)=0.
        enddo
     enddo
     do j=1-mrg-2,jj+mrg+2
        do i=1-mrg-1,ii+mrg+1
           fsv(i,j)=0.
        enddo
     enddo
     do j=1-mrg-2,jj+mrg+2
        do i=1-mrg-1,ii+mrg+1
           cu(i,j)=0.
        enddo
     enddo
     do j=1-mrg-2,jj+mrg+2
        do i=1-mrg-1,ii+mrg+1
           cv(i,j)=0.
        enddo
     enddo
     delta=delta+(wallclock()-t1)
  enddo
  write(*,*) "done"
  write(*,'(a,3(f14.6,x))') "timing", delta, delta/real(MAX_ITERATIONS), delta_orig/delta
#endif
     i = random_uniform(1,idm)
     j = random_uniform(1,jdm)
     write(*,'(a,2i8,10e14.6)') 'out', i,j, dp(i,j), pup(i,j),cu(i,j), cv(i,j)&
          , fdu(i,j) &
          , fdv(i,j) &
          , ftu(i,j) &
          , ftv(i,j) &
          , fsu(i,j) &
          , fsv(i,j)



end program remap_zero
