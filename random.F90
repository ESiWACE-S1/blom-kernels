! https://masuday.github.io/fortran_tutorial/random.html
module random_mod
  implicit none
  private

  public :: init_seed
  
  public :: random_uniform
  interface random_uniform
     module procedure random_uniform_int
     module procedure random_uniform_real
  end interface random_uniform
contains
  subroutine init_seed()
    implicit none
    integer :: n
    integer,allocatable :: seed(:)

    call random_seed(size=n)
    allocate(seed(n))
    seed = 123456789    ! putting arbitrary seed to all elements
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine init_seed

  subroutine random_stduniform(u)
    implicit none
    real,intent(out) :: u
    real :: r
    call random_number(r)
    u = 1 - r
  end subroutine random_stduniform
  ! assuming a<b
  function random_uniform_int(a,b) result(x)
    implicit none
    integer,intent(in) :: a,b
    integer :: x
    real :: u
    call random_stduniform(u)
    x = nint(real(b-a)*u + real(a))
  end function random_uniform_int
  function random_uniform_real(a,b) result(x)
    implicit none
    real,intent(in) :: a,b
    real :: x
    real :: u
    call random_stduniform(u)
    x = (b-a)*u + a
  end function random_uniform_real
end module random_mod
