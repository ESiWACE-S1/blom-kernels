! https://masuday.github.io/fortran_tutorial/random.html
module random_mod
  implicit none
  private

  public :: random_uniform
contains
  subroutine random_stduniform(u)
    implicit none
    real,intent(out) :: u
    real :: r
    call random_number(r)
    u = 1 - r
  end subroutine random_stduniform
  ! assuming a<b
  function random_uniform(a,b) result(x)
    implicit none
    integer,intent(in) :: a,b
    integer :: x
    real :: u
    call random_stduniform(u)
    x = nint(real(b-a)*u + real(a))
  end function random_uniform
end module random_mod
