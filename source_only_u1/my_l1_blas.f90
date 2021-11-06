module my_l1_blas
  implicit none
contains
  !
  subroutine my_zaxpy(n,a,x,y)
    integer, intent(in) :: n
    complex(8), intent(in) :: a, x(n)
    complex(8), intent(inout) :: y(n)
    integer :: i
    !$omp parallel do
    do i = 1, n
       y(i) = a*x(i) + y(i)
    end do
  end subroutine my_zaxpy
  !
  subroutine my_zdscal(n,a,x)
    integer, intent(in) :: n
    real(8), intent(in) :: a
    complex(8), intent(inout) :: x(n)
    integer :: i
    !$omp parallel do
    do i = 1, n
       x(i) = a*x(i)
    end do
  end subroutine my_zdscal
  !
end module my_l1_blas
