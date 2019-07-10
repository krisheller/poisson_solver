!
!  Unit:      linear_sys_solve_test_driver
!  Purpose:   Test the linear_sys_solve subroutine
!  Author:    F. D. Swesty
!  Date:      2/15/2018
!
program linear_sys_solve_test_driver
use linear_system_mod
implicit none
  
  integer, parameter :: n=3              ! Linear system size
  
  real(kind=dk) :: a(n,n)                ! Matrix

  real(kind=dk) :: b(n) = [2.0,3.0,4.0]  ! Right-hand-side of linear system

  real(kind=dk) :: x(n)                  ! Column vector of unknowns

                                         ! Rank 1 array used to init matrix
  real(kind=dk), parameter :: r1a(n*n) = [1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0]




  a = 2.0*reshape(r1a,[3,3]) ! Initialize A to 3x3 diagonal matrix
  b = [4.0,6.0,8.0]          ! Initialize b to a length three constant

                             ! Solve A*x = b for x
  call linear_sys_solve(a,b,x)

  write(*,*) ' x = ',x       ! Output x

  stop 0                     ! Stop with normal error code

end program linear_sys_solve_test_driver
