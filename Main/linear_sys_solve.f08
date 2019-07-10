!
!  Unit:      linear_system_mod
!  Purpose:   Module containing linear system solver
!  Author:    F. D. Swesty
!  Date:      2/15/2018
!
!  Usage:     use linear_system_mod
!
!  Contains:  linear_sys_solve   (subroutine)
!
module linear_system_mod
  implicit none

  integer, parameter :: dk=kind(1.0d0)    ! Real "kind" for double precision

  integer, allocatable :: ipivot(:)       ! Pivot vector needed by LAPACK

  real(kind=dk), allocatable :: lu(:,:)   ! matrix needed to hold A and 
                                          ! the LU decomposition of A

  contains

!   Subroutine:   linear_sys_solve
!   Purpose:      Solve linear system Ax=b for unknown column vector x
!   Call line:    linear_sys_solve(a,b,x)
!   Arguments: 
!                 a = input matrix (rank 2, double precision)
!                 b = input right-hand-side (rank 1, double precision)
!                 x = output solution column vector (rank 1, double precision)
!
!   Note:         arguments a, b, & x must be of same extent in all dimensions
    subroutine linear_sys_solve(a,b,x)
    implicit none
    real(kind=dk), intent(in) :: a(:,:)  ! Input matrix
    real(kind=dk), intent(in) :: b(:)    ! Input right-hand-side
    real(kind=dk), intent(out) :: x(:)   ! Output soution vector

    integer :: nadim1, nadim2, nb, nx    ! Extents of array arguments
    integer :: info                      ! Error code variable for LAPACK

    nadim1 = size(a,dim=1)               ! Extent of A in dimension 1
    nadim2 = size(a,dim=2)               ! Extent of A in dimension 2
    nb = size(b,dim=1)                   ! Extent of b
    nx = size(x,dim=1)                   ! Extent of x

                                         ! If extents are all OK then...
    if(nadim1 == nadim2 .and. nb == nadim1 .and. nx == nadim1 ) then

                                         ! Solve linear system
      write(*,*) ' solving ',nadim1,'x',nadim2,' linear system'

                                         ! Allocate pivot array for LAPACK if 
                                         ! it has not yet been done
      if(.not.allocated(ipivot)) allocate(ipivot(nadim1))

                                         ! Allocate LU array for LAPACK if 
                                         ! it has not yet been done
      if(.not.allocated(lu)) allocate(lu(nadim1,nadim2))

 
      x = b                              ! Copy b into x vector
      lu = a                             ! Copy a into LU matrix

                                         ! Solve the system
      call dgesv(nadim1,1,lu,nadim1,ipivot,x,nx,info)

      if(info /= 0) then                 ! If icode is not zero...
        write(*,*) ' LAPACK routine dgetrs unable to solve linear system'
        write(*,*) ' LAPACK routine dgetrs info code ',info
        stop 4                             ! Stop with abnormal condition
      endif

      return                             ! All successful so let's return!

    endif
 
    if(nadim1 /= nadim2) then            ! If matrix is not square...
      write(*,*) ' matrix must be square: dim-1 is ',nadim1,' dim=2 is ',nadim2
      stop 1
    endif

    if(nx /= nadim1) then                ! If A & x are not same size...
      write(*,*) ' x vector is length ',nx, ' but must be ',nadim1
      stop 2
    endif

    if(nb /= nadim1) then                ! If A & b are not same size...
      write(*,*) ' b vector is length ',nb, ' but must be ',nadim1
      stop 3
    endif

    return                               ! All done so let's return

  end subroutine linear_sys_solve

end module linear_system_mod
