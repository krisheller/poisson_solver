  !
  !Program: problem2d
  !Author:  Kristofer Heller
  !Purpose: Achieve a second order accurate isoalted solution to the
  !         Poisson equation.
  !Created: 12/1/18
  !Edited:  4/18/19
  !

program problem2d
  use linear_system_mod
  implicit none

  !Parameters & Constants
  integer, parameter :: n = 8                            !The number of cells on each side of the grid
  integer, parameter :: dr=kind(1.0d0)                   !Precision type
  real(kind=dr), parameter :: PI=3.141592653589793d0     !Pi
  real(kind=dr), parameter :: sigma = 1.0d0

  !Grid elements
  real(kind=dr) :: rhoH(0:n+1,0:n+1), rhoC(0:n+1, 0:n+1)                      !The charge distribution
  real(kind=dr) :: phiH(0:n+1,0:n+1), phiC(0:n+1, 0:n+1)  !The output potential of first step
  real(kind=dr) :: phi_iso(0:n+1, 0:n+1)
  real(kind=dr) :: phi_analytic(0:n+1,0:n+1)             !The discretized analytic solution
  real(kind=dr) :: phi_boundaryZero(n)
  real(kind=dr) :: phi_n(n), phi_s(n), phi_e(n), phi_w(n)
  real(kind=dr) :: phi_an(n), phi_as(n), phi_ae(n), phi_aw(n)
  real(kind=dr) :: xMax = 1.0d0, xMin = -1.0d0
  real(kind=dr) :: yMax = 1.0d0, yMIn = -1.0d0
  real(kind=dr) :: dx, dy
  real(kind=dr) :: iC, jC

  !Helpers and Misc. Variables
  integer :: i, j, k                                     !Looping indices
  integer :: LUN1, LUN2, LUN3, LUN4, LUN5, LUN6, LUN7, LUN8                !I/O LUNs
  real(kind=dr) :: sum
  real(kind=dr) :: x2,y2
  real(kind=dr) :: sC, squareCharge                                    !Surface Charge sum
  real(kind=dr) :: leftCharge(n),rightCharge(n),topCharge(n),bottomCharge(n)
  real(kind=dr) :: leftSum, rightSum, bottomSum, topSum

  !Functions
  real(kind=dr) :: analytic                              !Analytic Function
  integer :: kDel                                        !Kronecker Delta Function
  real :: analyticZero                                   !Zero Potential Function


  !Set our elements to zero.

  rhoH = 0.0d0
  phiH = 0.0d0
  phiC = 0.0d0
  rhoC = 0.0d0
  phi_analytic = 0.0d0
  leftCharge = 0.0d0
  rightCharge = 0.0d0
  bottomCharge = 0.0d0
  topCharge = 0.0d0
  phi_boundaryZero = 0.0d0
  phi_w = 0.0d0
  phi_e = 0.0d0
  phi_s = 0.0d0
  phi_n = 0.0d0
  phi_ae = 0.0d0
  phi_aw = 0.0d0
  phi_as = 0.0d0
  phi_an = 0.0d0
  squareCharge = 0.0d0

  !Define some constants

  dx = (xMax - xMin)/real(n)
  dy = (yMax - yMin)/real(n)

  iC = real(n)/2.0d0 + 0.5d0
  jC = real(n)/2.0d0 + 0.5d0

  !Place the charge inside the rho matrix.

  rhoH(1+n/4:n-n/4, 1+n/4:n-n/4) = sigma

  do i = 1,n
     do j =1,n
        squareCharge = squareCharge + rhoH(i,j)*dx*dy
     enddo
  enddo

  !Write the analytic solution to the phi_analytic matrix.

  do i = 1, n
     do j = 1, n
        x2 = (i - iC)*dx
        y2 = (j - jC)*dy
        phi_analytic(i,j) = analytic(x2,y2,sigma)
     enddo
  enddo


  !Write the analytic solution to their respective boundary vectors.
  do i = 1,n
     x2 = (i-iC)*dx
     phi_an(i) = analytic(x2,yMax,sigma)
     phi_as(i) = analytic(x2,yMin,sigma)
     phi_ae(i) = analytic(xMax,y2,sigma)
     phi_aw(i) = analytic(xMin,y2,sigma)
  enddo

  !Call the poisson solver

  call poisson(n, rhoH, phi_boundaryZero, phi_boundaryZero, &
       &       phi_boundaryZero, phi_boundaryZero, phiH)

  !The first call works perfectly, homogeneous solution is correct.


  !Work with the surface charges caused by the zero potential.

  sC = 0.0d0

  !We iterate over the boundaries and add the surface charges up.

  do i = 1, n
     !Left Boundary
     leftCharge(i) = phiH(i,1)/(0.5d0*dx)
     !Right Boundary
     rightCharge(i) = phiH(i,n)/(0.5d0*dx)
     !Top Boundary
     topCharge(i) = phiH(n,i)/(0.5d0*dy)
     !Bottom Boundary
     bottomCharge(i) = phiH(1,i)/(0.5d0*dy)
  enddo


  !Get total charge with midpoint rule
  leftSum = 0.0d0
  rightSum = 0.0d0
  topSum = 0.0d0
  bottomSum = 0.0d0

  do i = 1,n
     leftSum = leftSum + leftCharge(i)*dy
     rightSum = rightSum + rightCharge(i)*dy
     topSum = topSum + topCharge(i)*dx
     bottomSum = bottomSum + bottomCharge(i)*dx
  enddo

  sC = topSum + bottomSum + rightSum + leftSum

  !Total surface charge is correct.

  !Now that we have the surface charges, we can calculate the potentials on the boundaries
  !due to these surface charges.  Iterate over every other charge, using midpoint
  !rule to calculate.

  do i = 1, n
     do j = 1, n
        !Contributions from same side
        if(i.ne.j) then
           phi_n(i) = phi_n(i) + topCharge(j)*dx*LOG((real((i-j))*dy)**2)
           phi_s(i) = phi_s(i) + bottomCharge(j)*dx*LOG((real((i-j))*dy)**2)
           phi_w(i) = phi_w(i) + leftCharge(j)*dy*LOG((real((i-j))*dy)**2)
           phi_e(i) = phi_e(i) + rightCharge(j)*dy*LOG((real((i-j))*dy)**2)
        endif

        !Contributions from other sides
        phi_n(i) = phi_n(i) + leftCharge(j)*dy*LOG((real(0.5d0-i)*dx)**2 + (real(n+0.5d0-j)*dy)**2)
        phi_n(i) = phi_n(i) + rightCharge(j)*dy*LOG((real(n+0.5d0-i)*dx)**2 + (real(n+0.5d0-j)*dy)**2)
        phi_n(i) = phi_n(i) + bottomCharge(j)*dx*LOG((real(i-j)*dx)**2 + (real(n)*dy)**2)
        phi_s(i) = phi_s(i) + leftCharge(j)*dy*LOG((real(0.5d0-i)*dx)**2 + (real(j-0.5d0)*dy)**2)
        phi_s(i) = phi_s(i) + rightCharge(j)*dy*LOG((real(n+0.5d0-i)*dx)**2 + (real(j-0.5d0)*dy)**2)
        phi_s(i) = phi_s(i) + topCharge(j)*dx*LOG((real(i-j)*dx)**2 + (real(n)*dy)**2)
        phi_w(i) = phi_w(i) + topCharge(j)*dx*LOG((real(j-0.5d0)*dx)**2 + (real(n+0.5d0-i)*dy)**2)
        phi_w(i) = phi_w(i) + bottomCharge(j)*dx*LOG((real(j-0.5d0)*dx)**2 + (real(i-0.5d0)*dy)**2)
        phi_w(i) = phi_w(i) + rightCharge(j)*dy*LOG((real(n)*dx)**2 + (real(i-j)*dy)**2)
        phi_e(i) = phi_e(i) + topCharge(j)*dx*LOG((real(n+0.5d0-j)*dx)**2 + (real(n+0.5d0-i)*dy)**2)
        phi_e(i) = phi_e(i) + bottomCharge(j)*dx*LOG((real(n+0.5d0-j)*dx)**2 + (real(i-0.5d0)*dy)**2)
        phi_e(i) = phi_e(i) + leftCharge(j)*dy*LOG((real(n)*dx)**2 + (real(i-j)*dy)**2)
     enddo
  enddo

  phi_n = phi_n/(4.0d0*PI)
  phi_s = phi_s/(4.0d0*PI)
  phi_e = phi_e/(4.0d0*PI)
  phi_w = phi_w/(4.0d0*PI)

  !Boundary potentials are calculated correctly.


  !These potentials are then passed into the poisson solver again with rho=0 everywhere inside.

  call poisson(n, rhoC, phi_an, phi_as, phi_ae, phi_aw, phiC)

  do i = 1,n
     !     write(*,*)phi_an(i), phi_n(i)
  enddo


  !The isolated solution is the first phi minus the second.

  phi_iso = phiH - phiC

  !END MAIN PROGRAM
  !ERROR ANALYSIS / OUTPUT OF DATA TO FILES

  !Error analysis of first potential

  sum = 0.0d0

  do i = 1, n
     do j = 1, n
        sum = sum + (phi_iso(i,j) - phi_analytic(i,j))**2
     enddo
  enddo

  sum = sqrt(sum/dble(n**2))

  write(*,*)"Log of N", log(dble(n))
  write(*,*)"Log of L2 Norm of System:", log(sum)

  !Write rho to a file

  open(newunit=LUN1, file='rho.dat', status='replace')
  do i = 1,n
     do j = 1,n
        write(LUN1,*)i,j,rhoH(i,j)
     enddo
  enddo
  close(unit=LUN1)

  !Write the potential of the homogeneous solution to a file

  open(newunit=LUN2, file='phiH.dat', status='replace')
  do i =1,n
     do j=1,n
        write(LUN2,*)i,j,phiH(i,j)
     enddo
  enddo
  close(unit=LUN2)

  !Write the potential of the charge potential to a file

  open(newunit=LUN3, file='phiC.dat', status='replace')
  do i=1,n
     do j=1,n
        write(LUN3,*)i,j,phiC(i,j)
     enddo
  enddo
  close(unit=LUN3)

  !Write the potential of the analytic potential to a file

  open(newunit=LUN4, file='pot_analytic.dat', status='replace')
  do i=1,n
     do j=1,n
        write(LUN4,*)i,j,phi_analytic(i,j)
     enddo
  enddo
  close(unit=LUN4)

  !Write the isolated potential to a file

  open(newunit=LUN5, file='phi_iso.dat', status='replace')
  do i=1,n
     do j=1,n
        write(LUN5,*)i,j,phi_iso(i,j)
     enddo
  enddo
  close(unit=LUN5)


  !More intermediate testing
  !Plot potential along the center in the x-direction vs. analytic
  open(newunit=LUN6, file='test1.dat', status='replace')
  do i =1,n
     write(LUN6,*)i, phi_iso(i, int(n/2 +1)),phi_analytic(i, int(n/2 +1))
  enddo
  close(unit=LUN6)

  open(newunit=LUN7, file='test2.dat', status='replace')
  do i =1,n
     write(LUN7,*)i, phi_iso(i, i), phi_analytic(i, i)
  enddo
  close(unit=LUN7)

  open(newunit=LUN8, file='test3.dat', status='replace')
  do i = 1,n
     write(LUN8,*) i, phi_n(i), phi_an(i)
  enddo
  close(unit=LUN8)


  stop 0
end program problem2d

!Begin Poisson Solver Subroutine
subroutine poisson(n, rho, phi_n, phi_s, phi_e, phi_w, phi)
  use linear_system_mod
  implicit none

  !Take in values
  integer, intent(in):: n
  real(kind(1.0d0)), intent(in) :: phi_n(n), phi_s(n), phi_e(n), phi_w(n)
  real(kind(1.0d0)), intent(in) :: rho(0:n+1, 0:n+1)

  !Declare the output
  real(kind(1.0d0)), intent(out) :: phi(0:n+1,0:n+1)

  !Create some necessary variables
  real(kind(1.0d0)) :: a(n**2, n**2)
  real(kind(1.0d0)) :: b(n**2)
  real(kind(1.0d0)) :: x(n**2)

  !Bring in the functions necessary
  real(kind(1.0d0)) :: analytic
  integer :: kDel
  real :: analyticZero

  !Make some constants
  real(kind(1.0d0)), parameter :: PI = 3.141592653589793d0
  real(kind(1.0d0)), parameter :: sigma = 1.0d0
  real(kind(1.0d0)), parameter :: xMax = 1.0d0, xMin = -1.0d0
  real(kind(1.0d0)), parameter :: yMax = 1.0d0, yMin = -1.0d0
  real(kind(1.0d0)) :: dx, dy
  real(kind(1.0d0)) :: iC, jC

  !Bring in some helpers
  real(kind(1.0d0)) :: x2,y2
  integer :: i,j,k

  !Clear necessary variables
  a = 0.0d0
  b = 0.0d0
  x = 0.0d0

  !Define some others
  dx = (xMax - xMin)/real(n)
  dy = (yMax - yMin)/real(n)
  iC = real(n)/2.0d0 + 0.5d0
  jC = real(n)/2.0d0 + 0.5d0


  !Define the coefficient matrix
  do i = 1, n
     do j = 1, n
        k = i + n*(j-1)
        !Define the major diagonal
        a(k,k) = -4.0d0 - dble(kDel(i,n) + kDel(i,1) + kDel(j,1) + kDel(j,n))

        !Define the inner minor diagonal
        if (k < n**2) then
           a(k,k+1) = 1.0d0*(1-kDel(i,n))
        endif
        if (k > 1) then
           a(k,k-1) = 1.0d0*(1-kDel(i,1))
        endif

        !Define the outer diagonal
        if (k <= n**2 - n) then
           a(k,k+n) = 1.0d0*(1-kDel(j,n))
        endif
        if (k > n) then
           a(k,k-n) = 1.0d0*(1-kDel(j,1))
        endif
     enddo
  enddo

  !Define the b vector
  do i = 1, n
     do j = 1, n
        k = j + n*(i-1)
        b(k) = rho(i,j)*(dx*dy)
        !Add in the boundary conditions
        if (kDel(i,1) == 1) then
           b(k) = b(k) + 2.0d0*phi_w(j)
        endif
        if (kDel(i,n) == 1) then
           b(k) = b(k) + 2.0d0*phi_e(j)
        endif
        if (kDel(j,n) == 1) then
           b(k) = b(k) + 2.0d0*phi_n(i)
        endif
        if (kDel(j,1) == 1) then
           b(k) = b(k) + 2.0d0*phi_s(i)
        endif
     enddo
  enddo

  !Solve the system
  call linear_sys_solve(a,b,x)

  !Convert x into phi
  phi = 0.0d0
  do i = 1, n
     do j = 1, n
        k = j + n*(i-1)
        phi(i,j) = x(k)
     enddo
  enddo


  return
end subroutine poisson


!Begin Functions

!Analytic Solution
function analytic(u,v,sigma) result(a)
  implicit none
  real(kind(1.0d0)),intent(in) :: u,v,sigma
  real(kind(1.0d0)) :: x0,y0,x1,y1,a
  real(kind(1.0d0)), parameter :: PI = 3.141592653589793d0
  x0 = -.5d0 - u
  x1 = .5d0 - u
  y0 = -.5d0 - v
  y1 = .5d0 - v

  a = -(sigma/(4.0*PI))*(-3.0d0*x0*y0 - 3.0d0*x1*y1 + 3.0d0*x1*y0 + 3.0d0*x0*y1 &
       &       +x0*y0*LOG(x0**2 + y0**2) + x1*y1*LOG(x1**2 + y1**2) - &
       &       x0*y1*LOG(x0**2 + y1**2) - x1*y0*LOG(x1**2 + y0**2) + &
       &       (y0**2)*ATAN(x0/y0) + (y1**2)*ATAN(x1/y1) - (y0**2)*ATAN(x1/y0) - &
       &       (y1**2)*ATAN(x0/y1) + (x0**2)*ATAN(y0/x0) + (x1**2)*ATAN(y1/x1) - &
       &       (x0**2)*ATAN(y1/x0) - (x1**2)*ATAN(y0/x1))

  return
end function analytic

!Kronecker Delta Function
integer function kDel(i,j) result(a)
  implicit none
  integer, intent(in) :: i,j
  a = 0
  if (i == j) then
     a = 1
  endif
  return
end function kDel


!Zero Potential Boundary Function
real function analyticZero(u,v,sigma) result(a)
  implicit none
  real(kind(1.0d0)), intent (in) :: u,v,sigma
  a = 0.0
  return
end function analyticZero
