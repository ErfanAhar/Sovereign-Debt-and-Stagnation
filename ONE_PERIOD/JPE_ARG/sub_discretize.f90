subroutine sub_discretize

!subroutine sub_discretize(Pz,zvec)

! We discretize a normal distribution, z = sdz*e,  e~N(0,1)

use mod_parameters

implicit none

!----------------------------------------------------------------
! 1. Declaring Variables
!----------------------------------------------------------------

! ! Outputs
! double precision, dimension(Nz) :: zvec, Pz

! Other variables
double precision :: zvecmax, zvecmin
integer :: i5
double precision :: m1, m2

!----------------------------------------------------------------
! 2. Discretization
!----------------------------------------------------------------

zvecmax = mz*sdz
zvecmin = -mz*sdz

! Grids
do i5 = 1,Nz
	zvec(i5) = zvecmin + (zvecmax-zvecmin)*(i5-1)/(Nz-1)
end do

! Transition Matrix: Pz
Pz(1) = stdnormcdf((zvec(1)+zvec(2))/(2.0*sdz))
Pz(Nz) = 1 - stdnormcdf((zvec(Nz)+zvec(Nz-1))/(2.0*sdz))
do i5 = 2,(Nz-1)
	m1 = stdnormcdf((zvec(i5)+zvec(i5+1))/(2.0*sdz))
	m2 = stdnormcdf((zvec(i5)+zvec(i5-1))/(2.0*sdz))
	Pz(i5) = m1 - m2
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

function stdnormcdf(z)

implicit none

! 1. Declaring Variables

! Inputs
double precision :: z
! Output
double precision :: stdnormcdf

! 2. Function

stdnormcdf = 0.5*(1+erf(z/sqrt(2.0)))

if (z/sqrt(2.0) < -5.0*pi/2.0) stdnormcdf = 0.0
if (z/sqrt(2.0) > 5.0*pi/2.0) stdnormcdf = 1.0

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function erf(v)

implicit none

! 1. Declaring Variables

! Inputs
double precision :: v
! Output
double precision :: erf
! Other Variables
integer :: k5
double precision :: rk5

! 2. Function

erf = v/5.0
do k5 = 1,37
	rk5 = dfloat(k5)
	erf = erf + exp(-1.0*((rk5/5.0)**2.0))*sin(2.0*rk5*v/5.0)/rk5
end do
erf = (2.0/pi)*erf

end function


end subroutine sub_discretize
