module Toolbox
    
    use mod_parameters
    
contains

    
    subroutine sub_discretize(Pz,zvec)
    ! We discretize a normal distribution, z = sdz*e,  e~N(0,1)
    implicit none
    !----------------------------------------------------------------
    ! 1. Declaring Variables
    !----------------------------------------------------------------
    ! Outputs
    double precision, dimension(Nz) :: zvec, Pz
    
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
    Pz(1)  = stdnormcdf((zvec(1)+zvec(2))/(2.0*sdz))
    Pz(Nz) = 1 - stdnormcdf((zvec(Nz)+zvec(Nz-1))/(2.0*sdz))
    do i5 = 2,(Nz-1)
        m1 = stdnormcdf((zvec(i5)+zvec(i5+1))/(2.0*sdz))
        m2 = stdnormcdf((zvec(i5)+zvec(i5-1))/(2.0*sdz))
        Pz(i5) = m1 - m2
    end do
    
    end subroutine sub_discretize
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                   !!!!!    
!!!!!     Search points in a vector     !!!!!   
!!!!!                                   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! vsearch(X,v) returns an integer with the two closest points in v to X
!
! Input: 
!   X = Nx x 1 vector
!   v = p  x 1 vector
!
! Output:
!   ind = Nx x 1 vector
!
! Note:
! if X(i) <= v(1) : ind(i) = 1
! if X(i) >= v(p) : ind(i) = p-1
! else            : v(ind(i)) <= X(i) < v(ind(i)+1)
!
! Please contact us if you find any errors:
! Axelle Ferriere, European University Institute: axelle.ferriere@eui.eu
! Gaston Navarro, Federal Reserve Board: gaston.m.navarro@frb.gov
! 
! First version: July 2016
! NEEDS UPDATE

function vsearch_FET(X,v) result(ind)

    real(8), intent(in)           :: X(:), v(:)
    integer, dimension(size(X,1)) :: ind
    
    integer, dimension(size(X,1)) :: i_lb, i_ub, i_mid, flag
    real(8), dimension(size(X,1)) :: v_mid 
    integer :: Nx, p, ssss
    
    p = size(v,1)
    
    flag = 1
    ssss=0d0
    where( X >= v(p) )
        
   
        ind  = p  -1
        flag = 0 
    end where
    where(X <= v(1))
        ind  = 1
        flag = 0
    end where    
    
    i_lb = 1
    i_ub = p
    i_mid = floor( dble(i_lb + i_ub)/2.0d0 )    
    
    v_mid = v(i_mid)    
   ! ssss=0
    do while(maxval(flag)>0) 
    ! ssss=   ssss+1
    ssss=ssss+1
    Sdrama=0d0
    if (ssss >100) stop 'Something is wrong in vsearch, maybe NAN'   
    ! print*, ssss
        where(flag == 1 .and. X < v_mid) 
            i_ub  = i_mid
            i_mid = ceiling( dble(i_lb+i_ub)/2.0d0 ) 
            v_mid = v(i_mid)
        end where
        where(flag == 1 .and. X >= v_mid) 
            i_lb  = i_mid
            i_mid = floor( dble(i_lb+i_ub)/2.0d0)
            v_mid = v(i_mid)
        end where
        where(i_ub == i_lb + 1) 
            flag = 0
            ind = i_lb
        end where

    end do
    
end function vsearch_FET   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                   !!!!!    
!!!!!     		  LINSPACE     		    !!!!!   
!!!!!                                   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Oliko is to be balmed ...

subroutine linspace(x, x_start, x_end, x_len)
	implicit none
	real(8), dimension(:), intent(out) :: x
	real(8) :: x_start, x_end, xd
	integer :: x_len, i

	if (x_len<1) then
		stop 'Error linspace: cannot build the vector'
	elseif (x_len == 1) then
		if (x_start == x_end) then
			x = x_start
		else
			stop 'Error linspace: cannot build the vector'
		endif
	else
		xd=(x_end-x_start)/dble(x_len-1)
		x(1:x_len)=[(x_start+xd*dble(i-1), i=1 ,x_len)]
	endif
end subroutine


    
end module Toolbox
