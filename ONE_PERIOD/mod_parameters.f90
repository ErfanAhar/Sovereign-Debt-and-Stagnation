module mod_parameters

implicit none

!----------------------------------------------------------------
! 1. Parameters
!----------------------------------------------------------------    

real(8) :: dlta,bta,ybarL,ybarH,kppa,sdz,gL,gH,pgL,pgH,psB,psG,b_lb,b_ub
real(8), dimension(2) :: ybarg		
real(8), parameter :: Rs    = 1.035d0                           ! risk-free rate (gross)
real(8), parameter :: rm    = Rs-1.0D0                          ! risk-free rate
real(8), parameter :: gma   = 3.0d0	                            ! coefficient of relative risk aversion (CRRA)
real(8), parameter :: cmin  = 1.0e-06	                        ! lower bound for consumption
real(8), parameter :: pthta = 0.90D0                            ! re-entry probability
real(8), parameter :: rhoub = 0.35D0                            ! for borrowing limit: implied rate must be lower than rhoub	
real(8)  		   :: wsunspot = 0                              ! wsunspot = 0 no sunspot, wsunspot = 1 for sunspot
real(8), parameter :: pdub  = 0.65D0							! max prob of default next period
integer, parameter :: ws = 2									! size of rolling window for interest rate schedule
integer, parameter :: load_guess = 0							! load_guess=0 if no guess available, load_guess=1 if guess available	


!----------------------------------------------------------------
! 2. Grids and function
!----------------------------------------------------------------    

! Parameters for the state space
integer, parameter :: Ng   = 2  								! number of growth states
integer, parameter :: Ns   = 2  								! number of sunspot states
integer, parameter :: Nz   = 17									! size of grid for transitory shock
integer, parameter :: Nb   = 3000                               ! size of the grid for debt service
real(8), parameter :: mz   = 3.5d0  							! Tauchen parameter

! Grids
real(8) :: zvec(Nz), Pz(Nz)						    ! grid of transitory shocks and respective probabilities
real(8) :: gvec(Ng)                                 ! grid of growth shocks
real(8) :: Y(Ng,Nz)									! grid of endowment realization
real(8) :: Pg(Ng,Ng), Ps(Ns,Ns)                     ! transition matrices: Pg(g'|g) and Ps(s'|s)
real(8) :: bvec(Nb)    								! grid of debt service

! value and policies
real(8), dimension(Nb,Ng,Nz) :: Vnd              ! value function no-default
real(8), dimension(Nb,Ng,Nz) :: Vd               ! value function default
real(8), dimension(Nb,Ng,Nz) :: Vndk             ! Vnd(kappa*b,g,s,z)
real(8), dimension(Nb,Ng,Nz) :: Wv               ! Wv = max{Vnd , Vd}
real(8), dimension(Nb,Ng,Nz) :: Wtv              ! Wv = max{Vndk, Vd}
real(8), dimension(Nb,Ng)    :: EWv              ! EWv(b,g,s)  = sum_{g',s',z'}{Wv(b,g',s',z') p(z')p(g'|g)p(s'|s)}
real(8), dimension(Nb,Ng)    :: EVd              ! EVd(b,g,s)  = sum_{g',s',z'}{Vd(b,g',s',z') p(z')p(g'|g)p(s'|s)}
real(8), dimension(Nb,Ng)    :: EWtv             ! EWtv(b,g,s) = sum_{g',s',z'}{Wtv(b,g',s',z')p(z')p(g'|g)p(s'|s)}    

integer, dimension(Nb,Ng,Nz) :: bopt_loc

real(8), dimension(Nb,Ng,Nz) :: bpol             ! debt policy: bpol = bvec(bopt_loc)	
real(8), dimension(Nb,Ng,Nz) :: ppol             ! probability of default policy	
real(8), dimension(Nb,Ng,Nz) :: rpol             ! interest rate policy		
real(8), dimension(Nb,Ng,Nz) :: dpol             ! default policy
real(8), dimension(Nb,Ng,Nz) :: epol             ! entry policy

real(8), dimension(Nb,Ng,Nz) :: Qprice           ! price of no defaulted bond
real(8), dimension(Nb,Ng,Nz) :: Xprice           ! price of defaulted bond
real(8), dimension(Nb,Ng)    :: EQ,EXX           ! EQ(b,g,s) = sum_{g',s',z'}{Q(b,g',s',z')p(z')p(g'|g)p(s'|s)}

real(8), dimension(Nb,Ng)    :: Rsched           ! R(b',g,s) -> the "full" schedule // note that b doesn't matter given b'
real(8), dimension(Nb,Ng)    :: rhosched         ! interest rate implied by R(b',g,s)	
real(8), dimension(Nb,Ng)    :: nsched           ! n(b',b,g,s) -> gb' = (1-dlta)*b + R(b',g,s)n(b',b,g,s)
real(8), dimension(Nb,Ng)    :: Ednp				! expected dfault next period		
integer, dimension(Nb,Nb,Ng) :: lsch,isch,lisch,msch

real(8) :: Qub, Qlb, Xub, Xlb

integer :: btog_vind(Nb,Ng), btogdm_vind(Nb,Ng), bkpa_ind(Nb)    

real(8), parameter :: pi = 3.14159265359D0			! pi

real(8) :: Vlow = -900000.0D0	                   

real(8) :: QQub, QQlb, XXub, XXlb

!----------------------------------------------------------------
! 3. Read input and save outputs
!----------------------------------------------------------------    
integer, parameter :: NCALIB = 11    ! number of calibrations to use
integer, parameter :: Nprm   = 13    ! number of parameters we may change, prMAT = [dlta,ybarL,ybarH,kppa,pgL,pgH,bta,gL,gH,sdz,psB,b_lb,b_ub]	
real(8)            :: prmMAT(NCALIB,Nprm)
!---to read inputs and save
character(len=5)     :: iWxs
integer              :: iWx 
character(len = 250) :: wvS
	



end module mod_parameters