module ModuleSIMULATION

	use mod_parameters
	use Toolbox
	
contains

	subroutine SimulationRun
	
	real(8) :: uu(3), ue(1), CDFg(Ng,Ng), CDFs(Ns,Ns), CDFz(Nz), bbtp, Vbrx, Vbnx, nbrx
	integer :: tt, ig, iglp, igtp, is, islp, istp, iz, izlp, iztp, ib, ibtp, idef, idefnp, ibri, ibrx, ibdnv(1), ibdn
	integer :: count_re, count_nre, count_def
	integer, parameter :: Tsim = 2000000
	real(8), dimension(Tsim)   :: gsim, zsim, rhosim, nsim, ysim, esim, qsim, qbsim, psim
	real(8), dimension(Tsim+1) :: bsim, cnsim
	integer, dimension(Tsim)   :: igsim, issim, izsim, iresim
	integer, dimension(Tsim+1) :: ibsim, idsim
	integer :: aux_int
	real(8), dimension(Tsim-1000) :: aux_real_mom1, aux_real_mom2
	integer, dimension(Tsim-1000) :: aux_int_mom1, aux_int_mom2, aux_int_mom3	
	
	call RANDOM_SEED()
	
	!-------------------------------------------
	!---CDFs for draws
	!-------------------------------------------
	do iglp = 1, Ng
		CDFg(iglp,1) = Pg(iglp,1)		
		do igtp = 2, Ng
			CDFg(iglp,igtp) = Pg(iglp,igtp) + CDFg(iglp,igtp-1)	
		enddo
	enddo
	
	do islp = 1, Ns
		CDFs(islp,1) = Ps(islp,1)		
		do istp = 2, Ns
			CDFs(islp,istp) = Ps(islp,istp) + CDFs(islp,istp-1)	
		enddo	
	enddo
	
	CDFz(1) = Pz(1)
	do iz = 2, Nz
		CDFz(iz) = Pz(iz) + CDFz(iz-1) 		
	enddo
	
	!-------------------------------------------
	!---initial point
	!-------------------------------------------
	idef = 0   ! start w market access	
	
	!---initial b: lowest (positive) debt level on grid
	do ib = 1, Nb
		if (bvec(ib) > 0) then
			bbtp = bvec(ib)
			ibtp = ib
			exit		
		endif
	enddo
	
	nsim = 0.0D0; cnsim = 0.0D0; esim = 0.0D0
	
	iglp = 1; islp = 1; izlp = 2;
	bsim(1) = bbtp; ibsim(1) = ibtp; idsim(1) = idef; 
	rhosim(1) = rhosched(ibtp,iglp,islp)

	!-------------------------------------------
	!---simulation
	!-------------------------------------------	
	count_re = 0
	count_nre = 0
	count_def = 0	
	
	
	do tt = 1, Tsim
	
		!---draw g, s and z: uu(1)->g, uu(2)->s, uu(3)->z
		call RANDOM_NUMBER(uu)
		
		do ig = 1, Ng
			if (uu(1)<=CDFg(iglp,ig)) then
				igtp = ig
				exit
			endif
		enddo
		gsim(tt) = gvec(igtp); igsim(tt) = igtp;
		
		do is = 1, Ns
			if (uu(2)<=CDFs(islp,is)) then
				istp = is
				exit
			endif
		enddo 
		issim(tt) = istp
		
		do iz = 1, Nz
			if (uu(3)<=CDFz(iz)) then
				iztp = iz
				exit
			endif
		enddo
		zsim(tt) = zvec(iztp); izsim(tt) = iztp;
				
		!---compute policies
		if (idef == 0) then      ! have market access this period --> current debt is on grid
			! price of debt and debt stock at market value (debt before default/reenter decision)
			qsim(tt)  = Qprice(ibtp,igtp,istp,iztp);
			qbsim(tt) = qsim(tt)*bsim(tt);						
			if (dpol(ibtp,igtp,istp,iztp) == 0) then  ! no default
				! print *, 'idef = 0, no def'
				bsim(tt+1)  = bpol(ibtp,igtp,istp,iztp)
				ibsim(tt+1) = bopt_loc(ibtp,igtp,istp,iztp)					
				rhosim(tt)  = rhosched(ibsim(tt+1),igtp,istp)
				psim(tt)    = ppol(ibtp,igtp,istp,iztp)
				nsim(tt)    = nsched(ibsim(tt+1),ibtp,igtp,istp)
				!cnsim(tt+1) = (1.0D0-dlta)*cnsim(tt) + nsim(tt)
				cnsim(tt+1) = (1.0D0-dlta)*cnsim(tt)/gsim(tt) + nsim(tt)/gsim(tt)
				ysim(tt)    = Y(igtp,iztp)
				idefnp      = 0
			elseif (dpol(ibtp,igtp,istp,iztp) == 1) then ! default
				count_def = count_def+1
				! print *, 'idef = 0, def'
				bsim(tt+1)  = bsim(tt)/gsim(tt)				
				ibsim(tt+1) = -1
				rhosim(tt)  = (dlta*Rs/f_EX(bsim(tt)/gsim(tt),igsim(tt),issim(tt))) - dlta
				psim(tt)    = 0.0D0
				nsim(tt)    = 0.0D0
				cnsim(tt+1) = cnsim(tt)/gsim(tt)
				ysim(tt)    = ybarg(igtp)*Y(igtp,iztp)
				idefnp      = 1
			else
				STOP
			endif
			
		elseif (idef == 1) then  ! don't have market access this period 
			call RANDOM_NUMBER(ue)  ! draw for re-enter
			! price of debt and debt stock at market value (debt before default/reenter decision)
			qsim(tt)  = Xprice(ibtp,igtp,istp,iztp);
			qbsim(tt) = qsim(tt)*bsim(tt);									
			if (ue(1)<pthta) then  ! remain in exclusion
				! print *, 'idef = 1, no ent'
				bsim(tt+1)  = bsim(tt)/gsim(tt)
				ibsim(tt+1) = -1 
				rhosim(tt)  = (dlta*Rs/f_EX(bsim(tt)/gsim(tt),igsim(tt),issim(tt))) - dlta
				psim(tt)    = 0.0D0
				nsim(tt)    = 0.0D0
				cnsim(tt+1) = cnsim(tt)/gsim(tt)
				ysim(tt)    = ybarg(igtp)*Y(igtp,iztp)
				idefnp      = 1
			else   ! re-enter is possible
				!*value of re-enter*
				call Vnd_ent_grid(kppa*bsim(tt),igtp,istp,iztp,ibri,ibrx,Vbrx) 	
				
				!*value of no re-enter*				
				bbtp  = bsim(tt)
				ibdnv = vsearch_FET([bbtp],bvec); ibdn = ibdnv(1)
				Vbnx  = Vd(ibdn,igtp,istp,iztp) + (bbtp - bvec(ibdn)) * ( Vd(ibdn+1,igtp,istp,iztp) - Vd(ibdn,igtp,istp,iztp)  )/(bvec(ibdn+1) - bvec(ibdn))
				
				if (Vbrx >= Vbnx) then ! re-enter is prefered
					count_re = count_re+1
					bsim(tt)    = kppa*bsim(tt)
					bsim(tt+1)  = bvec(ibrx)
					ibsim(tt+1) = ibrx
					rhosim(tt)  = rhosched(ibsim(tt+1),igtp,istp)
					psim(tt)    = 0.0D0
					nsim(tt)    = nsched(ibsim(tt+1),ibri,igtp,istp)
					cnsim(tt)   = kppa*cnsim(tt)
					cnsim(tt+1) = (1.0D0-dlta)*cnsim(tt)/gsim(tt) + nsim(tt)/gsim(tt)
					ysim(tt)    = Y(igtp,iztp)
					esim(tt)    = 1.0D0 
				    idefnp      = 0
				else	! no re-enter 		
					count_nre = count_nre+1
					bsim(tt+1)  = bsim(tt)/gsim(tt)
					ibsim(tt+1) = -1
					rhosim(tt)  = (dlta*Rs/f_EX(bsim(tt)/gsim(tt),igsim(tt),issim(tt))) - dlta
					psim(tt)    = 0.0D0
					nsim(tt)    = 0.0D0
					cnsim(tt+1) = cnsim(tt)/gsim(tt)
					ysim(tt)    = ybarg(igtp)*Y(igtp,iztp)
					esim(tt)    = 2.0D0 
				    idefnp      = 1
				endif		
			
			endif
		endif	
	
		!---next period states		
		iglp = igtp; islp = istp;
		idef = idefnp
		ibtp = ibsim(tt+1)
		idsim(tt+1) = idefnp
	
	enddo  ! end loop tt
	
	!--- computing main moments (unconditional)
	print *, '*************************************'
	print *,'PARAMETERS  '
	print *,'b_lb        ', b_lb
	print *,'b_ub        ', b_ub
	print *,'psB         ', psB
	print *,'ybarL       ', ybarL
	print *,'ybarH       ', ybarH
	print *,'pgL         ', pgL
	print *,'pgH         ', pgH
	print *,'gL       	 ', gL
	print *,'gH       	 ', gH
	print *,'sdz       	 ', sdz
	print *,'dlta        ', dlta	
	print *,'bta         ', bta	
	print *,'kppa        ', kppa	

	aux_int_mom1 = idsim(1002:Tsim+1)
	aux_int_mom3 = igsim(1001:Tsim)

	print *, 'UNCONDITIONAL MOMENTS'
	aux_int_mom2 = 1
	where (aux_int_mom1>0.5) aux_int_mom2 = 0;
	aux_int = sum(aux_int_mom2)
 	aux_real_mom1 = rhosim(1001:Tsim);
 	where (aux_int_mom2<0.5) aux_real_mom1 = 0.0D0;
 	print*, 'avg(spread) ', 100*((sum(aux_real_mom1)/aux_int)-(Rs-1.0))

 	aux_real_mom1 = bsim(1001:Tsim)/ysim(1001:Tsim);
 	where (aux_int_mom2<0.5) aux_real_mom1 = 0.0D0;
 	print*, 'avg(b)      ', 100*sum(aux_real_mom1)/aux_int

 	aux_real_mom1 = qbsim(1001:Tsim)/dlta;	
 	where (aux_int_mom2<0.5) aux_real_mom1 = 0.0D0;
 	print*, 'avg(qb)     ', 100*sum(aux_real_mom1)/aux_int

 	aux_real_mom1 = cnsim(1001:Tsim)/ysim(1001:Tsim);	
 	where (aux_int_mom2<0.5) aux_real_mom1 = 0.0D0;
 	print*, 'avg(f)      ', 100*sum(aux_real_mom1)/aux_int

 	aux_real_mom1 = nsim(1001:Tsim)/ysim(1001:Tsim);	
 	where (aux_int_mom2<0.5) aux_real_mom1 = 0.0D0;
 	print*, 'avg(n)      ', 100*sum(aux_real_mom1)/aux_int

	!---save simulation
	open(1,  file = 'OUTPUT_'//trim(wvS)//'/bsim.txt'	  ,  status = 'unknown')	
	open(2,  file = 'OUTPUT_'//trim(wvS)//'/cnsim.txt'	  ,  status = 'unknown')	
	open(3,  file = 'OUTPUT_'//trim(wvS)//'/idsim.txt'	  ,  status = 'unknown')
	open(4,  file = 'OUTPUT_'//trim(wvS)//'/ibsim.txt'	  ,  status = 'unknown')
	do tt = 1, Tsim+1
		write(1, '(*(f18.6))') ( bsim(tt)   )		
		write(2, '(*(f18.6))') ( cnsim(tt)   )	
		write(3, '(I8.1)')     ( idsim(tt)  )
		write(4, '(I8.1)')     ( ibsim(tt)  )
	enddo
	close(1); close(2);	close(3); close(4);
	
	open(1,  file = 'OUTPUT_'//trim(wvS)//'/rhosim.txt'	  ,  status = 'unknown')
	open(2,  file = 'OUTPUT_'//trim(wvS)//'/nsim.txt'	  ,  status = 'unknown')
	open(3,  file = 'OUTPUT_'//trim(wvS)//'/gsim.txt'	  ,  status = 'unknown')
	open(4,  file = 'OUTPUT_'//trim(wvS)//'/zsim.txt'	  ,  status = 'unknown')	
	open(5,  file = 'OUTPUT_'//trim(wvS)//'/issim.txt'	  ,  status = 'unknown')	
	open(6,  file = 'OUTPUT_'//trim(wvS)//'/ysim.txt'	  ,  status = 'unknown')
	open(7,  file = 'OUTPUT_'//trim(wvS)//'/esim.txt'	  ,  status = 'unknown')
	open(8,  file = 'OUTPUT_'//trim(wvS)//'/qsim.txt'	  ,  status = 'unknown')	
	open(9,  file = 'OUTPUT_'//trim(wvS)//'/qbsim.txt'	  ,  status = 'unknown')		
	open(10, file = 'OUTPUT_'//trim(wvS)//'/psim.txt'	  ,  status = 'unknown')		
	open(11, file = 'OUTPUT_'//trim(wvS)//'/igsim.txt'	  ,  status = 'unknown')	
	open(12, file = 'OUTPUT_'//trim(wvS)//'/izsim.txt'	  ,  status = 'unknown')	

	do tt = 1, Tsim		
		write(1, '(*(f18.6))') ( rhosim(tt) )
		write(2, '(*(f18.6))') ( nsim(tt)   )
		write(3, '(*(f18.6))') ( gsim(tt)   )
		write(4, '(*(f18.6))') ( zsim(tt)   )	
		write(5, '(I8.1)')     ( issim(tt)  )		
		write(6, '(*(f18.6))') ( ysim(tt)   )
		write(7, '(*(f18.6))') ( esim(tt)   )
		write(8, '(*(f18.6))') ( qsim(tt) )
		write(9, '(*(f18.6))') ( qbsim(tt) )	
		write(10,'(*(f18.6))') ( psim(tt) )							
		write(11,'(I8.1)')     ( igsim(tt)  )
		write(12,'(I8.1)')     ( izsim(tt)  )		
	enddo

	close(1); close(2); close(3); close(4); close(5); close(6); close(7); close(8); close(9); close(10); close(11); close(12);
	
	end subroutine SimulationRun
	
	!-------------------------------------------------------
    !--- F1.grid. Vnd eval upon re-entry
    !-------------------------------------------------------
    subroutine Vnd_ent_grid(bb,ig,is,iz,ib_in,ib_out,Vb_out)            
	
	implicit none
	real(8) :: bb
	integer :: ig, is, iz
	
	integer :: ib_in, ib_out
	real(8) :: Vb_out
	real(8) :: gg, px(Nb)
	integer :: ib, ibx, ibv(1)
	
	gg   = gvec(ig); px = Ednp(:,ig,is);
	
	!---find closest b on grid
	if (bb<bvec(1)) then
		ib = 1
	elseif (bb>bvec(Nb)) then
		ib = Nb
	else
		ibv = vsearch_FET([bb],bvec);  ibx = ibv(1);
		if ( (bb-bvec(ibx)) .LE. (bvec(ibx+1)-bb) ) then
			ib = ibx
		else
			ib = ibx+1
		endif	
	
	endif
	
	ib_in  = ib
	Vb_out = Vnd(ib,ig,is,iz)
	ib_out = bopt_loc(ib,ig,is,iz)
	
	
	end subroutine Vnd_ent_grid
	
    !-------------------------------------------------------
    !--- 6. f_EX
    !-------------------------------------------------------
	! Interpolation/extrabolation of EX(b,ig,is), b = dbt service, ig = growth rate state, is = sunspot state
	function f_EX(x1f1,x2f1,x3f1)
	implicit none
	! 1. Declaring Variables
	! Inputs
	double precision :: x1f1
	integer :: x2f1, x3f1
	! New Variables
	integer :: if1, indf1
	! Output
	double precision :: f_EX
	! 2. Function
	if (x1f1 < bvec(1)) then
		f_EX = EXX(1,x2f1,x3f1) + ((EXX(2,x2f1,x3f1)-EXX(1,x2f1,x3f1))/(bvec(2)-bvec(1)))*(x1f1-bvec(1))
	else if (x1f1 > bvec(Nb)) then
		f_EX = EXX(Nb,x2f1,x3f1) + ((EXX(Nb,x2f1,x3f1)-EXX(Nb-1,x2f1,x3f1))/(bvec(Nb)-bvec(Nb-1)))*(x1f1-bvec(Nb))
	else
		indf1 = 0
		do if1 = 1,Nb
			if (bvec(if1)<=x1f1) indf1 = indf1+1
		end do
		f_EX = EXX(indf1,x2f1,x3f1) + ((EXX(indf1+1,x2f1,x3f1)-EXX(indf1,x2f1,x3f1))/(bvec(indf1+1)-bvec(indf1)))*(x1f1-bvec(indf1))
	end if

    f_EX = max(f_EX,XXlb)
    f_EX = min(f_EX,XXub)  

	end function		

end module ModuleSIMULATION