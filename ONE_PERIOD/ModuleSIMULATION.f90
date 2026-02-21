module ModuleSIMULATION

    use mod_parameters
    use Toolbox
   
contains 

    !---------------------------------------------------
    !--- runs simulation
    !---------------------------------------------------
    subroutine SimulationRUN    
        
    implicit none   

    real(8) :: uu(2), ue(1), CDFg(Ng,Ng), CDFz(Nz), bbtp, Vbrx, Vbnx, nbrx
    integer :: tt, ig, iglp, igtp, iz, izlp, iztp, ib, ibtp, idef, idefnp, ibri, ibrx, ibdnv(1), ibdn
    integer :: count_re, count_nre, count_def

    integer, parameter :: Tsim = 2000000
    real(8), dimension(Tsim)   :: gsim, zsim, rhosim, nsim, ysim, esim, qsim, qbsim, psim
	real(8), dimension(Tsim+1) :: bsim, cnsim
	integer, dimension(Tsim)   :: igsim, izsim, iresim
	integer, dimension(Tsim+1) :: ibsim, idsim    

    !integer, parameter :: Tburn = 10000
    !real(8), dimension(Tsim-Tburn) :: aux_real_mom
    !integer, dimension(Tsim-Tburn) :: aux_int_mom

    print *, 'can you see me?'  

    print *, '0) random seed'
    call RANDOM_SEED()

    !-------------------------------------------
	!---CDFs for draws
	!-------------------------------------------
    print *, '1) computing CDFs'
	do iglp = 1, Ng
		CDFg(iglp,1) = Pg(iglp,1)		
		do igtp = 2, Ng
			CDFg(iglp,igtp) = Pg(iglp,igtp) + CDFg(iglp,igtp-1)	
		enddo
	enddo	
	
	CDFz(1) = Pz(1)
	do iz = 2, Nz
		CDFz(iz) = Pz(iz) + CDFz(iz-1) 		
	enddo

    !-------------------------------------------
	!---initial point
	!-------------------------------------------
    print *, '2) initial point'
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
	
	iglp = 1; izlp = 2;
	bsim(1) = bbtp; ibsim(1) = ibtp; idsim(1) = idef; 
	rhosim(1) = rhosched(ibtp,iglp)

    
	!-------------------------------------------
	!---simulation
	!-------------------------------------------	
    print *, '3) starting sim'
	count_re = 0
	count_nre = 0
	count_def = 0	
	
	do tt = 1, Tsim

        print *, 'sim tt : ', tt

        !---draw g and z: uu(1)->g, uu(2)->z
        call RANDOM_NUMBER(uu)

        do ig = 1, Ng
			if (uu(1)<=CDFg(iglp,ig)) then
				igtp = ig
				exit
			endif
		enddo
		gsim(tt) = gvec(igtp); igsim(tt) = igtp;		
		
		do iz = 1, Nz
			if (uu(2)<=CDFz(iz)) then
				iztp = iz
				exit
			endif
		enddo
		zsim(tt) = zvec(iztp); izsim(tt) = iztp;

        if (idef == 0) then ! have market access this period --> current debt is on grid
            ! price of debt and debt stock at market value (debt before default/reenter decision)
			qsim(tt)  = Qprice(ibtp,igtp,iztp);
			qbsim(tt) = qsim(tt)*bsim(tt);						

            if (dpol(ibtp,igtp,iztp) == 0) then  ! no default
                bsim(tt+1)  = bpol(ibtp,igtp,iztp)
				ibsim(tt+1) = bopt_loc(ibtp,igtp,iztp)					
				rhosim(tt)  = rhosched(ibsim(tt+1),igtp)
				psim(tt)    = ppol(ibtp,igtp,iztp)
				nsim(tt)    = nsched(ibsim(tt+1),igtp)				
				ysim(tt)    = Y(igtp,iztp)
				idefnp      = 0
		    elseif (dpol(ibtp,igtp,iztp) == 1) then ! default
				count_def = count_def+1
				! print *, 'idef = 0, def'
				bsim(tt+1)  = bsim(tt)/gsim(tt)				
				ibsim(tt+1) = -1
				rhosim(tt)  = 189.0D0 ! should use Xprice here
				psim(tt)    = 0.0D0
				nsim(tt)    = 0.0D0				
				ysim(tt)    = ybarg(igtp)*Y(igtp,iztp)
				idefnp      = 1
			else
				STOP
			endif  
            
	    elseif (idef == 1) then  ! don't have market access this period 
			call RANDOM_NUMBER(ue)  ! draw for re-enter
			! price of debt and debt stock at market value (debt before default/reenter decision)
			qsim(tt)  = Xprice(ibtp,igtp,iztp);
			qbsim(tt) = qsim(tt)*bsim(tt);									
			if (ue(1)<pthta) then  ! remain in exclusion
				! print *, 'idef = 1, no ent'
				bsim(tt+1)  = bsim(tt)/gsim(tt)
				ibsim(tt+1) = -1 
				rhosim(tt)  = 189.0D0 ! should use Xprice here
				psim(tt)    = 0.0D0
				nsim(tt)    = 0.0D0
				!cnsim(tt+1) = cnsim(tt)/gsim(tt)
				ysim(tt)    = ybarg(igtp)*Y(igtp,iztp)
				idefnp      = 1
			else   ! re-enter is possible
				!*value of re-enter*    function: Vnd_ent_grid(bb,ig,iz,ib_in,ib_out,Vb_out)
				call Vnd_ent_grid(kppa*bsim(tt),igtp,iztp,ibri,ibrx,Vbrx) 	
				
				!*value of no re-enter*				
				bbtp  = bsim(tt)
				ibdnv = vsearch_FET([bbtp],bvec); ibdn = ibdnv(1)
				Vbnx  = Vd(ibdn,igtp,iztp) + (bbtp - bvec(ibdn)) * ( Vd(ibdn+1,igtp,iztp) - Vd(ibdn,igtp,iztp)  )/(bvec(ibdn+1) - bvec(ibdn))
				
				if (Vbrx >= Vbnx) then ! re-enter is prefered
					count_re = count_re+1
					bsim(tt)    = kppa*bsim(tt)
					bsim(tt+1)  = bvec(ibrx)
					ibsim(tt+1) = ibrx
					rhosim(tt)  = rhosched(ibsim(tt+1),igtp)
					psim(tt)    = 0.0D0
					nsim(tt)    = nsched(ibsim(tt+1),igtp)					
					ysim(tt)    = Y(igtp,iztp)
					esim(tt)    = 1.0D0 
				    idefnp      = 0
				else	! no re-enter 		
					count_nre = count_nre+1
					bsim(tt+1)  = bsim(tt)/gsim(tt)
					ibsim(tt+1) = -1
					rhosim(tt)  = 189.0D0 ! should use Xprice here 
					psim(tt)    = 0.0D0
					nsim(tt)    = 0.0D0					
					ysim(tt)    = ybarg(igtp)*Y(igtp,iztp)
					esim(tt)    = 2.0D0 
				    idefnp      = 1
				endif		
			
			endif           
                
        endif

        !---next period states		
		iglp = igtp; 
		idef = idefnp
		ibtp = ibsim(tt+1)
		idsim(tt+1) = idefnp

    enddo !end loop tt


    ! !--- computing main moments (unconditional)
	! ! print *, '*************************************'
	! ! print *,'PARAMETERS  '
	! ! print *,'b_lb        ', b_lb
	! ! print *,'b_ub        ', b_ub
	! ! print *,'psB         ', psB
	! ! print *,'ybarL       ', ybarL
	! ! print *,'ybarH       ', ybarH
	! ! print *,'pgL         ', pgL
	! ! print *,'pgH         ', pgH
	! ! print *,'gL       	 ', gL
	! ! print *,'gH       	 ', gH
	! ! print *,'sdz       	 ', sdz
	! ! print *,'dlta        ', dlta	
	! ! print *,'bta         ', bta	
	! ! print *,'kppa        ', kppa	

    ! !---save simulation
	! open(1,  file = 'OUTPUT_'//trim(wvS)//'/bsim.txt'	  ,  status = 'unknown')	
	! open(2,  file = 'OUTPUT_'//trim(wvS)//'/cnsim.txt'	  ,  status = 'unknown')	
	! open(3,  file = 'OUTPUT_'//trim(wvS)//'/idsim.txt'	  ,  status = 'unknown')
	! open(4,  file = 'OUTPUT_'//trim(wvS)//'/ibsim.txt'	  ,  status = 'unknown')
	! do tt = 1, Tsim+1
	! 	write(1, '(*(f18.6))') ( bsim(tt)   )		
	! 	!write(2, '(*(f18.6))') ( cnsim(tt)   )	
	! 	write(3, '(I8.1)')     ( idsim(tt)  )
	! 	write(4, '(I8.1)')     ( ibsim(tt)  )
	! enddo
	! close(1); close(2);	close(3); close(4);
	
	! open(1,  file = 'OUTPUT_'//trim(wvS)//'/rhosim.txt'	  ,  status = 'unknown')
	! open(2,  file = 'OUTPUT_'//trim(wvS)//'/nsim.txt'	  ,  status = 'unknown')
	! open(3,  file = 'OUTPUT_'//trim(wvS)//'/gsim.txt'	  ,  status = 'unknown')
	! open(4,  file = 'OUTPUT_'//trim(wvS)//'/zsim.txt'	  ,  status = 'unknown')	
	! open(50,  file = 'OUTPUT_'//trim(wvS)//'/issim.txt'	  ,  status = 'unknown')	
	! open(60,  file = 'OUTPUT_'//trim(wvS)//'/ysim.txt'	  ,  status = 'unknown')
	! open(7,  file = 'OUTPUT_'//trim(wvS)//'/esim.txt'	  ,  status = 'unknown')
	! open(8,  file = 'OUTPUT_'//trim(wvS)//'/qsim.txt'	  ,  status = 'unknown')	
	! open(9,  file = 'OUTPUT_'//trim(wvS)//'/qbsim.txt'	  ,  status = 'unknown')		
	! open(10, file = 'OUTPUT_'//trim(wvS)//'/psim.txt'	  ,  status = 'unknown')		
	! open(11, file = 'OUTPUT_'//trim(wvS)//'/igsim.txt'	  ,  status = 'unknown')	
	! open(12, file = 'OUTPUT_'//trim(wvS)//'/izsim.txt'	  ,  status = 'unknown')	

	! do tt = 1, Tsim		
	! 	write(1, '(*(f18.6))') ( rhosim(tt) )
	! 	write(2, '(*(f18.6))') ( nsim(tt)   )
	! 	write(3, '(*(f18.6))') ( gsim(tt)   )
	! 	write(4, '(*(f18.6))') ( zsim(tt)   )	
	! 	!write(50, '(I8.1)')     ( issim(tt)  )		
	! 	write(60, '(*(f18.6))') ( ysim(tt)   )
	! 	write(7, '(*(f18.6))') ( esim(tt)   )
	! 	write(8, '(*(f18.6))') ( qsim(tt) )
	! 	write(9, '(*(f18.6))') ( qbsim(tt) )	
	! 	write(10,'(*(f18.6))') ( psim(tt) )							
	! 	write(11,'(I8.1)')     ( igsim(tt)  )
	! 	write(12,'(I8.1)')     ( izsim(tt)  )		
	! enddo

	! close(1); close(2); close(3); close(4); close(50); close(60); close(7); close(8); close(9); close(10); close(11); close(12);

    end subroutine SimulationRUN


    !-------------------------------------------------------
    !--- F1.grid. Vnd eval upon re-entry
    !-------------------------------------------------------
    subroutine Vnd_ent_grid(bb,ig,iz,ib_in,ib_out,Vb_out)            
	
	implicit none
	real(8) :: bb
	integer :: ig, iz
	
	integer :: ib_in, ib_out
	real(8) :: Vb_out
	real(8) :: gg, px(Nb)
	integer :: ib, ibx, ibv(1)
	
	gg   = gvec(ig); px = Ednp(:,ig);
	
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
	Vb_out = Vnd(ib,ig,iz)
	ib_out = bopt_loc(ib,ig,iz)
	
	
	end subroutine Vnd_ent_grid

end module ModuleSIMULATION