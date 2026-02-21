program MAIN
    
    use mod_parameters
    use FUNCTIONS
	use ModuleSIMULATION	
    use Toolbox
	use toolbox_CE
    use ModuleSAVE
    
    
    implicit none
    real(8) :: cxb(Nb), uxb(Nb), gg, btog(Nb), bx(Nb), bkpa(Nb), vz
	real(8) :: Vndnew(Nb,Ng,Ns,Nz), Wvold(Nb,Ng,Ns,Nz), Wtvold(Nb,Ng,Ns,Nz), epolold(Nb,Ng,Ns,Nz), dpolold(Nb,Ng,Ns,Nz)
	real(8) :: bpolold(Nb,Ng,Ns,Nz), aux_real(Nb,Ng,Ns,Nz)	
	integer :: bopt_old(Nb,Ng,Ns,Nz), aux_int(Nb,Ng,Ns,Nz)	
    integer :: ib, ig, is, iz
    real(8) :: vec_errW(5000)
		
    real(8) :: errW(1), errWt(1), wwVnd, tolVnd
    integer :: itVnd, maxitVnd, itsave
			
	!----------------------------------------------------------------
    ! 0. Parameters grid and allocate values
    !----------------------------------------------------------------
	! prMAT =      [dlta  , ybarL , ybarH , kppa, pgL   , pgH   , bta  , gL , gH   ,  sdz  , psB  , b_lb  , b_ub]	
	prmMAT(1,:)  = [0.40D0,0.971D0,0.90D0,0.75D0,0.60D0,0.75D0,0.75d0,0.96d0,1.04d0,0.023d0,0.25d0,-0.05d0,0.50d0] ! Benchmark
	prmMAT(2,:)  = [0.40D0,0.971D0,0.90D0,0.75D0,0.60D0,0.75D0,0.75d0,0.96d0,1.04d0,0.023d0,0.01d0,-0.05d0,0.50d0] ! low pB
	prmMAT(3,:)  = [0.15D0,0.971D0,0.90D0,0.75D0,0.60D0,0.75D0,0.75d0,0.96d0,1.04d0,0.023d0,0.25d0,-0.05d0,0.50d0] ! low delta
	prmMAT(4,:)  = [0.60D0,0.971D0,0.90D0,0.75D0,0.60D0,0.75D0,0.75d0,0.96d0,1.04d0,0.023d0,0.25d0,-0.05d0,0.60d0] ! high delta
	prmMAT(5,:)  = [0.40D0,0.971D0,0.90D0,0.75D0,0.65D0,0.75D0,0.75d0,0.96d0,1.04d0,0.023d0,0.25d0,-0.05d0,0.50d0] ! pL
	prmMAT(6,:)  = [0.40D0,0.971D0,0.90D0,0.70D0,0.60D0,0.75D0,0.75d0,0.96d0,1.04d0,0.023d0,0.25d0,-0.05d0,0.50d0] ! kappa
	prmMAT(7,:)  = [0.40D0,0.971D0,0.90D0,0.75D0,0.60D0,0.75D0,0.75d0,0.96d0,1.04d0,0.040d0,0.25d0,-0.05d0,0.50d0] ! sigma
	prmMAT(8,:)  = [0.40D0,0.971D0,0.92D0,0.75D0,0.60D0,0.75D0,0.75d0,0.96d0,1.04d0,0.023d0,0.25d0,-0.05d0,0.50d0] ! ybarH
	prmMAT(9,:)  = [0.40D0,0.971D0,0.90D0,0.75D0,0.60D0,0.70D0,0.75d0,0.96d0,1.04d0,0.023d0,0.25d0,-0.05d0,0.50d0] ! pH
	prmMAT(10,:) = [0.40D0,0.971D0,0.90D0,0.75D0,0.60D0,0.75D0,0.80d0,0.96d0,1.04d0,0.023d0,0.25d0,-0.05d0,0.50d0] ! beta
	prmMAT(11,:) = [0.40D0,0.971D0,0.90D0,0.75D0,0.60D0,0.75D0,0.75d0,0.9477d0,1.0477d0,0.023d0,0.25d0,-0.05d0,0.50d0] ! gH-gL

	!---read case
	if ( COMMAND_ARGUMENT_COUNT() /= 1) then		
		print *, 'ERROR: A command line argument specifying the parameter is required'
		STOP
	endif 
	call get_command_argument(1, iWxs)
	read(iWxs,*) iWx	
	
	wvS = iWxs

	! prMAT =      [dlta  , ybarL , ybarH , kppa, pgL   , pgH   , bta  , gL , gH   ,  sdz  , psB  ,psG   , b_lb  , b_ub]		
    dlta  = prmMAT(iWx,1);
	ybarL = prmMAT(iWx,2);
	ybarH = prmMAT(iWx,3);
	ybarg = [ybarL,ybarH]
	kppa  = prmMAT(iWx,4);
	pgL   = prmMAT(iWx,5);
	pgH   = prmMAT(iWx,6);
	bta   = prmMAT(iWx,7);
	gL    = prmMAT(iWx,8);
	gH    = prmMAT(iWx,9);		
	gvec = [gL, gH];
	sdz   = prmMAT(iWx,10);	
	psB   = prmMAT(iWx,11);
	psG   = 1.0d0-psB;
	b_lb  = prmMAT(iWx,12);	
	b_ub  = prmMAT(iWx,13);	
	vz    = sdz**2.0D0 
    
    !----------------------------------------------------------------
    ! 1. State Space and Transition Matrix
    !----------------------------------------------------------------
    
    !---discretizing the distribution of transitory shocks
	call sub_discretize(Pz,zvec)     
  
    !---bvec
    call linspace(bvec, b_lb, b_ub, Nb)	

    !---transition matrix
    Pg(1,1) = pgL; Pg(1,2) = 1.0-pgL; Pg(2,1) = 1.0-pgH; Pg(2,2) = pgH;	! growth
    Ps(1,1) = psB; Ps(1,2) = 1.0-psB; Ps(2,1) = 1.0-psG; Ps(2,2) = psG	! sunspot    
    
    !---output   
    ! gvec = [gL, gH]
    do iz=1,Nz
	do ig=1,Ng
		Y(ig,iz) = exp(zvec(iz))*gvec(ig)
	enddo
    enddo
    
    !---bounds on prices  (not imposing bounds right now)
    Qub = dlta*(Rs/(Rs-1.0D0+dlta)); Qlb = 0.0
    Xub = (kppa/Rs)*Qub;  Xlb = 0.0D0
    
    !---btog and bkpa indexes
    do ig = 1, Ng
         gg  = gvec(ig);  btog = bvec/gg
         btog_vind(:,ig) = vsearch_FET(btog,bvec); 
		 
		 bx = bvec/gg
		 btogdm_vind(:,ig) = vsearch_FET(bx,bvec); 
    enddo
    bkpa = kppa*bvec; bkpa_ind = vsearch_FET(bkpa,bvec);   	
	
	call SaveGrids	
    
    !----------------------------------------------------------------
    ! 2. Load initial guess
    !----------------------------------------------------------------
    
	if (load_guess==0) then
	
		!---Initial guess for Vnd
		do ig = 1, Ng
		do is = 1, Ns
			do iz = 1, Nz
				cxb = Y(ig,iz) - 0.10D0*bvec
				cxb = max(cxb, cmin)
				uxb = (cxb**(1.0D0-gma))/(1.0D0-gma)
				Vnd(:,ig,is,iz) = uxb/(1.0D0-bta)
			enddo
		enddo
		enddo
		!---Initial guess for bopt_loc
		do ig = 1, Ng
		do is = 1, Ns
		do iz = 1, Nz
			do ib = 1, Nb
				bopt_loc(ib,ig,is,iz) = ib
			enddo
		enddo
		enddo
		enddo

		!---Initial guess for Vd
		Vd = 0.0D0
		!---Initial guess for prices
		Qprice = Qub
		Xprice = Xub

	else	
	    
		!---load initial guess
		open(1,  file = 'OUTPUT_'//trim(wvS)//'/Vnd_in.txt',        status = 'unknown')
		open(2,  file = 'OUTPUT_'//trim(wvS)//'/bopt_loc_in.txt' ,  status = 'unknown')
		open(3,  file = 'OUTPUT_'//trim(wvS)//'/Vd_in.txt' ,  status = 'unknown')
		open(4,  file = 'OUTPUT_'//trim(wvS)//'/Qprice_in.txt' ,  status = 'unknown')
		open(5,  file = 'OUTPUT_'//trim(wvS)//'/Xprice_in.txt' ,  status = 'unknown')
		open(6,  file = 'OUTPUT_'//trim(wvS)//'/Wv_in.txt' ,  status = 'unknown')			
		read(1, *) Vnd
		read(2, *) bopt_loc
		read(3, *) Vd
		read(4, *) Qprice
		read(5, *) Xprice
		read(6, *) Wv	  
		close(1); close(2); close(3); close(4); close(5); close(6)
	
		wsunspot = 1
	
	endif
	
    !----------------------------------------------------------------
    ! 3. Iteration
    !----------------------------------------------------------------    

    maxitVnd = 4500; tolVnd = 1e-5; wwVnd = 0.45D0; itsave = 100 
	vec_errW = 1.0e3;
	
    do itVnd = 1, maxitVnd

		if (itVnd > 85) wsunspot = 1

		! Updating weights
		if (itVnd > 500) wwVnd = 0.850D0
		if (itVnd > 1000) wwVnd = 0.900D0
		if (itVnd > 1500) wwVnd = 0.950D0
		if (itVnd > 2000) wwVnd = 0.990D0
		if (itVnd > 2500) wwVnd = 0.995D0
		if (itVnd > 3500) wwVnd = 0.999D0
        
		Wvold    = Wv
		Wtvold   = Wtv
		dpolold  = dpol
		bpolold  = bpol
		epolold  = epol
		bopt_old = bopt_loc
		
        call COMPUTE_Vd
        call COMPUTE_QandX
        call COMPUTE_Rschedule
        
        do ig = 1, Ng
        do is = 1, Ns
        do iz = 1, Nz
            !$OMP PARALLEL DO
            do ib = 1, Nb
                Vndnew(ib,ig,is,iz) = Vnd_eval(ib,ig,is,iz)
            enddo
            !$OMP END PARALLEL DO 
        enddo    
        enddo    
        enddo     		
		
		errW  = maxval(abs(Wv-Wvold))
			
		vec_errW(itVnd) = errW(1);
		
        print *, '**********************************************'
        print *, '*** itVnd    = ', itVnd
		print *, '*** errW     =', errW
		print *, '*** min_errW =', minval(vec_errW(10:maxitVnd))
        print *, '**********************************************'
        
		Vnd = wwVnd*Vnd + (1.0D0-wwVnd)*Vndnew
		
		!---exit if Wv and Wtv converged
		if ( errW(1)<tolVnd .and. errWt(1)<tolVnd) then            
			print *, 'Wv and Wtv converged -- saving and exiting'	
			exit
		endif
        
        if (itsave <= 0 ) then
            print *, 'saving policies and schedule'
            call SaveVFandPolicies
            call SaveSCHEDULE
            itsave = 100
            print *, 'done saving'
        endif
        itsave = itsave -1
        
                
    enddo  ! end itVnd    
    
	call SaveVFandPolicies
    call SaveSCHEDULE	
	
	!----------------------------------------------------------------
    ! 4. Simulation
    !----------------------------------------------------------------    
	call SimulationRun	
	
end program MAIN

    
    