module FUNCTIONS
    
    use mod_parameters
    use Toolbox
    
contains


    !-------------------------------------------------------
    !--- 1. compute Vd (given Vnd)
    !-------------------------------------------------------
    subroutine COMPUTE_Vd
    
    implicit none
    real(8) :: gg, bnp(Nb), Vd_np(Nb), Wtv_np(Nb), Vx_np(Nb), bkpa(Nb), cxb(Nb), uxb(Nb)
    integer :: ib, ig, is, iz, ig_np, is_np, iz_np, bind_np(Nb) 
    real(8) :: Vdnew(Nb,Ng,Ns,Nz), errVd(1), tolVd, wwVd
    integer :: itVd, maxitVd
    
	maxitVd = 500; tolVd = 1e-8; wwVd = 0.25D0
	
	!---compute Vnd(kppa(b,g,s,z)
	bkpa = kppa*bvec   ! bkpa_ind = vsearch_FET(bkpa,bvec);   
	do ig = 1, Ng
	do is = 1, Ns
	do iz = 1, Nz
		Vndk(:,ig,is,iz) = Vnd(bkpa_ind,ig,is,iz) + (bkpa - bvec(bkpa_ind)) * ( Vnd(bkpa_ind+1,ig,is,iz) - Vnd(bkpa_ind,ig,is,iz) )/(bvec(bkpa_ind+1) - bvec(bkpa_ind))
	enddo
	enddo
	enddo
    
    
    do itVd = 1, maxitVd
        
        Wv  = max(Vnd , Vd)	
		Wtv = max(Vndk, Vd)    
		
        call COMPUTE_EWv
        
        do ig = 1, Ng
        do is = 1, Ns            
            
			gg  = gvec(ig)
			bnp = bvec/gg
			bind_np = btogdm_vind(:,ig)   ! = vsearch_FET(bx,bvec); with bx = bvec/gg
			
            Vd_np  = EVd(bind_np,ig,is)  + (bnp - bvec(bind_np)) * ( EVd(bind_np+1,ig,is)  - EVd(bind_np,ig,is)  )/(bvec(bind_np+1) - bvec(bind_np))
            Wtv_np = EWtv(bind_np,ig,is) + (bnp - bvec(bind_np)) * ( EWtv(bind_np+1,ig,is) - EWtv(bind_np,ig,is) )/(bvec(bind_np+1) - bvec(bind_np))			
			
            Vx_np = bta * (gg**(1.0D0-gma)) * (pthta*Vd_np + (1.0D0-pthta)*Wtv_np)                        
            
            do iz = 1, Nz
                cxb = ybarg(ig)*Y(ig,iz)
                cxb = max(cxb, cmin)
                uxb = (cxb**(1.0D0-gma))/(1.0D0-gma)
                
                Vdnew(:,ig,is,iz) = uxb + Vx_np
            enddo
            
        enddo        
        enddo
        
        errVd = maxval(abs(Vdnew - Vd))        

        if (errVd(1) < tolVd) then
            exit
        else
            Vd = wwVd*Vd + (1.0D0-wwVd)*Vdnew
        endif
        
    enddo    
    
    !---compue dpol and epol
    dpol = 0.0D0; epol = 0.0D0; 
    do ib = 1, Nb
    do ig = 1, Ng
    do is = 1, Ns
    do iz = 1, Nz
        if (Vd(ib,ig,is,iz) > Vnd(ib,ig,is,iz) ) then
            dpol(ib,ig,is,iz) = 1.0D0
        endif
        
        if (Vndk(ib,ig,is,iz) > Vd(ib,ig,is,iz) ) then
            epol(ib,ig,is,iz) = 1.0D0
        endif
        
    enddo    
    enddo    
    enddo    
    enddo
	
	!---computed expected default
    Ednp = 0.0D0
	do ig = 1, Ng
	do is = 1, Ns
		do ig_np = 1, Ng
		do is_np = 1, Ns
		do iz_np = 1, Nz
			Ednp(:,ig,is) = Ednp(:,ig,is) + dpol(:,ig_np,is_np,iz_np)*Pz(iz_np)*Pg(ig,ig_np)*Ps(is,is_np)
		enddo
		enddo
		enddo		
	enddo
	enddo
	
    end subroutine COMPUTE_Vd
    
    
    !-------------------------------------------------------
    !--- 2. compute expected Wd 
    !-------------------------------------------------------
    subroutine COMPUTE_EWv
    
    implicit none
    integer :: ig, is, ig_np, is_np, iz_np
     
    EWv = 0.0D0; EWtv = 0.0D0; EVd = 0.0D0;
    do ig = 1, Ng
    do is = 1, Ns
        do ig_np = 1, Ng
        do is_np = 1, Ns
        do iz_np = 1, Nz    
            EWv(:,ig,is)  = EWv(:,ig,is)  + Wv(:,ig_np,is_np,iz_np) *Pz(iz_np)*Pg(ig,ig_np)*Ps(is,is_np)
            EVd(:,ig,is)  = EVd(:,ig,is)  + Vd(:,ig_np,is_np,iz_np) *Pz(iz_np)*Pg(ig,ig_np)*Ps(is,is_np)     
			EWtv(:,ig,is) = EWtv(:,ig,is) + Wtv(:,ig_np,is_np,iz_np)*Pz(iz_np)*Pg(ig,ig_np)*Ps(is,is_np)     
        enddo        
        enddo        
        enddo        
    enddo    
    enddo
    
    end subroutine COMPUTE_EWv
    
    
    !----------------------------------
    !--- 3. compute Qprice and Xprice
    !----------------------------------
    subroutine COMPUTE_QandX
    
    implicit none
    real(8) :: gg, btog(Nb), bx(Nb), bkpa(Nb), Xfi(Nb,Ng,Ns,Nz), omexXfi(Nb,Ng,Ns,Nz), exQfi(Nb,Ng,Ns,Nz), omexX(Nb), Qkpa(Nb), exQ(Nb), EQ_np(Nb), EX_np(Nb), ax(Nb), Qnew(Nb,Ng,Ns,Nz), Xnew(Nb,Ng,Ns,Nz)
    integer :: ig, is, iz, ig_np, is_np, iz_np, btog_ind(Nb), btogdm_ind(Nb)
    
    real(8) :: errQ(1), errX(1), tolQandX, wwQ, wwX
    integer :: itQandX, maxitQandX, showQnadX
	integer :: QXclck_counts_beg, QXclck_counts_end, QXclck_rate
    
	maxitQandX = 500; tolQandX = 1e-6; wwQ = 0.40D0; wwX = 0.40D0; showQnadX = 0;

	QQub = dlta*Rs/(Rs-1.0d0+dlta); QQlb = 1.0e-7;
	XXub = (((1-pthta)*kppa)/(Rs-pthta))*QQub; XXlb = 1.0e-7;
			
    do itQandX = 1, maxitQandX
        
        !---compute objects for integration
        do ig = 1, Ng
        do is = 1, Ns
        do iz = 1, Nz
            
            !--- (b/g)*((1-m*dlta) index
            gg = gvec(ig);  btog = bvec/gg         ! btog_ind = btog_vind(:,ig) 		! btog_ind = vsearch_FET(btog,bvec);   	
			bx = bvec/gg			
			btogdm_ind = btogdm_vind(:,ig)  ! btogdm_vind(:,ig) = vsearch_FET(bx,bvec); 
						
			!---------------------------------------------------------------------------------------
	        !---T1: evluate Xprice at b/g
            Xfi(:,ig,is,iz) = Xprice(btogdm_ind,ig,is,iz) + (bx - bvec(btogdm_ind)) * ( Xprice(btogdm_ind+1,ig,is,iz) - Xprice(btogdm_ind,ig,is,iz) )/(bvec(btogdm_ind+1) - bvec(btogdm_ind))
            Xfi(:,ig,is,iz) = max(Xfi(:,ig,is,iz),XXlb)
            Xfi(:,ig,is,iz) = min(Xfi(:,ig,is,iz),XXub)            
			
            !---T2: (1-e(b,ig,is,iz))*X(b,ig,is,iz) -> then evaluate at b/g
            omexX   = (1-epol(:,ig,is,iz)) * Xprice(:,ig,is,iz)
            omexXfi(:,ig,is,iz) = omexX(btogdm_ind) + (bx - bvec(btogdm_ind)) * ( omexX(btogdm_ind+1) - omexX(btogdm_ind) )/(bvec(btogdm_ind+1) - bvec(btogdm_ind))
            
            !---T3: compute Q(kpa b,ig,is,iz), then e(b,ig,is,iz)*kpa*Q(kpa b,ig,is,iz) -> then evaluate at b/g
            bkpa = kppa*bvec; ! bkpa_ind = vsearch_FET(bkpa,bvec);   
            Qkpa = Qprice(bkpa_ind,ig,is,iz) + (bkpa - bvec(bkpa_ind)) * ( Qprice(bkpa_ind+1,ig,is,iz) - Qprice(bkpa_ind,ig,is,iz) )/(bvec(bkpa_ind+1) - bvec(bkpa_ind))
            Qkpa = max(Qkpa,QQlb)
            Qkpa = min(Qkpa,QQub)
			            
            exQ   = epol(:,ig,is,iz) * kppa * Qkpa
            exQfi(:,ig,is,iz) = exQ(btogdm_ind) + (bx - bvec(btogdm_ind)) * ( exQ(btogdm_ind+1) - exQ(btogdm_ind) )/(bvec(btogdm_ind+1) - bvec(btogdm_ind))
			!---------------------------------------------------------------------------------------
			
        enddo
        enddo
        enddo    
        
        !---compute integration
        do ig = 1, Ng
        do is = 1, Ns
        do iz = 1, Nz
            
            EQ_np = 0.0D0 
            EX_np = 0.0D0 
            
            do ig_np = 1, Ng
            do is_np = 1, Ns
            do iz_np = 1, Nz
                EQ_np  = EQ_np  + Qprice(bopt_loc(:,ig,is,iz),ig_np,is_np,iz_np) *Pz(iz_np)*Pg(ig,ig_np)*Ps(is,is_np)                                
                
				ax = pthta*Xfi(:,ig_np,is_np,iz_np) + (1.0D0-pthta)*( omexXfi(:,ig_np,is_np,iz_np) +  exQfi(:,ig_np,is_np,iz_np) )
				EX_np = EX_np + ax*Pz(iz_np)*Pg(ig,ig_np)*Ps(is,is_np)                                 
            enddo            
            enddo            
            enddo
            
            Qnew(:,ig,is,iz) = (1.0D0-dpol(:,ig,is,iz)) * ( dlta + ((1.0D0-dlta)/Rs)*EQ_np )  + dpol(:,ig,is,iz)*Xprice(:,ig,is,iz)
            
            Xnew(:,ig,is,iz) = (1.0D0/Rs) * EX_np
            
        enddo        
        enddo        
        enddo
        
        errQ = maxval(abs(Qnew - Qprice))
        errX = maxval(abs(Xnew - Xprice))
        
        if (showQnadX <= 0) then
            showQnadX = 20+1
        endif
        showQnadX = showQnadX-1
        
        if ( (errQ(1)<tolQandX) .and. (errX(1)<tolQandX) ) then
            exit
        else
            Qprice = wwQ*Qprice + (1.0D0-wwQ)*Qnew
            Xprice = wwX*Xprice + (1.0D0-wwX)*Xnew
        endif        
        
    enddo  ! end itQandX   	
    
    !---compue EQ = sum_{}
    EQ = 0.0D0
    EXX = 0.0D0	
    do ig = 1, Ng
    do is = 1, Ns    
        do ig_np = 1, Ng
        do is_np = 1, Ns
        do iz_np = 1, Nz
            EQ(:,ig,is) = EQ(:,ig,is) + Qprice(:,ig_np,is_np,iz_np)*Pz(iz_np)*Pg(ig,ig_np)*Ps(is,is_np)   
            EXX(:,ig,is) = EXX(:,ig,is) + Xprice(:,ig_np,is_np,iz_np)*Pz(iz_np)*Pg(ig,ig_np)*Ps(is,is_np)   			             
        enddo        
        enddo        
        enddo    
    enddo    
    enddo
    
    end subroutine COMPUTE_QandX
    
    !-------------------------------------------------------
    !--- 4. compute R schedule
    !-------------------------------------------------------
    subroutine COMPUTE_Rschedule
    
    implicit none
    real(8) :: gg
    integer :: ig, is, ib, ibp
	real(8) :: aux(Nb),aux1(Nb),aux2(Nb)

	lsch  = 0
	isch  = 0
	lisch = 0
	msch = 1    

	do ig = 1, Ng
	do is = 1, Ns
		gg = gvec(ig)
		Rsched(:,ig,is)   = Rs/EQ(:,ig,is)
		rhosched(:,ig,is) =  dlta*Rsched(:,ig,is) - dlta
		do ib = 1, Nb
			nsched(:,ib,ig,is) = (gg*bvec - (1-dlta)*bvec(ib) )/(dlta*Rsched(:,ig,is))
			if ((wsunspot==1).AND.(is==1)) then				
				where (Ednp(:,ig,is)<0.999*pdub) lsch(:,ib,ig,is)=1 
				! identifying increasing parts
				!$OMP PARALLEL DO
				do ibp=(ws+1),(Nb-ws)
					aux1(ibp) = sum(nsched(ibp-ws:ibp,ib,ig,is))/(ws+1)
			 		aux2(ibp) = sum(nsched(ibp:ibp+ws,ib,ig,is))/(ws+1)
			 		if (aux1(ibp)<aux2(ibp)) isch(ibp,ib,ig,is)=1
			 	enddo
				!$OMP END PARALLEL DO
				isch(1:ws,ib,ig,is)=isch((ws+1),ib,ig,is)
				isch((Nb-(ws-1)):Nb,ib,ig,is)=isch(Nb-ws,ib,ig,is)
				! lisch (increasing and below the limit)
				lisch(:,ib,ig,is) = lsch(:,ib,ig,is)+isch(:,ib,ig,is)
				!identifying multiplicity: there must be a lower nb, in the increasing part, with respective R below 0.99*Rmax
				!if msch = 0, it is not selected when bad sunspot occurs
				!$OMP PARALLEL DO            
				do ibp = 1,Nb
					aux(ibp) = minval(nsched(ibp:Nb,ib,ig,is),MASK=(lisch(ibp:Nb,ib,ig,is)==2));
					if (aux(ibp)<nsched(ibp,ib,ig,is)) msch(ibp,ib,ig,is)=0
				enddo
				!$OMP END PARALLEL DO
			end if !((wsunspot==1).AND.(is==1))
		enddo ! end ib   
	enddo  ! end is    
	enddo  ! end ig    
    
    end subroutine COMPUTE_Rschedule
    
    !-------------------------------------------------------
    !--- 5. Vnd eval and max
    !-------------------------------------------------------
    function Vnd_eval(ib,ig,is,iz) result(V_out)
    
    integer :: ib, ig, is, iz, ibmax(1)    
    real(8) :: V_out
    
    real(8) :: gg, rx(Nb), px(Nb), cxb(Nb), uxb(Nb), V_np(Nb), Vxb(Nb)
    
    gg   = gvec(ig); rx = rhosched(:,ig,is); px = Ednp(:,ig,is);
    
    cxb  = Y(ig,iz) + nsched(:,ib,ig,is) - bvec(ib)
    cxb  = max(cxb, cmin)
    uxb  = (cxb**(1.0D0-gma))/(1.0D0-gma)
    
    V_np = EWv(:,ig,is)
    
    Vxb  = uxb + bta * (gg**(1.0D0-gma)) * V_np	
	where (px>pdub) Vxb = -99999999999999999999999.D0		
	if ((wsunspot==1).AND.(is==1)) then
		where (msch(:,ib,ig,is)==0) Vxb  =  -99999999999999999999999.D0; !Vd(ib) -1.0e-04;
	endif
    
    ibmax = maxloc(Vxb)
    
    bopt_loc(ib,ig,is,iz) = ibmax(1)
	bpol(ib,ig,is,iz)     = bvec(bopt_loc(ib,ig,is,iz))
	rpol(ib,ig,is,iz)     = rx(bopt_loc(ib,ig,is,iz))
	ppol(ib,ig,is,iz)     = px(bopt_loc(ib,ig,is,iz))	
    V_out                 = Vxb(ibmax(1))
    
    end function Vnd_eval

end module FUNCTIONS