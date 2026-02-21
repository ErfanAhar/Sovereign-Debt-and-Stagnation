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
    integer :: ib, ig, iz, ig_np, iz_np, bind_np(Nb) 
    real(8) :: Vdnew(Nb,Ng,Nz), errVd(1), tolVd, wwVd
    integer :: itVd, maxitVd
    
	maxitVd = 500; tolVd = 1e-8; wwVd = 0.25D0
	
	!---compute Vnd(kppa(b,g,s,z)
	bkpa = kppa*bvec   ! bkpa_ind = vsearch_FET(bkpa,bvec);   
	do ig = 1, Ng	
	do iz = 1, Nz
		Vndk(:,ig,iz) = Vnd(bkpa_ind,ig,iz) + (bkpa - bvec(bkpa_ind)) * ( Vnd(bkpa_ind+1,ig,iz) - Vnd(bkpa_ind,ig,iz) )/(bvec(bkpa_ind+1) - bvec(bkpa_ind))
	enddo	
	enddo
    
    
    do itVd = 1, maxitVd
        
        Wv  = max(Vnd , Vd)	
		Wtv = max(Vndk, Vd)    
		
        call COMPUTE_EWv
        
        do ig = 1, Ng                  
            
			gg  = gvec(ig)
			bnp = bvec/gg
			bind_np = btogdm_vind(:,ig)   ! = vsearch_FET(bx,bvec); with bx = bvec/gg
			
            Vd_np  = EVd(bind_np,ig)  + (bnp - bvec(bind_np)) * ( EVd(bind_np+1,ig)  - EVd(bind_np,ig)  )/(bvec(bind_np+1) - bvec(bind_np))
            Wtv_np = EWtv(bind_np,ig) + (bnp - bvec(bind_np)) * ( EWtv(bind_np+1,ig) - EWtv(bind_np,ig) )/(bvec(bind_np+1) - bvec(bind_np))			
			
            Vx_np = bta * (gg**(1.0D0-gma)) * (pthta*Vd_np + (1.0D0-pthta)*Wtv_np)                        
            
            do iz = 1, Nz
                cxb = ybarg(ig)*Y(ig,iz)
                cxb = max(cxb, cmin)
                uxb = (cxb**(1.0D0-gma))/(1.0D0-gma)
                
                Vdnew(:,ig,iz) = uxb + Vx_np
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
    do iz = 1, Nz
        if (Vd(ib,ig,iz) > Vnd(ib,ig,iz) ) then
            dpol(ib,ig,iz) = 1.0D0
        endif
        
        if (Vndk(ib,ig,iz) > Vd(ib,ig,iz) ) then
            epol(ib,ig,iz) = 1.0D0
        endif
        
    enddo    
    enddo        
    enddo
	
	!---computed expected default
    Ednp = 0.0D0
	do ig = 1, Ng	
		do ig_np = 1, Ng		
		do iz_np = 1, Nz
			Ednp(:,ig) = Ednp(:,ig) + dpol(:,ig_np,iz_np)*Pz(iz_np)*Pg(ig,ig_np)
		enddo		
		enddo			
	enddo
	
    end subroutine COMPUTE_Vd

    !-------------------------------------------------------
    !--- 2. compute expected Wd 
    !-------------------------------------------------------
    subroutine COMPUTE_EWv
    
    implicit none
    integer :: ig, ig_np, iz_np
     
    EWv = 0.0D0; EWtv = 0.0D0; EVd = 0.0D0;
    do ig = 1, Ng    
        do ig_np = 1, Ng        
        do iz_np = 1, Nz    
            EWv(:,ig)  = EWv(:,ig)  + Wv(:,ig_np,iz_np) *Pz(iz_np)*Pg(ig,ig_np)
            EVd(:,ig)  = EVd(:,ig)  + Vd(:,ig_np,iz_np) *Pz(iz_np)*Pg(ig,ig_np)
			EWtv(:,ig) = EWtv(:,ig) + Wtv(:,ig_np,iz_np)*Pz(iz_np)*Pg(ig,ig_np) 
        enddo                
        enddo            
    enddo
    
    end subroutine COMPUTE_EWv

    !----------------------------------
    !--- 3. compute Qprice and Xprice
    !----------------------------------
    subroutine COMPUTE_QandX
    
    implicit none
    real(8) :: gg, btog(Nb), bx(Nb), bkpa(Nb), Xfi(Nb,Ng,Nz), omexXfi(Nb,Ng,Nz), exQfi(Nb,Ng,Nz), omexX(Nb), Qkpa(Nb), exQ(Nb), EX_np(Nb), ax(Nb), Xnew(Nb,Ng,Nz)
    integer :: ig, iz, ig_np, iz_np, btog_ind(Nb), btogdm_ind(Nb)
    
    real(8) :: errQ(1), errX(1), tolQandX, wwQ, wwX
    integer :: itQandX, maxitQandX, showQnadX
	integer :: QXclck_counts_beg, QXclck_counts_end, QXclck_rate
    
	maxitQandX = 500; tolQandX = 1e-6; wwQ = 0.40D0; wwX = 0.40D0; showQnadX = 0;

	QQub = Qub; QQlb = 1.0e-7; !QQub = dlta*Rs/(Rs-1.0d0+dlta); QQlb = 1.0e-7;
	XXub = Xub; XXlb = 1.0e-7; !XXub = (((1-pthta)*kppa)/(Rs-pthta))*QQub; XXlb = 1.0e-7;
			
    do itQandX = 1, maxitQandX

        !--update Q given X
        Qprice = (1.0D0 - dpol) + dpol*Xprice
        
        !---compute objects for integration
        do ig = 1, Ng        
        do iz = 1, Nz            
            
            !--- (b/g)*((1-m*dlta) index
            gg = gvec(ig);  btog = bvec/gg         ! btog_ind = btog_vind(:,ig) 		! btog_ind = vsearch_FET(btog,bvec);   	
			bx = bvec/gg			
			btogdm_ind = btogdm_vind(:,ig)  ! btogdm_vind(:,ig) = vsearch_FET(bx,bvec); 
						
			!---------------------------------------------------------------------------------------
	        !---T1: evluate Xprice at b/g
            Xfi(:,ig,iz) = Xprice(btogdm_ind,ig,iz) + (bx - bvec(btogdm_ind)) * ( Xprice(btogdm_ind+1,ig,iz) - Xprice(btogdm_ind,ig,iz) )/(bvec(btogdm_ind+1) - bvec(btogdm_ind))
            Xfi(:,ig,iz) = max(Xfi(:,ig,iz),XXlb)
            Xfi(:,ig,iz) = min(Xfi(:,ig,iz),XXub)            
			
            !---T2: (1-e(b,ig,is,iz))*X(b,ig,is,iz) -> then evaluate at b/g
            omexX   = (1-epol(:,ig,iz)) * Xprice(:,ig,iz)
            omexXfi(:,ig,iz) = omexX(btogdm_ind) + (bx - bvec(btogdm_ind)) * ( omexX(btogdm_ind+1) - omexX(btogdm_ind) )/(bvec(btogdm_ind+1) - bvec(btogdm_ind))
            
            !---T3: compute Q(kpa b,ig,is,iz), then e(b,ig,is,iz)*kpa*Q(kpa b,ig,is,iz) -> then evaluate at b/g
            bkpa = kppa*bvec; ! bkpa_ind = vsearch_FET(bkpa,bvec);   
            Qkpa = Qprice(bkpa_ind,ig,iz) + (bkpa - bvec(bkpa_ind)) * ( Qprice(bkpa_ind+1,ig,iz) - Qprice(bkpa_ind,ig,iz) )/(bvec(bkpa_ind+1) - bvec(bkpa_ind))
            Qkpa = max(Qkpa,QQlb)
            Qkpa = min(Qkpa,QQub)
			            
            exQ   = epol(:,ig,iz) * kppa * Qkpa
            exQfi(:,ig,iz) = exQ(btogdm_ind) + (bx - bvec(btogdm_ind)) * ( exQ(btogdm_ind+1) - exQ(btogdm_ind) )/(bvec(btogdm_ind+1) - bvec(btogdm_ind))
			!---------------------------------------------------------------------------------------
			
        enddo
        enddo
          
        
        !---compute integration
        do ig = 1, Ng        
        do iz = 1, Nz            
            
            EX_np = 0.0D0 
            
            do ig_np = 1, Ng            
            do iz_np = 1, Nz                
                
				ax = pthta*Xfi(:,ig_np,iz_np) + (1.0D0-pthta)*( omexXfi(:,ig_np,iz_np) +  exQfi(:,ig_np,iz_np) )
				EX_np = EX_np + ax*Pz(iz_np)*Pg(ig,ig_np)
            enddo                                
            enddo
            
            !Qnew(:,ig,is,iz) = (1.0D0-dpol(:,ig,is,iz)) * ( dlta + ((1.0D0-dlta)/Rs)*EQ_np )  + dpol(:,ig,is,iz)*Xprice(:,ig,is,iz)
            
            Xnew(:,ig,iz) = (1.0D0/Rs) * EX_np
            
        enddo                
        enddo
        
        !errQ = maxval(abs(Qnew - Qprice))
        errX = maxval(abs(Xnew - Xprice))
        
        if (showQnadX <= 0) then
            showQnadX = 20+1
        endif
        showQnadX = showQnadX-1
        
        if ( errX(1)<tolQandX ) then     !if ( (errQ(1)<tolQandX) .and. (errX(1)<tolQandX) ) then
            exit
        else
            !Qprice = wwQ*Qprice + (1.0D0-wwQ)*Qnew
            Xprice = wwX*Xprice + (1.0D0-wwX)*Xnew
        endif        
        
    enddo  ! end itQandX   	
    
    !---compue EQ = sum_{}
    EQ = 0.0D0
    EXX = 0.0D0	
    do ig = 1, Ng    
        do ig_np = 1, Ng        
        do iz_np = 1, Nz
            EQ(:,ig)  = EQ(:,ig)  + Qprice(:,ig_np,iz_np)*Pz(iz_np)*Pg(ig,ig_np)
            EXX(:,ig) = EXX(:,ig) + Xprice(:,ig_np,iz_np)*Pz(iz_np)*Pg(ig,ig_np)
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
    integer :: ig, ib, ibp
	real(8) :: aux(Nb),aux1(Nb),aux2(Nb)

	lsch  = 0
	isch  = 0
	lisch = 0
	msch = 1    

	do ig = 1, Ng	
		gg = gvec(ig)
		Rsched(:,ig)   = Rs/EQ(:,ig)
        nsched(:,ig)   = (gg*bvec )/(Rsched(:,ig))
		
        !rhosched(:,ig,is) =  dlta*Rsched(:,ig,is) - dlta		
	enddo  ! end ig 
    
    rhosched = Rsched - 1.0D0
    
    end subroutine COMPUTE_Rschedule
    
    !-------------------------------------------------------
    !--- 5. Vnd eval and max
    !-------------------------------------------------------
    function Vnd_eval(ib,ig,iz) result(V_out)
    
    integer :: ib, ig, iz, ibmax(1)    
    real(8) :: V_out
    
    real(8) :: gg, rx(Nb), px(Nb), cxb(Nb), uxb(Nb), V_np(Nb), Vxb(Nb)
    
    gg   = gvec(ig); rx = rhosched(:,ig); px = Ednp(:,ig);
    
    cxb  = Y(ig,iz) + nsched(:,ig) - bvec(ib)
    cxb  = max(cxb, cmin)
    uxb  = (cxb**(1.0D0-gma))/(1.0D0-gma)
    
    V_np = EWv(:,ig)
    
    Vxb  = uxb + bta * (gg**(1.0D0-gma)) * V_np	

    where (px>pdub) Vxb = -99999999999999999999999.D0		
    
    ibmax = maxloc(Vxb)
    
    bopt_loc(ib,ig,iz) = ibmax(1)
	bpol(ib,ig,iz)     = bvec(bopt_loc(ib,ig,iz))
	rpol(ib,ig,iz)     = rx(bopt_loc(ib,ig,iz))
	ppol(ib,ig,iz)     = px(bopt_loc(ib,ig,iz))	
    V_out              = Vxb(ibmax(1))
    
    end function Vnd_eval




end module FUNCTIONS