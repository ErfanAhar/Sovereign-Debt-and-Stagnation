module ModuleSAVE

use mod_parameters

contains

    !------------------------------------!
    !--- save grids
    !------------------------------------!
    subroutine SaveGRIDS
    implicit none 
    real(8) :: prms(12)
    integer :: ib, ig, iz, ii, ig_np
    
    open(1,  file = 'OUTPUT_'//trim(wvS)//'/bvec.txt',    status = 'unknown')
    open(2,  file = 'OUTPUT_'//trim(wvS)//'/zvec.txt',    status = 'unknown')
    open(3,  file = 'OUTPUT_'//trim(wvS)//'/Pzvec.txt',   status = 'unknown')
    do ib = 1, Nb;     write(1, *) bvec(ib); enddo
    do iz = 1, Nz;     write(2, *) zvec(iz); enddo	
    do iz = 1, Nz;     write(3, *) Pz(iz)  ; enddo		
    close(1); close(2); close(3)
    
    open(1,  file = 'OUTPUT_'//trim(wvS)//'/Pg.txt',    status = 'unknown')
    do ig = 1, Ng
        write(1, '(*(f18.6))') ( Pg(ig,ig_np), ig_np= 1,Ng )	
    enddo
    close(1)
    
    ! open(1,  file = 'OUTPUT_'//trim(wvS)//'/Ps.txt',    status = 'unknown')
    ! do is = 1, Ns
    ! 	write(1, '(*(f18.6))') ( Ps(is,is_np), is_np= 1,Ns )	
    ! enddo
    ! close(1)	
    
    open(1,  file = 'OUTPUT_'//trim(wvS)//'/Ygz.txt',    status = 'unknown')
    do iz = 1, Nz
        write(1, '(*(f18.6))') ( Y(ig,iz), ig=1,Ng )
    enddo
    close(1)
    
    open(1,  file = 'OUTPUT_'//trim(wvS)//'/gvec.txt',    status = 'unknown')
    write(1, '(*(f18.6))') ( gvec(ig), ig= 1,Ng )
    close(1)
    
    open(1,  file = 'OUTPUT_'//trim(wvS)//'/prms.txt',    status = 'unknown')
    prms = [Rs,dlta,ybarL,ybarH,kppa,pgL,pgH,bta,gL,gH,sdz,psB]
    write(1, '(*(f18.6))') ( prms(ii), ii= 1,12 )
    close(1)

    end subroutine SaveGRIDS

    !--------------------------------------!
    !---   save value and policies      ---!
    !--------------------------------------!
    subroutine SaveVFandPolicies
    
    implicit none
    integer :: ib, ig, iz   
    

    !print *, '1) about to start saving VF and policies'    

    !---save values
    open(1,  file = 'OUTPUT_'//trim(wvS)//'/Vnd.txt'	  ,  status = 'unknown'); ! print *, 'open 1';
    open(2,  file = 'OUTPUT_'//trim(wvS)//'/Vd.txt' 	  ,  status = 'unknown'); ! print *, 'open 2';
    open(3,  file = 'OUTPUT_'//trim(wvS)//'/dpol.txt' 	  ,  status = 'unknown'); ! print *, 'open 3';
    open(4,  file = 'OUTPUT_'//trim(wvS)//'/epol.txt'     ,  status = 'unknown'); ! print *, 'open 4';
	open(5,  file = 'OUTPUT_'//trim(wvS)//'/bopt_loc.txt' ,  status = 'unknown'); ! print *, 'open 5';
	open(70,  file = 'OUTPUT_'//trim(wvS)//'/bpol.txt' ,  status = 'unknown'); ! print *, 'open 6'; 
	open(7,  file = 'OUTPUT_'//trim(wvS)//'/rpol.txt' ,  status = 'unknown'); ! print *, 'open 7';
	open(8,  file = 'OUTPUT_'//trim(wvS)//'/ppol.txt' ,  status = 'unknown'); ! print *, 'open 8';

    !print *, '1.1) still here'    
	
    do ig = 1, Ng            
    do iz = 1, Nz
        do ib = 1, Nb
            write(1, '(*(f18.6))') ( max( Vnd(ib,ig,iz) , Vlow) )
            write(2, '(*(f18.6))') ( max( Vd(ib,ig,iz)  , Vlow) )
            write(3, '(*(f18.6))') ( dpol(ib,ig,iz)  )
            write(4, '(*(f18.6))') ( epol(ib,ig,iz)  )
			write(5, '(I8.1)')     ( bopt_loc(ib,ig,iz)  )
            write(70, '(*(f18.6))') ( bpol(ib,ig,iz)  )
            write(7, '(*(f18.6))') ( rpol(ib,ig,iz)  )
            write(8, '(*(f18.6))') ( ppol(ib,ig,iz)  )
        enddo        
    enddo        
    enddo

    !print *, '1.x) I made it here '

    close(1); close(2); close(3); close(4); close(5); close(70); close(7); close(8)   ; 
    
    ! print *, '2) starting with Ws'
	
	open(1,  file = 'OUTPUT_'//trim(wvS)//'/Vndk.txt'	 ,  status = 'unknown')
    open(2,  file = 'OUTPUT_'//trim(wvS)//'/Wv.txt' 	 ,  status = 'unknown')
	open(3,  file = 'OUTPUT_'//trim(wvS)//'/Wtv.txt' 	 ,  status = 'unknown')
	do ig = 1, Ng    
    do iz = 1, Nz
		do ib = 1, Nb
			write(1, '(*(f18.6))') (  Vndk(ib,ig,iz)  )
			write(2, '(*(f18.6))') (  Wv(ib,ig,iz)  )
			write(3, '(*(f18.6))') (  Wtv(ib,ig,iz)  )
		enddo
	enddo	
	enddo
	close(1); close(2); close(3);

    !print *, '3) to do prices'
	
	open(1,  file = 'OUTPUT_'//trim(wvS)//'/EQ.txt'	 ,  status = 'unknown')
    open(2,  file = 'OUTPUT_'//trim(wvS)//'/EXX.txt' 	 ,  status = 'unknown')
    do ig = 1, Ng        
		do ib = 1, Nb
			write(1, '(*(f18.6))') (  EQ(ib,ig)  )
			write(2, '(*(f18.6))') (  EXX(ib,ig)  )
        enddo            
    enddo
	
    !print *, '4) to do expected values'
	open(1,  file = 'OUTPUT_'//trim(wvS)//'/EWv.txt'	 ,  status = 'unknown')
	open(2,  file = 'OUTPUT_'//trim(wvS)//'/EWtv.txt'	 ,  status = 'unknown')
    open(3,  file = 'OUTPUT_'//trim(wvS)//'/EVd.txt' 	 ,  status = 'unknown')
	open(4,  file = 'OUTPUT_'//trim(wvS)//'/Ednp.txt' 	 ,  status = 'unknown')	
	do ig = 1, Ng    
		do ib = 1, Nb
			write(1, '(*(f18.6))') (  EWv(ib,ig)  )
			write(2, '(*(f18.6))') (  EWtv(ib,ig) )
			write(3, '(*(f18.6))') (  EVd(ib,ig)  )
			write(4, '(*(f18.6))') (  Ednp(ib,ig) )
		enddo    
	enddo
	close(1); close(2); close(3); close(4);
	
	
    !print *, 'about to save initial guess'
	!---save for intial guess
    open(1,  file = 'OUTPUT_'//trim(wvS)//'/Vnd_in.txt',        status = 'unknown')
    open(2,  file = 'OUTPUT_'//trim(wvS)//'/bopt_loc_in.txt' ,  status = 'unknown')
    open(3,  file = 'OUTPUT_'//trim(wvS)//'/Vd_in.txt' ,  status = 'unknown')	
    open(4,  file = 'OUTPUT_'//trim(wvS)//'/Qprice_in.txt' ,  status = 'unknown')	
    open(5,  file = 'OUTPUT_'//trim(wvS)//'/Xprice_in.txt' ,  status = 'unknown')
	open(70,  file = 'OUTPUT_'//trim(wvS)//'/Wv_in.txt' ,  status = 'unknown')			
    write(1, *) Vnd
    write(2, *) bopt_loc
    write(3, *) Vd
    write(4, *) Qprice
    write(5, *) Xprice
    write(70, *) Wv	
    close(1); close(2); close(3); close(4); close(5); close(70)   

    ! print *, 'at the end of SaveVFandPolicies'
    
    end subroutine SaveVFandPolicies
    
    
    !--------------------------------------!
    !---   save value and policies      ---!
    !--------------------------------------!
    subroutine SaveSCHEDULE
    
    implicit none
    integer :: ib, ib_np, ig, iz 
    
    open(1,  file = 'OUTPUT_'//trim(wvS)//'/Rsched.txt',   status = 'unknown')
    do ig = 1, Ng    
        do ib_np = 1, Nb
			write(1, *) ( min( Rsched(ib_np,ig) , 10000.01D0 ) )        
        enddo    
    enddo        
    close(1)
    
	open(1,  file = 'OUTPUT_'//trim(wvS)//'/Qprice.txt',   status = 'unknown')
	open(2,  file = 'OUTPUT_'//trim(wvS)//'/Xprice.txt',   status = 'unknown')
	do ig = 1, Ng    
    do ib = 1, Nb 
	do iz = 1, Nz
		write(1, '(*(f18.6))') ( Qprice(ib,ig,iz) )
		write(2, '(*(f18.6))') ( Xprice(ib,ig,iz) )
	enddo
	enddo	
	enddo
	close(1); close(2);
	
    end subroutine SaveSCHEDULE


end module ModuleSAVE