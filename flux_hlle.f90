subroutine flux_hlle

  use arrays
  use global_numbers

  implicit none

  integer i,j

  real(kind=8) idx

  real(kind=8) DL,DR
  real(kind=8) SxL,SxR
  real(kind=8) TauL,TauR
  real(kind=8) u4L,u4R
  real(kind=8) u5L,u5R
  real(kind=8) pressL,pressR
  real(kind=8) denL,denR
  real(kind=8) vxL,vxR
  real(kind=8) CSL,CSR
  real(kind=8) ErL,ErR
  real(kind=8) FrL,FrR
  real(kind=8) PrL,PrR
  real(kind=8) f_L,f_R
  real(kind=8) fcomovil_L,fcomovil_R
  real(kind=8) eeL,eeR
  real(kind=8) hL,hR
  real(kind=8) WL,WR
  real(kind=8) VelL,VelR
  real(kind=8) lamb1L,lamb1R
  real(kind=8) lamb2L,lamb2R
  real(kind=8) lamb3L,lamb3R
  real(kind=8) lamb4L,lamb4R
  real(kind=8) lamb5L,lamb5R
  real(kind=8) lamb6L,lamb6R
  real(kind=8) lamb_max,lamb_min
  real(kind=8) f1L,f1R
  real(kind=8) f2L,f2R
  real(kind=8) f3L,f3R 
  real(kind=8) f4L,f4R
  real(kind=8) f5L,f5R  
  real(kind=8) reconstructor
  real(kind=8) sigma, Cr
  real(kind=8) zetaL,zetaR
  real(kind=8) dzetaL,dzetaR
  real(Kind=8) thick_PL,thick_PR
  real(Kind=8) thick_ML,thick_MR
  real(kind=8) thin_PL,thin_PR
  real(kind=8) thin_ML,thin_MR
  real(kind=8) limL,limR
  real(kind=8) ReL,ReR,Re1L,Re1R
  real(kind=8) partialPr_ErL,partialPr_ErR
  real(kind=8) partialPr_FrL,partialPr_FrR
  real(kind=8) ErComovilL,ErComovilR
  real(kind=8) FrComovilL,FrComovilR
  real(kind=8) chi_L,chi_R
  real(kind=8) vel_advL,vel_advR
  real(kind=8) limxComovilL,limxComovilR

   idx  = 1.0D0/dx

  do i = -g_pts,Nx+g_pts-1
   do j=1,ngroups
     if (Limiter.eq.'minmod_TVD' .or. Limiter.eq.'MC_TVD') then 

        call lpm(i,press(i-1),press(i),press(i+1),press(i+2),pressL,pressR)
        call lpm(i,rho(i-1),rho(i),rho(i+1),rho(i+2),denL,denR)
        call lpm(i,vx(i-1),vx(i),vx(i+1),vx(i+2),vxL,vxR)
        call lpm(i,CS(i-1),CS(i),CS(i+1),CS(i+2),CSL,CSR)
        call lpm(i,Er(i-1,j),Er(i,j),Er(i+1,j),Er(i+2,j),ErL,ErR)
        call lpm(i,Fr(i-1,j),Fr(i,j),Fr(i+1,j),Fr(i+2,j),FrL,FrR)
        call lpm(i,f(i-1,j),f(i,j),f(i+1,j),f(i+2,j),f_L,f_R)
        call lpm(i,Pr(i-1,j),Pr(i,j),Pr(i+1,j),Pr(i+2,j),PrL,PrR)
        call lpm(i,ErComovil(i-1,j),ErComovil(i,j),ErComovil(i+1,j),ErComovil(i+2,j),ErComovilL,ErComovilR)
        call lpm(i,FrComovil(i-1,j),FrComovil(i,j),FrComovil(i+1,j),FrComovil(i+2,j),FrComovilL,FrComovilR)
        call lpm(i,chi(i-1,j),chi(i,j),chi(i+1,j),chi(i+2,j),chi_L,chi_R)

     end if

     ! ============================================================ 

     if ( (vxL**2) .ge. (1.0d0 - 1.e-6) ) then
        pressL = press(i)
        denL   = rho(i)
        vxL    = vx(i)       
        CSL    = CS(i)       
        ErL    = Er(i,j)       
        FrL    = Fr(i,j)       
        f_L    = f(i,j)       
        PrL    = Pr(i,j)       
        ErComovilL = ErComovil(i,j)       
        FrComovilL = FrComovil(i,j)
        chi_L = chi(i,j)        
     else
     end if

     if ( (vxR**2) .ge. (1.0d0 - 1.e-6) ) then
        pressR = press(i+1)
        denR   = rho(i+1)
        vxR    = vx(i+1)       
        CSR    = CS(i+1)       
        ErR    = Er(i+1,j)       
        FrR    = Fr(i+1,j)       
        f_R    = f(i+1,j)       
        PrR    = Pr(i+1,j)      
        ErComovilR = ErComovil(i+1,j)       
        FrComovilR = FrComovil(i+1,j)  
        chi_R = chi(i+1,j)          
     else
     end if
 end do
     ! ============================================================ 

if(eos_type.eq.'id') then
     eeL = pressL/((gamma - one)*denL)
     eeR = pressR/((gamma - one)*denR)

     hL = one + eeL + pressL/denL 
     hR = one + eeR + pressR/denR 

     CSL = sqrt(gamma*pressL/hL/denL)
     CSR = sqrt(gamma*pressR/hR/denR)
else if(eos_type.eq.'tm') then
     hL = five*half*(pressL/denL) +sqrt(9.0d0/four*(pressL/denL)**2 + one)
     hR = five*half*(pressR/denR) +sqrt(9.0d0/four*(pressR/denR)**2 + one)
  
     CSL = sqrt(third*(pressL/denL)/hL*(five*hL -8.0d0*pressL/denL)/(hL -pressL/denL))
     CSR = sqrt(third*(pressR/denR)/hR*(five*hR -8.0d0*pressR/denR)/(hR -pressR/denR))

!     CSL = sqrt( (five* (pressL/denL)*sqrt((pressL/denL)**2 + four/9.0d0) +three*(pressL/denL)**2)/&
!            (12.0d0*(pressL/denL)*sqrt((pressL/denL)**2 + four/9.0d0) +12.0d0*(pressL/denL)**2 + two) )
!     CSR = sqrt( (five* (pressR/denR)*sqrt((pressR/denR)**2 + four/9.0d0) +three*(pressR/denR)**2)/&
!            (12.0d0*(pressR/denR)*sqrt((pressR/denR)**2 + four/9.0d0) +12.0d0*(pressR/denR)**2 + two) )
end if

     VelL = sqrt(vxL**2)
     VelR = sqrt(vxR**2)

     WL   = one/sqrt(one - VelL**2)
     WR   = one/sqrt(one - VelR**2)
 


     DL = denL*WL     
     DR = denR*WR     

     SxL = denL*hL*WL**2*velL
     SxR = denR*hR*WR**2*velR

     TauL = denL*hL*WL**2 - pressL
     TauR = denR*hR*WR**2 - pressR 

     u4L = ErL
     u4R = ErR

     u5L = FrL
     u5R = FrR

     ! ***************************************************************************
     ! ***************************************************************************
     ! ***************************************************************************

     
  fcomovil_L = max((-vxL*ERL + (one + vxL**2)*FRL - vxL*PrL)/((ERL - two*vxL*FRL + vxL**2*PrL)),floorrad)
  
   if (fcomovil_L.gt.one)then
    limxComovilL = one/fcomovil_L
         fcomovil_L = max((-vxL*ERL + (one + vxL**2)*FRL - vxL*PrL)*limxComovilL/((ERL - two*vxL*FRL + vxL**2*PrL)),floorrad)
    end if

  fcomovil_R = max((-vxR*ERR + (one + vxR**2)*FRR - vxR*PrR)/((ERR - two*vxR*FRR + vxR**2*PrR)),floorrad)
  
  if (fcomovil_R.gt.one)then
   limxComovilR = one/fcomovil_R
        fcomovil_R = max((-vxL*ERR + (one + vxR**2)*FRR - vxR*PrR)*limxComovilR/((ERR - two*vxR*FRR + vxR**2*PrR)),floorrad)
   end if
   
     Cr = sqrt(1.0d0/3.0d0)
     thick_PL = (two*WL*vxL + sqrt(two*WL**2+one-two*VelL**2))/(two*WL**2+one)
     thick_PR = (two*WR*vxR + sqrt(two*WR**2+one-two*VelR**2))/(two*WR**2+one)

     thick_ML = (two*WL*vxL - sqrt(two*WL**2+one-two*VelL**2))/(two*WL**2+one)
     thick_MR = (two*WR*vxR - sqrt(two*WR**2+one-two*VelR**2))/(two*WR**2+one)

     thin_PL = abs(FrL)/sqrt(FrL*FrL)
     thin_PR = abs(FrR)/sqrt(FrR*FrR)

     thin_ML = -abs(FrL)/sqrt(FrL*FrL)
     thin_MR = -abs(FrR)/sqrt(FrR*FrR)
	 

     ! Eigenvalues along x-direction

     ! triply degenerate eigenvalue

     lamb1L = vxL 
     lamb1R = vxR 

     ! the other two eigenvalues are : 

     lamb2L = (one/(one-VelL**2*CSL**2))*(vxL*(one - CSL**2) &
          + sqrt(CSL**2*(one-VelL**2)*((one - VelL**2*CSL**2) &
          - vxL*vxL*(one - CSL**2))))

     lamb2R = (one/(one-VelR**2*CSR**2))*(vxR*(one - CSR**2) &
          + sqrt(CSR**2*(one-VelR**2)*((one - VelR**2*CSR**2) &
          - vxR*vxR*(one - CSR**2))))

     lamb3L = (one/(one-VelL**2*CSL**2))*(vxL*(one - CSL**2) &
          - sqrt(CSL**2*(one-VelL**2)*((one - VelL**2*CSL**2) &
          - vxL*vxL*(one - CSL**2))))

     lamb3R = (one/(one-VelR**2*CSR**2))*(vxR*(one - CSR**2) &
          - sqrt(CSR**2*(one-VelR**2)*((one - VelR**2*CSR**2) &
          - vxR*vxR*(one - CSR**2))))
  if(close_type.eq.'dif')then
     lamb4L = thick_ML  
     lamb4R = thick_MR  

     lamb5L = thick_PL 
     lamb5R = thick_PR 
	 
	      lamb4L = max(lamb4L,cr)
	      lamb4R = max(lamb4R,cr)

	      lamb5L = min(lamb5L,-cr)
	      lamb5R = min(lamb5R,-cr)
  else if(close_type.eq.'M1')then


!takahashi & Ohsuga 2013 
	 
	 zetaL = (three + four*(fcomovil_L)**2)/(five + two*sqrt(four -three*(fcomovil_L)**2))
     zetaR = (three + four*(fcomovil_R)**2)/(five + two*sqrt(four -three*(fcomovil_R)**2))

 	 
	dzetaL = two*(fcomovil_L)/sqrt(four - three*(fcomovil_L)**2)
	dzetaR = two*(fcomovil_R)/sqrt(four - three*(fcomovil_R)**2)
	
    lamb4L = half*(dzetaL + sqrt(dzetaL**2 + four*zetaL - four*(fcomovil_L)*dzetaL))
    lamb4R = half*(dzetaR + sqrt(dzetaR**2 + four*zetaR - four*(fcomovil_R)*dzetaR))

    lamb5L = half*(dzetaL - sqrt(dzetaL**2 + four*zetaL - four*(fcomovil_L)*dzetaL))
    lamb5R = half*(dzetaR - sqrt(dzetaR**2 + four*zetaR - four*(fcomovil_R)*dzetaR))

	 lamb4L = min(lamb4L,third*four*idx/sqrt(chi_L*chi_L ))
	 lamb4R = min(lamb4R,third*four*idx/sqrt(chi_R*chi_R ))
	 
	 lamb5L = max(lamb5L,-third*four*idx/sqrt(chi_L*chi_L ))
	 lamb5R = max(lamb5R,-third*four*idx/sqrt(chi_R*chi_R ))

!		 print *, '-----------------------------------------------------------------------------------------'
 ! This doesn't work according to  Fouart et al. 2015, so commented out for now.  
   !  thick_PL = min(thick_PL,vxL/WL)
   !  thick_PR = min(thick_PR,vxR/WR)

   !  thick_ML = min(thick_ML,vxL/WL)
   !  thick_MR = min(thick_ML,vxR/WR)

!     lamb4L = half*(three*zetaL - one)*thin_PL + three*half*(one - zetaL)*thick_PL
!     lamb4R = half*(three*zetaR - one)*thin_PR + three*half*(one - zetaR)*thick_PR

!     lamb5L = half*(three*zetaL - one)*thin_ML + three*half*(one - zetaL)*thick_ML
!     lamb5R = half*(three*zetaR - one)*thin_MR + three*half*(one - zetaR)*thick_MR

	 
  end if
     lamb_max = max(0.0d0,lamb1R,lamb1L,lamb2R,lamb2L,lamb3R,lamb3L,lamb4R,lamb4L,lamb5R,lamb5L)

     lamb_min = min(0.0d0,lamb1R,lamb1L,lamb2R,lamb2L,lamb3R,lamb3L,lamb4R,lamb4L,lamb5R,lamb5L)

     ! ::::::::::::::::::::::::::::::::

     ! Fluxes along the x-direction

     f1L = vxL*DL
     f1R = vxR*DR

     f2L = vxL*SxL + pressL
     f2R = vxR*SxR + pressR

     f3L = SxL
     f3R = SxR

     f4L = u5L 
     f4R = u5R

     f5L = PrL 
     f5R = PrR

     ! Fluxes at the each intercell
     ! ============================

     FluxHyd(1,i) = (lamb_max*f1L - lamb_min*f1R + lamb_max*lamb_min*(DR - DL))/(lamb_max-lamb_min)
     FluxHyd(2,i) = (lamb_max*f2L - lamb_min*f2R + lamb_max*lamb_min*(SxR - SxL))/(lamb_max-lamb_min)
     FluxHyd(3,i) = (lamb_max*f3L - lamb_min*f3R + lamb_max*lamb_min*(TauR - TauL))/(lamb_max-lamb_min)
do j=1,ngroups
     FluxRad(1,i,j) = (lamb_max*f4L - lamb_min*f4R + lamb_max*lamb_min*(u4R - u4L))/(lamb_max-lamb_min)
     FluxRad(2,i,j) = (lamb_max*f5L - lamb_min*f5R + lamb_max*lamb_min*(u5R - u5L))/(lamb_max-lamb_min)
 if(close_type.eq.'M1')then
	 !correccion al flujo correspondiente a la densidad de energia radiada, segun  Fouart et al. 2015

vel_advL = four*WL**2/(two*WL**2 + one)*vxL
vel_advR = four*WR**2/(two*WR**2 + one)*vxR

if((vel_advL.ge.zero).and.(vel_advR.ge.zero)) then
	
	FluxRad(1,i,j) = tanh(idx/sqrt(chi(i,j)*chi(i+1,j) ))*FluxRad(1,i,j) + &
	            (one - tanh(idx/sqrt(chi(i,j)*chi(i+1,j) )))*(four*third*WL**2*vxL*ErComovilL - &
			    third*sqrt(W(i)*W(i+1))/sqrt(chi(i,j)*chi(i+1,j) )*(one + &
				(sqrt(vx(i)*vx(i+1)))**2)*(ErComovil(i+1,j) - ErComovil(i,j))*idx)
			  
	else if((vel_advL.le.zero).and.(vel_advR.le.zero))then
		
	FluxRad(1,i,j) = tanh(idx/sqrt(chi(i,j)*chi(i+1,j) ))*FluxRad(1,i,j) + &
		        (one - tanh(idx/sqrt(chi(i,j)*chi(i+1,j) )))*(four*third*WR**2*vxR*ErComovilR - &
				third*sqrt(W(i)*W(i+1))/sqrt(chi(i,j)*chi(i+1,j) )*(one + &
				(sqrt(vx(i)*vx(i+1)))**2)*(ErComovil(i+1,j) - ErComovil(i,j))*idx)
				  	
    else if((vel_advL.ge.zero.and.vel_advR.le.zero).or.(vel_advL.le.zero.and.vel_advR.ge.zero))then
		
	FluxRad(1,i,j) = tanh(idx/sqrt(chi(i,j)*chi(i+1,j) ))*FluxRad(1,i,j) + &
		        (one - tanh(idx/sqrt(chi(i,j)*chi(i+1,j) )))*( &
			  - third*sqrt(W(i)*W(i+1))/sqrt(chi(i,j)*chi(i+1,j) )*(one + &
			  (sqrt(vx(i)*vx(i+1)))**2)*(ErComovil(i+1,j) - ErComovil(i,j))*idx)
end if
		
				  
	  FluxRad(2,i,j) = tanh( idx/sqrt(chi(i,j)*chi(i+1,j)) )*FluxRad(2,i,j) + &
	                   (one - tanh( idx/sqrt(chi(i,j)*chi(i+1,j)) ))*(f5R + f5L)*half
  end if
end do  
  end do
  
end subroutine flux_hlle


 
