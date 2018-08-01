
subroutine initial

  use arrays
  use global_numbers

  implicit none

  !se tomara el sistema de ecuaciones que proponen en el paper de Takahashi et al 2012, esta formulacion sera util al momento de implementar el metodo multi-grupos
  !en este trabajo, las ecuaciones de la radiacion se escribiran en el sistema de referencia mixto. mientra la ecuacion de densidad de momento y de energia se escribiran en el sistema de laboratorio, los coeficientes de absorcion/emision y dispersion se escriben en el sistema de referencia comovil.
  ! el tensor de energia-momento para la radiacion se escribe como:
  !       Er  Fxr
  !Tmunu= Fxr Pxxr
  ! y las variables conservativas, flujos y terminos fuente son (una dimension):
  !     D	  rho*W                    D*vx              0     0
  !     Sx	  rho*h*W**2*vx            Sx*vx + P        Gx    -ar*T**4*chi_t*W*vx + (chi_t + chi_s)*(W*Fxr - W*vx*Pxxr)- chi_s*W*vx*(W**2*Er-two*W**2*vx*Fxr + W**2vx**2*Pxxr)
  ! U = tau	= rho*h*W**2 - P      ;F = Sx          ;S = G0  = -chi_t*(ar*T**4*W - W*Er + W*vx*Fxr) - chi_s*(W**3*vx**2*Er + W**3*vx**2*Pxxr - (W**2 + W**2*vx**2)*W*vx*Fxr)
  !     Er	  Er                       Fxr             -G0    -G0
  !     Fxr	  Fxr                      Pxxr            -Gx    -Gx
  
  ! sin embargo las ecuacion que cierra el sistema completo de ecuaciones de la radiacion obedece una relacion de la forma P'xxr = P'xxr(E'r,F'xr), donde las variables primadas son
  ! medidas en el sistema de referancia comovil. Para obtener la presion de radiacion en el sistema de referencia del observador se le hace una transformacion de Lorentz al tensor de energia-momento de la radiacion, que combinado con la aproximacion de Eddington (1D ==> Pxxr = Er/3) se tiene:
  
  !       third*(W**2*Er - two*W*vx*Fxr) - W**2*vx**2Er + two*W*vx*Fxr + two*W**3*vx**3*Fxr/(one + W)
  ! Pxxr =--------------------------------------------------------------------------------------------
  !                one + (-third + W**2*vx**2/(one + W)**2)*W**2*vx + two*W**2*vx**2/(1 + W)
  
  
  !las transformaciones de Lorentz entre los sistemas de referencias del observador y el comovil para la densidad de energia radiada, el flujo radiado y la presion por radiacion son:
  
  ! ErComovil = W**2*(Er - two*vx*Fr + Vx**2*Pr)
  ! FrComovil = -W**2*vx*Er + W**2*(one + vx**2)*Fr - W**2*vx*Pr
  ! PrComovil = W**2*vx**2*Er - two*W**2*vx*Fr + W**2*Pr
  
  !reemplazando v --> -v se tiene la transformacion inversa
  integer i,j,k,l
  real(kind=8) funPr, funprimePr
  real(kind=8) tol,tolPr
  real(kind=8) fcomovil_iniL, fcomovil_iniR,zeta1Comovil_iniL,zeta1Comovil_iniR
  real(kind=8) Re_iniL,Re_iniR,Re1_iniL,Re1_iniR,PrThick_iniL,PrThick_iniR
  
  !los datos iniciales para las pruebas de los tubos de choques estan dados en el sistema de referencia comovil
  cte_rad = Er_iniL*(rho_iniL/p_iniL)**4
  Er_iniR = cte_rad*(p_iniR/rho_iniR)**4
  Fr_iniL = zero  !1.e-2*Er_iniL
  Fr_iniR = zero  !1.e-2*Er_iniR

  W_iniL = one/sqrt(one - vx_iniL**2)
  W_iniR = one/sqrt(one - vx_iniR**2)
  
  !cambiar los datos iniciales del sistema de referencia comovil al sistema de referencia del observador  
	  !debido a que el sistema comienza en una region opticamente gruesa, una buena aproximacion para recuperar los datos iniciales en el sistema de referencia del observador es considerear la aproximacion de Eddington
	  
 	 Pr_iniL = third*Er_iniL
	 Pr_iniR = third*Er_iniR
	 
     Er_iniObsL = W_iniL**2*(Er_iniL + two*vx_iniL*Fr_iniL + vx_iniL**2*Pr_iniL)
     Er_iniObsR = W_iniR**2*(Er_iniR + two*vx_iniR*Fr_iniR + vx_iniR**2*Pr_iniR)
     Fr_iniObsL = W_iniL**2*vx_iniL*Er_iniL + W_iniL**2*(one + vx_iniL**2)*Fr_iniL + W_iniL**2*vx_iniL*Pr_iniL
     Fr_iniObsR = W_iniR**2*vx_iniR*Er_iniR + W_iniR**2*(one + vx_iniR**2)*Fr_iniR + W_iniR**2*vx_iniR*Pr_iniR



  

  
  
  !datos iniciales en el sistema de referencia del observador
  do i = -g_pts,Nx+g_pts
     if (x(i).le.x0) then
        rho(i)   = rho_iniL
        press(i) = p_iniL
        vx(i)    = vx_iniL + 1.e-20
        Er(i,:)    = Er_iniObsL
        Fr(i,:)    = Fr_iniObsL 

        ErComovil(i,:)    = Er_iniL
        FrComovil(i,:)    = Fr_iniL 
     else
        rho(i)   = rho_iniR 
        press(i) = p_iniR
        vx(i)    = vx_iniR+ 1.e-20
        Er(i,:)    = Er_iniObsR
        Fr(i,:)    = Fr_iniObsR 
		
        ErComovil(i,:)    = Er_iniR
        FrComovil(i,:)    = Fr_iniR 
     end if
  end do

  W  = one/sqrt(one - vx**2)
 
 !actually is -gradP ;-) 
  do i = -g_pts+1,Nx+g_pts-1
  gradPress(i) = -(press(i+1) - press(i-1))/dx*half
  end do
  
  gradPress(0)  = -(-0.5d0*(3.0D0*gradPress(0)  - 4.0D0*gradPress(1)    + gradPress(2)   ) / dx)
  gradPress(Nx) = -( 0.5d0*(3.0D0*gradPress(Nx) - 4.0D0*gradPress(Nx-1) + gradPress(Nx-2)) / dx)
  
if(eos_type.eq.'id') then
  ee = press/(rho*(gamma - one))
  hh = one + ee + press/rho
  CS = sqrt((gamma*press)/(hh*rho))
else if(eos_type.eq.'tm') then
	if (Er_iniL/p_iniL.gt.three)then
        stop "Why are you doing the tm approach within the radiation-pressure-dominated scenario?"
     endif
  hh = five*half*(press/rho) + sqrt(9.0d0/four*(press/rho)**2 + one)
  CS = sqrt( third*(press/rho)/hh*(five*hh - 8.0d0*press/rho)/(hh - press/rho) )
!  CS = sqrt( (five* (press/rho)*sqrt((press/rho)**2 + four/9.0d0) +three*(press/rho)**2)/&
!            (12.0d0*(press/rho)*sqrt((press/rho)**2 + four/9.0d0) +12.0d0*(press/rho)**2 + two) )

end if

  
  uHyd(1,:) = max(Floor,rho*W)
  uHyd(2,:) = rho*hh*W**2*vx
  uHyd(3,:) = max(Floor,rho*hh*W**2 - press)
  uRad(1,:,:) = max(Er,floorrad)
  uRad(2,:,:) = Fr

  
  ErComovil = max(ErComovil,floorrad)
  FrComovil = FrComovil
 

  !aproximacion de Eddington
   if(close_type.eq.'dif') then
	   
	   !definir la presion de radiacion en el sistema de referencia de un observador euleriano
	   do j=1,ngroups 
	   Re(:,j)  = third*(W(:)**2*uRad(1,:,j) - two*W(:)**2*vx(:)*uRad(2,:,j)) - W(:)**2*vx(:)**2*uRad(1,:,j) + &
	              two*W(:)*vx(:)*uRad(2,:,j) + two*W(:)**3*vx(:)**3*uRad(2,:,j)/(one + W(:))
	   end do
	   Re1 = one + (-third + W**2*vx**2/(one + W)**2)*W**2*vx**2 + two*W**2*vx**2/(one + W)

	   do j=1,ngroups
	   PrThick(:,j) = Re(:,j)/Re1(:)
	   Pr(:,j) = max(PrThick(:,j),floor)
	   end do
	   
	   !aproximacion M1
  else if(close_type.eq.'M1') then

  
    do i = -g_pts,Nx+g_pts
		do j= 1,ngroups
    f(i,j) = uRad(2,i,j)/uRad(1,i,j)	  
  	  if (f(i,j).gt.one)then
  		   print *, "Fx/cE > 1 : i,group,fx=", i,j,f(i,j)
  		  limx(i,j) = one/f(i,j)
            uRad(2,i,j) = uRad(2,i,j) * limx(i,j)
           f(i,j) = uRad(2,i,j)/uRad(1,i,j)
        end if
	    end do
    end do
   

 	! takahashi & Ohsuga 2013  
 	 zeta1 = (three + four*f**2)/(five + two*sqrt(four - three*f**2))	

	pr = zeta1*uRad(1,:,:)	

	  
   end if
   
  
   !recuperar variables de la radiacion en el sistema de referencia comovil
  do i = -g_pts,Nx+g_pts
   do j=1,ngroups
   PrComovil(i,j) = W(i)**2*vx(i)**2*uRad(1,i,j) - two*W(i)**2*vx(i)*uRad(2,i,j) + W(i)**2*Pr(i,j)
   end do
  end do
  
  !opacidades cuerpo gris, opacidad constante 
  do i = -g_pts,Nx+g_pts
   do j=1,ngroups
  chi_t(i,j) = rho(i)*opacity_a
  chi_s(i,j) = rho(i)*opacity_s
   end do
  end do  
  chi   = chi_t + chi_s
  
  !profundidad optica medida a lo largo de la direccion x (Abramowicz, Novikov & Paczynski (1991), Eq. (3.1)) 
  !tao = int( chi*W*( 1-vx ) )*dx
    tao = zero
    do i=1-g_pts,Nx+g_pts
		do j=1,ngroups
       tao = tao + half*( dabs(chi(i-1,j)*W(i-1)*(one - vx(i-1)) + chi(i,j)*W(i)*(one - vx(i))) )*dx
        end do
    end do
  
    optical_depth(0,:) = tao
    do i=1-g_pts,Nx+g_pts
  	 optical_depth(i,:) = optical_depth(i-1,:) + half*( dabs(chi(i-1,:)*W(i-1)*(one - vx(i-1)) + chi(i,:)*W(i)*(one - vx(i))) )*dx
    end do

    do i = -g_pts,Nx+g_pts
     do j=1,ngroups	
 Temp(i,j)  = press(i)/rho(j)
     end do
   end do 
   Temp_r  = sqrt(sqrt(ErComovil/cte_rad))
   
  !fuerza de radiacion
  do i = -g_pts,Nx+g_pts
   do j=1,ngroups	
  force_r(i,j) = -chi_t(i,j)*(cte_rad*Temp(i,j)**4*W(i) - W(i)*uRad(1,i,j) + W(i)*vx(i)*uRad(2,i,j)) - &
			      chi_s(i,j)*(W(j)**3*vx(i)**2*uRad(1,i,j) + W(i)**3*vx(i)**2*Pr(i,j) - &
				  (W(i)**2 + W(i)**2*vx(i)**2)*W(i)*vx(i)*uRad(2,i,j)) 
   end do
  end do 
  
  !funcion de Planck
  do i = -g_pts,Nx+g_pts
   do j=1,ngroups	
  PlanckFunction(i,j) = freq(j)/(exp(freq(j)/Temp(i,j)) - one)
 end do
end do 
  
  !Third moment
  
  
  
!detectors
 do j=1,num_det
     Det(j) =  initial_det + dble(j-1)*space_between_det
  end do

  print *, '----------------------------------------'
  do j = 1,num_det
     print *, 'detector', j, 'at =',det(j)
  end do

  print *, 'WL = ', one/sqrt(one - vx_iniL**2)
  print *, 'ErR=', Er_iniR
  print *, 'Pr_iniL/Pg_iniL=', Er_iniL/(three*p_iniL)
  print *, 'optical depth(0) =', tao
  print *, 'optical depth(Nx,1) =', optical_depth(Nx,1)
  print *, 'ErObsL=', Er_iniObsL
  print *, 'ErObsR=', Er_iniObsR
  print *, 'cte_rad=', cte_rad
  print *, "Done with initial data :-)"

  print *, "Begin time integration loop:"
  print *, '----------------------------------------'
end subroutine initial
