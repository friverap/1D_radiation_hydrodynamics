
subroutine recv_primitives

  use arrays
  use global_numbers

  implicit none

   integer i,l,j
  integer before_det, after_det
  real(kind=8) x_i,x_i1
  real(kind=8) fun, funprime,funprad, dfunprad
  real(kind=8) funPr, dfunPr,Pr1,Pr2, xl, xh
  real(kind=8) funPr1, funPr2
  real(kind=8) tol,tolPr,pradguess
  
  real(kind=8):: W_,Re_,Re1_,PrThick_
  real(kind=8):: fComovil_,vx_,Er_,Fr_,Pr_
  real(kind=8):: zeta1Comovil_,zetaComovil_ 
  real(kind=8):: dfComovil_,dzeta1Comovil_,dzetaComovil_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(eos_type.eq.'id') then

  SS = sqrt(uHyd(2,:)**2)
!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(tol) SCHEDULE(STATIC,chunk_size)
  do i = -g_pts,Nx+g_pts

     P_GUESS(i) = 0.5d0 * (p_iniL + p_iniR)
     P_STARK(i) = P_GUESS(i)

     do l = 1,1000

        P_STARK_p(i) = P_STARK(i)

        fun = - P_STARK_p(i) + (gamma - one)*uHyd(3,i) &
             - (SS(i)**2*(gamma - one))/( P_STARK_p(i) + uHyd(3,i)) &
             - uHyd(1,i)*(gamma - one)*( sqrt(one - SS(i)**2/(( &
              P_STARK_p(i) + uHyd(3,i))**2)))

        funprime = - one + (SS(i)**2*(gamma - one))/(( P_STARK_p(i) + uHyd(3,i))**2) &
             - (uHyd(1,i)*SS(i)**2*(gamma - one))/(( P_STARK_p(i) &
             + uHyd(3,i))**3*sqrt(one - SS(i)**2/(( P_STARK_p(i) + uHyd(3,i))**2)))


        P_STARK(i) = P_STARK_p(i) - fun/funprime

        tol = two*abs(P_STARK(i) - P_STARK_p(i))/(P_STARK(i) + P_STARK_p(i))

        if (tol.lt.1.e-6) then
           press(i) = P_STARK(i)
           exit
        else
        end if
        if (l.gt.1000) then
           print *, 'The Newton_Raphson routine did not converge'
           print *, 'l=',l
           print *, 'i=',i
           print *, 'x=',x(i)
           stop
        end if

     end do
  end do
!$OMP END PARALLEL DO
  press = max(press,Floor)

  vx = uHyd(2,:)/( press + uHyd(3,:)) 
  W  = one/sqrt(one - vx**2)

else if(eos_type.eq.'tm') then!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

  SS = sqrt(uHyd(2,:)**2)
!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j,tol) SCHEDULE(STATIC,chunk_size)
  do i = -g_pts,Nx+g_pts

     P_GUESS(i) = 0.5d0 * (p_iniL + p_iniR)
     P_STARK(i) = P_GUESS(i)

     do l = 1,1000

        P_STARK_p(i) = P_STARK(i)

gam(i)  = one/sqrt( one - SS(i)**2/(uHyd(3,i)+P_STARK_p(i))**2 )
rhou(i) = gam(i)/uHyd(1,i)

        fun = uHyd(1,i)*( five/two*P_STARK_p(i)*rhou(i) + sqrt( 9.0d0/four*P_STARK_p(i)**2*rhou(i)**2 + one) )*gam(i)&
              - P_STARK_p(i) - uHyd(3,i)

        funprime = - one + uHyd(1,i)*(uHyd(3,i)+P_STARK_p(i))**2/(sqrt((uHyd(3,i)+P_STARK_p(i))**2 -SS(i)**2))**3*&
    ( sqrt(9.0d0*P_STARK_p(i)**2*(uHyd(3,i)+P_STARK_p(i))**2/(four*uHyd(1,i)**2*((uHyd(3,i)+P_STARK_p(i))**2 -SS(i)**2)) + one ) +&
            five*P_STARK_p(i)*   (uHyd(3,i)+P_STARK_p(i))   /(two*uHyd(1,i)*sqrt((uHyd(3,i)+P_STARK_p(i))**2 -SS(i)**2)) ) +&!primer termino
                 uHyd(1,i)*(uHyd(3,i)+P_STARK_p(i))/sqrt((uHyd(3,i)+P_STARK_p(i))**2 - SS(i)**2)*&
   ( (- 9.0d0*P_STARK_p(i)**2*( uHyd(3,i)+P_STARK_p(i) )**3/( two*uHyd(1,i)**2*(( uHyd(3,i)+P_STARK_p(i) )**2 - SS(i)**2)**2 )+ &
        9.0d0*P_STARK_p(i)**2*( uHyd(3,i)+P_STARK_p(i) )   /( two*uHyd(1,i)**2*(( uHyd(3,i)+P_STARK_p(i) )**2 - SS(i)**2)    )+ &
        9.0d0*P_STARK_p(i)   *( uHyd(3,i)+P_STARK_p(i) )**2/( two*uHyd(1,i)**2*(( uHyd(3,i)+P_STARK_p(i) )**2 - SS(i)**2)    ))/&
(two*sqrt(9.0d0*P_STARK_p(i)**2*(uHyd(3,i)+P_STARK_p(i))**2/(four*uHyd(1,i)**2*(( uHyd(3,i)+P_STARK_p(i) )**2 - SS(i)**2) )+ one))&
 -five*P_STARK_p(i)*( uHyd(3,i)+P_STARK_p(i) )**2/( two*uHyd(1,i)*(sqrt(( uHyd(3,i)+P_STARK_p(i) )**2 - SS(i)**2) )**3)&
 +five*             ( uHyd(3,i)+P_STARK_p(i) )   /( two*uHyd(1,i)* sqrt(( uHyd(3,i)+P_STARK_p(i) )**2 - SS(i)**2) )   &
 +five*P_STARK_p(i)                           /( two*uHyd(1,i)* sqrt(( uHyd(3,i)+P_STARK_p(i) )**2 - SS(i)**2) )  )&!segundo termino   
 +uHyd(1,i)/sqrt(( uHyd(3,i)+P_STARK_p(i) )**2 - SS(i)**2)*( &
  sqrt(9.0d0*P_STARK_p(i)**2*(uHyd(3,i)+P_STARK_p(i))**2/(four*uHyd(1,i)**2*((uHyd(3,i)+P_STARK_p(i))**2 -SS(i)**2)) + one ) +&
           five*P_STARK_p(i)*(uHyd(3,i)+P_STARK_p(i))/(two*uHyd(1,i)*sqrt((uHyd(3,i)+P_STARK_p(i))**2 -SS(i)**2)) )               

          
        P_STARK(i) = P_STARK_p(i) - fun/funprime

        tol = two*abs(P_STARK(i) - P_STARK_p(i))/(P_STARK(i) + P_STARK_p(i))

        if (tol.lt.1.e-6) then
           press(i) = P_STARK(i)
           exit
        else
        end if
        if (l.gt.1000) then
           print *, 'The Newton_Raphson routine did not converge'
           print *, 'l=',l
           print *, 'i=',i
           print *, 'x=',x(i)
           stop
        end if

     end do
  end do
!$OMP END PARALLEL DO
  press = max(press,Floor)
  
  
  vx = uHyd(2,:)/( press + uHyd(3,:)) 
  W  = one/sqrt(one - vx**2)

end if

!$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(STATIC,chunk_size)
do i = -g_pts+1,Nx+g_pts-1
gradPress(i) = -(press(i+1) - press(i-1))/dx*half
end do
!$OMP END PARALLEL DO
gradPress(0)  = -(-0.5d0*(3.0D0*gradPress(0)  - 4.0D0*gradPress(1)    + gradPress(2)   ) / dx)
gradPress(Nx) = -( 0.5d0*(3.0D0*gradPress(Nx) - 4.0D0*gradPress(Nx-1) + gradPress(Nx-2)) / dx)

  rho = uHyd(1,:)/W
  rho = max(rho,Floor) 
if(eos_type.eq.'id') then
  ee = press/(rho*(gamma - one))
  hh = one + ee + press/rho
  CS = sqrt((gamma*press)/(hh*rho))
else if(eos_type.eq.'tm') then
  hh = five*half*(press/rho) + sqrt(9.0d0/four*(press/rho)**2+one)
!  CS = sqrt( third*(press/rho)/hh*(five*hh - 8.0d0*press/rho)/(hh - press/rho) )
  CS = sqrt( (five* (press/rho)*sqrt((press/rho)**2 + four/9.0d0) +three*(press/rho)**2)/&
            (12.0d0*(press/rho)*sqrt((press/rho)**2 + four/9.0d0) +12.0d0*(press/rho)**2 + two) )

end if

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)
do i = -g_pts,Nx+g_pts
	do j=1,ngroups
    if (uRad(1,i,j) <= zero) then
       print *, '** cal_ray **'
       print *, 'i,group,Er(i)= ',i,j,uRad(1,i,j)
       uRad(1,i,j) = floorrad
    end if
    end do
end do
!$OMP END PARALLEL DO 
     Er  = uRad(1,:,:)
     Fr  = uRad(2,:,:)
	 



 if(close_type.eq.'dif') then
	 !definir la presion de radiacion en el sistema del observador para el medio opticamente grueso
	 do j=1,ngroups 
	 Re(:,j)  = third*(W(:)**2*uRad(1,:,j) - two*W(:)**2*vx(:)*uRad(2,:,j)) - W(:)**2*vx(:)**2*uRad(1,:,j) + &
	            two*W(:)*vx(:)*uRad(2,:,j) + two*W(:)**3*vx(:)**3*uRad(2,:,j)/(one + W(:))
	 end do
	 Re1 = one + (-third + W**2*vx**2/(one + W)**2)*W**2*vx**2 + two*W**2*vx**2/(one + W)

	 do j=1,ngroups
	     PrThick(:,j) = Re(:,j)/Re1(:)
	     Pr(:,j) = max(PrThick(:,j),floor)
	 end do	
	 
 else if(close_type.eq.'M1') then  
 !$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)
     do i = -g_pts,Nx+g_pts
		 do j=1,ngroups
	   f(i,j) = max(uRad(2,i,j)/uRad(1,i,j),floorrad)
   	  if (f(i,j).gt.one)then
   		   print *, "Fx/cE > 1 : i,group,fx=", i,j,f(i,j)
   		  limx(i,j) = one/f(i,j)
             uRad(2,i,j) = uRad(2,i,j) * limx(i,j)
            f(i,j) = max(uRad(2,i,j)/uRad(1,i,j),floorrad)
         end if
	     end do
     end do
!$OMP END PARALLEL DO 
 	! takahashi & Ohsuga 2013  
 	 zeta1 = (three + four*f**2)/(five + two*sqrt(four - three*f**2))
     pr = zeta1*uRad(1,:,:)	


 end if
 !$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)
 do i = -g_pts,Nx+g_pts
	 do j=1,ngroups
     if (pr(i,j) <= zero) then
         print *, 'Warning'
         print *, 'pr =', pr(i,j)
         print *, 'i=',i
         print *, 'group=',j
        pr(i,j) =floorrad
     end if
     end do
 end do
 !$OMP END PARALLEL DO 
 
 !recuperar variables de la radiacion en el sistema de referencia comovil
  !$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)
 do i = -g_pts,Nx+g_pts
  do j=1,ngroups
ErComovil(i,j) = W(i)**2*(uRad(1,i,j) - two*vx(i)*uRad(2,i,j) + vx(i)**2*Pr(i,j))
  end do
 end do

 !$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)
do i = -g_pts,Nx+g_pts
	do j=1,ngroups
    if (ErComovil(i,j) <= zero) then
       print *, '** cal_ray **'
       print *, 'i,group,ErComovil= ',i,j,ErComovil(i,j)
       ErComovil(i,j) = floorrad
    end if
    end do
end do
 !$OMP END PARALLEL DO 
 !$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)  
do i = -g_pts,Nx+g_pts
	do j=1,ngroups
FrComovil(i,j) = (-W(i)**2*vx(i)*uRad(1,i,j) + W(i)**2*(one + vx(i)**2)*uRad(2,i,j) - W(i)**2*vx(i)*Pr(i,j))
PrComovil(i,j) = W(i)**2*vx(i)**2*uRad(1,i,j) - two*W(i)**2*vx(i)*uRad(2,i,j) + W(i)**2*Pr(i,j)
    end do
end do

 !$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)
do i = -g_pts,Nx+g_pts
	do j=1,ngroups
    if (prComovil(i,j) <= zero) then
        print *, 'Warning'
        print *, 'prComovil =', prComovil(i,j)
        print *, 'i=',i
        print *, 'group=',j
       prComovil(i,j) = floorrad
    end if
   end do
end do
 !$OMP END PARALLEL DO 
 !calcular las opacidades
 
 !opacidades cuerpo gris, opacidad constante 
  !$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)
 do i = -g_pts,Nx+g_pts
  do j=1,ngroups
 chi_t(i,j) = rho(i)*opacity_a
 chi_s(i,j) = rho(i)*opacity_s
  end do
 end do 
  !$OMP END PARALLEL DO 
 chi   = chi_t + chi_s
 

Temp_r  = sqrt(sqrt(ErComovil/cte_rad))
  
  !funcion de Planck 
  do i = -g_pts,Nx+g_pts
   do j=1,ngroups	
  PlanckFunction(i,j) = freq(j)/(exp(freq(j)/Temp(i,j)) - one)
 end do
end do 


Er_p(0,:) = zero!3.0D0*Er(1) - 3.0D0*Er(2) + Er(3) 
lum(0,:) = 3.0D0*lum(1,:) - 3.0D0*lum(2,:) + lum(3,:) 
  do i=1-g_pts,Nx+g_pts
	do j=1,ngroups
     Er_p(i,j) = Er_p(i-1,j) !- Er(l-1)
     Lum(i,j)  = (Er(i,j)-Er_p(i,j))/dt
    end do
  end do

  before_det = 0
  after_det  = 0

  do l=num_det,1,-1
	  do j=1,ngroups
    if (det(l).lt.x(Nx)) then

      do i=-g_pts,Nx+g_pts
        if (det(l).le.x(i)) then
          before_det = i-1
          after_det  = i
          exit
        end if
      end do
     x_i = x(after_det)
     x_i1 = x(before_det)
     Lum_aux(l,j) = Lum(after_det,j) + (det(l) - x_i)*((Lum(after_det,j) - Lum(before_det,j))/(x(after_det) - x_i1))
    else
     Lum_aux(l,j) = 3.0D0*Lum(l+1,j) - 3.0D0*Lum(l+2,j) + Lum(l+3,j)
    end if
	end do
   end do
   
  ! print *, '----'
end subroutine recv_primitives
