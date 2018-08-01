
subroutine sources

  use global_numbers
  use arrays

  implicit none
  integer i,l,j
  real(kind=8) fun, funprime
  real(kind=8) tol



  !profundidad optica medida a lo largo de la direccion x (Abramowicz, Novikov & Paczynski (1991), Eq. (3.1)) 
  !tao = int( chi*W*( 1-vx ) )*dx

tao = zero
!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)
do i=1-g_pts,Nx+g_pts
	do j=1,ngroups
   tao = tao + half*( dabs(chi(i-1,j)*W(i-1)*(one - vx(i-1)) + chi(i,j)*W(i)*(one - vx(i))) )*dx
    end do
end do
!$OMP END PARALLEL DO

optical_depth(0,:) = tao
!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)
do i=1-g_pts,Nx+g_pts
 optical_depth(i,:) = optical_depth(i-1,:) + half*( dabs(chi(i-1,:)*W(i-1)*(one - vx(i-1)) + chi(i,:)*W(i)*(one - vx(i))) )*dx
end do
!$OMP END PARALLEL DO
  
!calcular la temperatura como en el paper de cuesta-martinez 2014 ecuacion A6

!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j,tol) SCHEDULE(STATIC,chunk_size)
do i = -g_pts,Nx+g_pts
	do j=1,ngroups
		
     T_GUESS(i,j) = p_iniL/rho_iniL
     T_STARK(i,j) = T_GUESS(i,j)

     do l = 1,1000

        T_STARK_p(i,j) = T_STARK(i,j)

         fun = press(i) + PrComovil(i,j)*( one-exp(-chi(i,j)*W(i)*(one - vx(i))*dx) ) - rho(i)*T_STARK_p(i,j) - &
		       third*cte_rad*T_STARK_p(i,j)**4*( one-exp(-chi(i,j)*W(i)*(one - vx(i))*dx) )

         funprime = - rho(i) - four*third*cte_rad*T_STARK_p(i,j)**3*( one-exp(-chi(i,j)*W(i)*(one - vx(i))*dx) )


!          fun =  press(i) + PrComovil(i,j)*( one-exp(-chi(i,j)*W(i)*(one - vx(i))*dx) ) - rho(i)*T_STARK_p(i,j) - &
!                 ( one-exp(-chi(i,j)*W(i)*(one - vx(i))*dx) )*( three*cte_rad*T_STARK_p(i,j)**4 + &
!                 four*FrComovil(i,j)**2/(cte_rad*T_STARK_p(i,j)**4) )/&
!                  ( five + two*sqrt(four - three*(FrComovil(i,j)/(cte_rad*T_STARK_p(i,j)**4))**2) )
				  
!         funprime = -rho(i)-(one-exp(-chi(i,j)*W(i)*(one - vx(i))*dx))*( (12.0d0*cte_rad*T_STARK_p(i,j)**3 - &
!		       16.0d0*FrComovil(i,j)**2/cte_rad/T_STARK_p(i,j)**5)/&
!               (five + two*sqrt(four - three*(FrComovil(i,j)/(cte_rad*T_STARK_p(i,j)**4))**2)) - &
!               24.0d0*FrComovil(i,j)**2*(four*FrComovil(i,j)**2/cte_rad/T_STARK_p(i,j)**4 + three*cte_rad*T_STARK_p(i,j)**4)/&
!               (cte_rad**2*T_STARK_p(i,j)**9*(five + two*sqrt(four - three*(FrComovil(i,j)/(cte_rad*T_STARK_p(i,j)**4))**2))**2*&
!                sqrt(four - three*(FrComovil(i,j)/(cte_rad*T_STARK_p(i,j)**4))**2)) )
           
        T_STARK(i,j) = T_STARK_p(i,j) - fun/funprime

        tol = two*abs(T_STARK(i,j) - T_STARK_p(i,j))/(T_STARK(i,j) + T_STARK_p(i,j))

        if (tol.lt.1.e-6) then
           Temp(i,j) = T_STARK(i,j)
           exit
        else
        end if
        if (l.gt.1000) then
           print *, 'The Newton_Raphson routine for T did not converge'
           print *, 'l=',l
           print *, 'i=',i
           print *, 'x=',x(i)
           stop
        end if

     end do
    end do
  end do
!$OMP END PARALLEL DO

  !fuerza de radiacion
  !$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)
  do i = -g_pts,Nx+g_pts
   do j=1,ngroups	
  force_r(i,j) = -chi_t(i,j)*(cte_rad*Temp(i,j)**4*W(i) - W(i)*uRad(1,i,j) + W(i)*vx(i)*uRad(2,i,j)) - &
			      chi_s(i,j)*(W(j)**3*vx(i)**2*uRad(1,i,j) + W(i)**3*vx(i)**2*Pr(i,j) - &
				  (W(i)**2 + W(i)**2*vx(i)**2)*W(i)*vx(i)*uRad(2,i,j)) 
   end do
  end do 
  !$OMP END PARALLEL DO
  
!  maxT  = Temp(0)
!  maxTr = Temp_r(0)

!  do i=1,Nx
!     if (Temp(i)>maxT) maxT = Temp(i)
!     if (Temp_r(i)>maxTr) maxTr = Temp_r(i)
!  end do

!terminos fuente
!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)
do i = -g_pts,Nx+g_pts
 do j=1,ngroups	
  sHyd(1,i,j) = zero

  sHyd(2,i,j) =  -cte_rad*Temp(i,j)**4*chi_t(i,j)*W(i)*vx(i) + (chi_t(i,j) + chi_s(i,j))*(W(i)*uRad(2,i,j) - &
                 W(i)*vx(i)*Pr(i,j)) - chi_s(i,j)*W(i)*vx(i)*(W(i)**2*uRad(1,i,j)-two*W(i)**2*vx(i)*uRad(2,i,j) + &
				 W(i)**2*vx(i)**2*Pr(i,j))

  sHyd(3,i,j) = -chi_t(i,j)*(cte_rad*Temp(i,j)**4*W(i) - W(i)*uRad(1,i,j) + W(i)*vx(i)*uRad(2,i,j)) - &
                 chi_s(i,j)*(W(i)**3*vx(i)**2*uRad(1,i,j) + W(i)**3*vx(i)**2*Pr(i,j) - (W(i)**2 + &
				 W(i)**2*vx(i)**2)*W(i)*vx(i)*uRad(2,i,j)) 
  
  sRad(1,i,j) = -sHyd(3,i,j)

  sRad(2,i,j) = -sHyd(2,i,j)
 end do
end do   
!$OMP END PARALLEL DO


end subroutine sources
