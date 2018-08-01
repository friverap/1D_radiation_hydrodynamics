
subroutine evolve

  use arrays
  use global_numbers

  implicit none

  integer i,j,k,l
  integer num_steps

  real(kind=8) dt_temp

 ! call allocate
  call grid

  print *,'================================='
  print *, "Your job name: ", trim(adjustl(jobname))
  print *, "Using ", Nx, " radial zones and ", g_pts, " ghosts zones"
  print *, 'xmin=',xmin
  print *, 'xmax=',xmax
  print *, 'dx=',dx
  print *, 'dt=',dt
  print *, 'courant=',dt/dx
  print *, "Using ",trim(adjustl(Limiter))," reconstruction"
  print *, "Using ",trim(adjustl(Riemann)), " for Rad-hydro"
  print *, "Using ",trim(adjustl(close_type)), " closure for radiation"
  print *, "Using ",ngroups, " groups for radiation"
  print *, 'Final time=', dt*Nt
  print *,'=================================='

  t = zero
  uHyd_p   = zero  
  uRad_p   = zero  

  print *,'----------------------------'
  print *,'|  Time step  |    Time    |'
  print *,'----------------------------'

  write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',0,'    | ',t,'  |'

  call initial
  call diagnostics

!   call saveLum(Lum_aux,'Lum',0)
!  call save0Ddata(tao,'tao',0)
!  call save0Ddata(maxT,'maxT',0)
!  call save0Ddata(maxTr,'maxTr',0)
!  call save0Ddata(ar,'ar',0)

  call save1Ddata(rho,'rho',0)
  call save1Ddata(press,'press',0)
  call save1Ddata(gradPress,'gradPress',0)
  call save1Ddata(vx,'vx',0)
  call save1DRaddata(temp,'T',0)
  
  call save1DRaddata(Er,'Er',0)
  call save1DRaddata(Fr,'Fr',0)
  call save1DRaddata(Pr,'Pr',0)
  call save1DRaddata(Prthick,'Prthick',0)
  call save1DRaddata(ErComovil,'ErComovil',0)
  call save1DRaddata(FrComovil,'FrComovil',0)
  call save1DRaddata(FrComovil/ErComovil,'FComovil',0)
  call save1DRaddata(PrComovil,'PrComovil',0)
  call save1DRaddata(PrComovil/ErComovil,'fac_eddComovil',0)
  call save1DRaddata(optical_depth,'optical_depth',0)
  call save1DRaddata(temp_r,'Tr',0)
  call save1DRaddata(force_r,'radiativeForce',0)

!  call save1Ddata(PrComovil/press,'ratioPress',0)
  
  do l=1,Nt

     t = t + dt

     if (mod(l,every_1D).eq.0) then
        write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',l,'    | ',t,'  |'
     end if

     uHyd_p  = uHyd
	 uRad_p  = uRad

     do k=1,3

        call rhs

        if (k.eq.1) then
           dt_temp = dt
           uHyd = uHyd_p + dt_temp*half*rhs_uHyd
		   uRad = uRad_p + dt_temp*half*rhs_uRad
        else if (k.eq.2) then
           dt_temp = 0.5d0*dt
           uRad = uRad + dt_temp*half*rhs_uRad
        else if (k.eq.3) then
           dt_temp = 0.5d0*dt
           uHyd = uHyd + dt_temp*rhs_uHyd 
           uRad = uRad + dt_temp*half*rhs_uRad 
        end if

        call BC            

     end do

     call diagnostics

     if (mod(l,every_0D).eq.0) then

!        call saveLum(Lum_aux,'Lum',1)
!        call save0Ddata(tao,'tao',1)
!        call save0Ddata(maxT,'maxT',1)
!        call save0Ddata(maxTr,'maxTr',1)
!        call save0Ddata(ar,'ar',1)

     end if

     if (mod(l,every_1D).eq.0) then
        call save1Ddata(rho,'rho',1)
        call save1Ddata(press,'press',1)
        call save1Ddata(gradPress,'gradPress',1) 
        call save1Ddata(vx,'vx',1)
        call save1DRaddata(temp,'T',1)
		
        call save1DRaddata(Er,'Er',1)
        call save1DRaddata(Fr,'Fr',1)
        call save1DRaddata(Pr,'Pr',1)
        call save1DRaddata(Prthick,'Prthick',1)
        call save1DRaddata(ErComovil,'ErComovil',1)
        call save1DRaddata(FrComovil,'FrComovil',1)
        call save1DRaddata(FrComovil/ErComovil,'FComovil',1)
        call save1DRaddata(PrComovil,'PrComovil',1)
        call save1DRaddata(PrComovil/ErComovil,'fac_eddComovil',1) 
        call save1DRaddata(optical_depth,'optical_depth',1) 
        call save1DRaddata(temp_r,'Tr',1) 
        call save1DRaddata(force_r,'radiativeForce',1)

!        call save1Ddata(PrComovil/press,'ratioPress',1) 
     end if

  end do

end subroutine evolve

