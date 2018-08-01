
subroutine grid

  use arrays
  use global_numbers

  implicit none

  integer i,j
  
  !espacio

  dx  = (xmax - xmin)/dble(Nx)

  do i=-g_pts,Nx+g_pts
     x(i)  = xmin + dble(i) * dx
  end do

  dt = courant * dx

  !frecuencias
  
  deltaf = (fmax - fmin)/dble(ngroups)

  do j =0,ngroups
	  freq(j) = fmin + dble(j)*deltaf
  end do

end subroutine grid
