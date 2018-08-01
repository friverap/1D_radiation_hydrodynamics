
subroutine BC

  use arrays
  use global_numbers

  implicit none


  integer i

  ! Outflow boundary conditions

  uHyd(:,-g_pts)   = uHyd(:,-g_pts+1)
  uHyd(:,Nx+g_pts) = uHyd(:,Nx+g_pts-1)

  uRad(:,-g_pts,:)   = uRad(:,-g_pts+1,:)
  uRad(:,Nx+g_pts,:) = uRad(:,Nx+g_pts-1,:)

end subroutine BC

