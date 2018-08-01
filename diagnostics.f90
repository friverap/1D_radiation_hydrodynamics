
  subroutine diagnostics

  use arrays
  use global_numbers

  implicit none

  integer i,j
  integer before2_det
  integer after2_det

  real(kind=8) idx

  idx = one/dx
! NaN checker

  do i=1,Nx-1

    if (vx(i).ne.vx(i)) then
      print *, 'There are NaNs on v at x=',x(i)
      print *, 'Time=', t
      stop
    end if

    if (rho(i).ne.rho(i)) then
      print *, 'There are NaNs on rho at x=',x(i)
      print *, 'Time=', t
      stop
    end if

    if (press(i).ne.press(i)) then
      print *, 'There are NaNs on press at x=',x(i)
      print *, 'Time=', t
      stop
    end if
	
do j=1,ngroups
    if (Fr(i,j).ne.Fr(i,j)) then
      print *, 'There are NaNs on Fr at x=',x(i)
      print *, 'Time=', t
      stop
    end if

    if (Er(i,j).ne.Er(i,j)) then
      print *, 'There are NaNs on Er at x=',x(i)
      print *, 'Time=', t
      stop
    end if
end do
  end do

  end subroutine diagnostics
