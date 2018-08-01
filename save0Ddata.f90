
  subroutine save0Ddata(fval,base_name,first_index)

  use arrays
  use global_numbers

  implicit none

  character(len=20) filestatus
  logical firstcall
  data firstcall / .true. /
  save firstcall

  character(len=*), intent(IN) :: base_name
  real(kind=8), dimension(-g_pts:Nx+g_pts), intent(IN) :: fval

  character(len=256) :: filename

  integer i,first_index
  real(kind=8) min, max, nm1, nm2, dspace

  if (res_num.eq.1) then
    filename = base_name // '_1.tl'
  else if (res_num.eq.2) then
    filename = base_name // '_2.tl'
  else if (res_num.eq.3) then
    filename = base_name // '_3.tl'
  else if (res_num.eq.4) then
    filename = base_name // '_4.tl'
  else if (res_num.eq.5) then
    filename = base_name // '_5.tl'
  end if

  if (first_index.eq.0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if


! --->   Calculating scalars

  max = fval(0)
  min = fval(0)
  nm1 = 0.0D0
  nm2 = 0.0D0

!  do i=1,Nx
!     if (fval(i)>max) max = fval(i)
!     if (fval(i)<min) min = fval(i)
!  end do

!  do i=1,Nx
!     nm1 = nm1 + 0.5D0*(dabs(fval(i-1)) + dabs(fval(i)))*dspace
!     nm2 = nm2 + 0.5D0*(fval(i-1)**2 + fval(i)**2)*dspace
!     if (fval(i)>max) max = fval(i)
!     if (fval(i)<min) min = fval(i)
!  end do

  nm2 = dsqrt(nm2)


  if (filestatus=='replace') then
     open(1,file=filename,form='formatted',status=filestatus)
  write(1,*) '#Time                fval'
  else
     open(1,file=filename,form='formatted',status=filestatus,position='append')
  end if
     write(1,*) t,fval(i)
  close(1)

  end subroutine save0Ddata
