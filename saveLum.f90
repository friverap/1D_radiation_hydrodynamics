
  subroutine saveLum(fval,base_name,first_index)

  use arrays
  use global_numbers

  implicit none

  character(len=20) filestatus
  logical firstcall
  data firstcall / .true. /
  save firstcall

  character(len=*), intent(IN) :: base_name
  real(kind=8), dimension(1:num_det), intent(IN) :: fval

  character(len=256) :: filename

  integer i, first_index

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


  if (filestatus=='replace') then
     open(1,file=filename,form='formatted',status=filestatus)
  write(1,*) '#Time                fval'
  else
     open(1,file=filename,form='formatted',status=filestatus,position='append')
  end if
     write(1,*) t,fval
  close(1)

  end subroutine saveLum

