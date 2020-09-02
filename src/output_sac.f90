subroutine output_sac(index, npts, dt, t_shift, ur, uz)
  use param
  implicit none 
  character(clen_max), intent(in) :: index
  integer, intent(in) :: npts
  real(8), intent(in) :: dt, t_shift
  real(8), intent(in) :: ur(npts_max), uz(npts_max)
  integer :: icmp, ierr, i, npre
  character(clen_max) :: r_file, z_file
  character(clen_max), dimension(2) :: files
  integer, parameter :: nhdr = 158
  real :: data(2, npts_max)

  r_file = trim(index) // ".r"
  z_file = trim(index) // ".z"
  files = (/r_file, z_file/)
  
  npre = nint((t_shift ) / dt) + 1
  
  if (npre >= 0) then
     do i = 1, npre
        data(1, i) = real(ur(npts_max - i + 1))
        data(2, i) = real(uz(npts_max - i + 1))
     end do
     do i = 1, npts_max - npre
        data(1, i + npre) = real(ur(i))
        data(2, i + npre) = real(uz(i))
     end do
  else
     write(0,*)"WARNING: time shift is short"
     do i = 1, npts_max + npre
        data(1, i) = real(ur(i - npre))
        data(2, i) = real(uz(i - npre))
     end do
     do i = 1, npre
        data(1, i) = real(ur(-npre-i))
        data(2, i) = real(uz(-npre-i))
     end do
  end if

  do icmp = 1, 2
     
     open(10, file = files(icmp), status = "unknown", iostat = ierr, &
          & access = "direct", recl = 4)
     if (ierr /= 0) then
        write(0,*)"ERROR: cannot create ", trim(files(icmp))
        stop
     end if
     do i = 1, nhdr
        write(10, rec = i) -12345
     end do
     write(10, rec = 1) real(dt)
     write(10, rec = 6) real(-t_shift)
     write(10, rec = 7) real(-t_shift + npts * dt)
     write(10, rec = 77) 6
     write(10, rec = 80) npts
     write(10, rec = 86) 1
     write(10, rec = 106) 1
     do i = 1, npts
        write(10, rec = nhdr + i) data(icmp, i)
     end do
     close(10)
  end do
  return 
end subroutine output_sac
