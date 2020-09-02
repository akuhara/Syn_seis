!-----------------------------------------------------------------------
subroutine read_input(index, npts, dt, t_shift, nlay ,rayp, ipha, &
     & alpha, beta, rho, h)
  use param
  implicit none 
  character(clen_max)  :: index
  integer, intent(out) :: ipha, nlay, npts
  real(8), intent(out) :: dt, rayp, t_shift
  real(8), intent(out) :: alpha(nlay_max), beta(nlay_max)
  real(8), intent(out) :: rho(nlay_max), h(nlay_max)
  integer :: i, ierr
  
  read(5, *, iostat = ierr) index
  if (ierr /= 0) then
     write(0,*)"Line 1: output file index"
     stop
  end if
  read(5, *, iostat = ierr) npts, dt, t_shift
  if (ierr /= 0) then
     write(0,*)"Line 2: [# time element] [sampling interval (s)] " // &
          & "[Time shift (s)]"
     stop
  end if
  read(5, *, iostat = ierr) nlay, rayp, ipha
  if (ierr /= 0) then
     write(0,*)"Line 3: [# layers] [ray parameter (s/km)] " // &
          & "[phase type: 1 for P, -1 for S]"
     stop
  end if
  do i = 1, nlay
     read(5, *, iostat = ierr) alpha(i), beta(i), rho(i), h(i)
     if (ierr /= 0) then
        write(0,*)"Line 4~: [Vp (km/s)] [Vs (km/s)] " // &
             & "[Density (g/cm^3)] [Thickness (km)]"
        stop
     end if
  end do
  return 
end subroutine read_input
