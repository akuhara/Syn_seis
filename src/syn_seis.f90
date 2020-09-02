!-----------------------------------------------------------------------
! Progaram to calculate seismograms by propagator matrix method
!   written by Takeshi Akuhara (akuhara@eri.u-tokyo.ac.jp)
!-----------------------------------------------------------------------
program syn_seis
  use param
  implicit none 
  integer :: ipha, nlay, npts
  real(8) :: alpha(nlay_max), beta(nlay_max), rho(nlay_max), h(nlay_max)
  real(8) :: dt, rayp, t_shift
  character(clen_max) :: index
  real(8) :: ur(npts_max), uz(npts_max)
  

  call read_input(index, npts, dt, t_shift, nlay, rayp, ipha, &
       & alpha, beta, rho, h)
  call fwd_seis(nlay, rayp, dt, ipha, &
       & alpha, beta, rho, h, ur, uz)
  call output_sac(index, npts, dt, t_shift, ur, uz)
  
  
end program syn_seis

