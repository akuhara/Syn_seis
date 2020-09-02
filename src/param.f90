module param
  integer, parameter :: clen_max = 100
  integer, parameter :: nlay_max = 100
  integer, parameter :: npts_max = 16384
  integer, parameter :: nhalf    = npts_max / 2 + 1
  complex(kind(0d0)), parameter :: ei = (0.d0, 1.d0)
  real(8), parameter :: pi2 = 2 * 3.1415926535897931
end module param
