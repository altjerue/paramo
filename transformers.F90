module transformers
use data_types
use constants
implicit none

contains

!> Convert distance in parsec to centimeter
function pc2cm(dpc) result(dcm)
    real(dp), intent(in) :: dpc
    real(dp) :: dcm
    dcm = 3.08567758149137d18 * dpc
end function pc2cm

!> Convert distance in centimeter to parsec
function cm2pc(dcm) result(dpc)
    real(dp), intent(in) :: dcm
    real(dp) :: dpc
    dpc = dcm / 3.08567758149137d18
end function cm2pc

!> Convert time in seconds to days
!! @parameter tsec real time in seconds
function sec2dy(tsec) result(tdy)
    real(dp), intent(in) :: tsec
    real(dp) :: tdy
    tdy = tsec / 86400d0
end function sec2dy

!> Convert time in days to seconds
!! @parameter tdy real time in days
function dy2sec(tdy) result(tsec)
    real(dp), intent(in) :: tdy
    real(dp) :: tsec
    tsec = tdy * 86400d0
end function dy2sec

!> Transfrom frequency in eV to Hz
function eV2Hz(veV) result(vHz)
    implicit none
    real(dp), intent(in) :: veV
    real(dp) :: vHz
    vHz = veV * 2.4179937422321953d14
end function eV2Hz

!> Transfrom frequency in Hz to eV
function Hz2eV(vHz) result(veV)
    implicit none
    real(dp), intent(in) :: vHz
    real(dp) :: veV
    veV = vHz / 2.4179937422321953d14
end function Hz2eV

!> Convert flux density (in egs cm^{-2} s^{-1} Hz^{-1}) to janskys
!! @parameter flux real energy flux in ergs
function erg2Jy(flux) result(res)
   implicit none
   real(dp), intent(in) :: flux
   real(dp) :: res
   res = flux * 1d23
end function erg2Jy

!> Convert flux density in jansky to egs cm^{-2} s^{-1} Hz^{-1}.
!! @parameter flux real energy flux in janskys
function Jy2erg(flux) result(res)
   implicit none
   real(dp), intent(in) :: flux
   real(dp) :: res
   res = flux * 1d-23
end function Jy2erg

end module transformers
