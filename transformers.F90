module transformers
use data_types
use constants
implicit none

contains

function pc2cm(dpc) result(dcm)
    !
    !  Description:
    !    Convert distance in parsec to centimeter
    !
    real(dp), intent(in) :: dpc
    real(dp) :: dcm
    dcm = 3.08567758149137d18 * dpc
end function pc2cm

function cm2pc(dcm) result(dpc)
    !
    !  Description:
    !    Convert distance in centimeter to parsec
    !
    real(dp), intent(in) :: dcm
    real(dp) :: dpc
    dpc = dcm / 3.08567758149137d18
end function cm2pc

function sec2dy(tsec) result(tdy)
    !
    !  Description:
    !    Convert time in seconds to days
    !
    real(dp), intent(in) :: tsec
    real(dp) :: tdy
    tdy = tsec / 86400d0
end function sec2dy

function dy2sec(tdy) result(tsec)
    !
    !  Description:
    !    Convert time in seconds to days
    !
    real(dp), intent(in) :: tdy
    real(dp) :: tsec
    tsec = tdy * 86400d0
end function dy2sec

function eV2Hz(veV) result(vHz)
    !  Description:
    !    Transfrom frequency in eV to Hz
    !
    implicit none
    real(dp) :: veV
    real(dp) :: vHz
    vHz = veV * 2.4179937422321953d14
end function eV2Hz

function Hz2eV(vHz) result(veV)
    !  Description:
    !    Transfrom frequency in Hz to eV
    !
    implicit none
    real(dp) :: vHz
    real(dp) :: veV
    veV = vHz / 2.4179937422321953d14
end function Hz2eV

end module transformers