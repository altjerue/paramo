module constants
   implicit none

   double precision, parameter :: cspeed = 2.99792458d10, &
        pi = 3.141592653589793238462643d0, &
        sqrtpi = 1.772453850905516027298167d0, &
        mp = 1.67262158d-24, &
        me = 9.10938188d-28, &
        echarge = 4.803204d-10, &
        sigmaT = 6.6524586d-25, &
        Planckh = 6.629069d-27, &
        nuconst = 2.7992491077281560779657886d6, & ! eCharge / 2 * pi * m_e * cLight
        jmbconst = 4.8352765842234213535181620d-29, & ! 2 * pi * eCharge**2 / cLight
        ambconst = 2.6540088545187588392244135d-2, & ! pi * eCharge**2 / (me * cLight)
        lzero = dlog(1d-200), &
        lzero100 = dlog(1d-100), &
        chunche_c100g100 = 2.2619939050180366385d-6, &
        chunche_c100g20  = 2.1157699720918349273d-1

end module constants
