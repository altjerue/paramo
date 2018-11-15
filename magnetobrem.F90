module magnetobrem
   use constants
   use misc, only:chebev
   implicit none
   
   public
   
   integer :: numXX, numL, numSS
   double precision :: globgmax, globgmin, LLmax!, chunche
   double precision, allocatable, dimension(:) :: XX, SS, dXX, dSS, jtable, atable, jRMAtable, aRMAtable, LLmin
   
contains
   
   !============================================================================
   !  SUBROUTINE: load_mb_table
   !
   !============================================================================
   subroutine load_mb_table(file)
      use hdfman
      use hdf5
      implicit none
      character :: file*(*)
      
      integer(HID_T) :: file_id,group_id
      
      write(*,*) "Reading table: ",file
      
      call hm_init
      call hm_open(file, file_id)
      call hm_gopen(file_id, "Params", group_id)
      
      call hm_read0_int(group_id,    numXX,    "num_chi")
      call hm_read0_int(group_id,    numL,     "num_gam")
      call hm_read0_int(group_id,    numSS,    "num_q")
      call hm_read0_double(group_id, globgmin, "g_min")
      call hm_read0_double(group_id, globgmax, "g_max")
      
      call hm_gclose(group_id)
      
      allocate(XX(numXX), SS(numSS), dXX(numXX - 1), dSS(numSS - 1), LLmin(numXX))
      allocate(jtable(numXX*numL*numSS), atable(numXX*numL*numSS), jRMAtable(numXX*numL*numSS), aRMAtable(numXX*numL*numSS))
      
      ! call hm_read0_double(file_id, chunche,  "chunche")
      
      call hm_read1_double(file_id, numXX, XX,    "chi")
      call hm_read1_double(file_id, numSS, SS,    "q")
      call hm_read1_double(file_id, numXX, LLmin, "xi_min")
      call hm_read1_double(file_id, numXX * numL * numSS, jtable,    "disTable")
      call hm_read1_double(file_id, numXX * numL * numSS, atable,    "adisTable")
      call hm_read1_double(file_id, numXX * numL * numSS, jRMAtable, "RMATable")
      call hm_read1_double(file_id, numXX * numL * numSS, aRMAtable, "aRMATable")
      
      dXX(1:numXX - 1) = 1d0 / (XX(2 : numXX) - XX(1:numXX - 1))
      dSS(1:numSS - 1) = 1d0 / (SS(2 : numSS) - SS(1:numSS - 1))
      
      call hm_close(file_id)
      call hm_finalize
      
      LLmax = 0d0
      
   end subroutine load_mb_table

   
   !============================================================================
   !  chtab_int
   !  ---------
   !
   !  Interpolation routine from Chebychev coefficients
   !
   !  chi  : scalar <- normalized frequency
   !  lor  : scalar <- normalized Lorentz factor \xi
   !  q    : scalar <- power-law index
   !  chtab: array  <- Chebychev coefficients
   !============================================================================
   function chtab_int(chi,lor,q,chtab) result(res)
      use misc, only: chebev
      implicit none
      
      double precision, dimension(:), intent(in):: chtab
      double precision, intent(in) :: chi,lor,q
      integer :: i,j
      double precision :: cx,cl,cq,u,v,u1,v1,res
      double precision :: valij,valipj,valijp,valipjp
      double precision, dimension(numL) :: coefs
      
      cx = dlog(chi)
      cl = dlog(lor)
      cq = q
      
      i = max( 1, min( idnint(dble(numXX - 1) * (cx - XX(1)) / (XX(numXX) - XX(1))) + 1, numXX - 1 ) )
      if (cx < XX(i    )) i = max(1       , i - 1)
      if (cx > XX(i + 1)) i = min(numXX - 1, i + 1)
      
      j = max( 1, min( idnint(dble(numSS - 1) * (cq - SS(1)) / (SS(numSS) - SS(1))) + 1, numSS - 1 ) )
      if (cq < SS(j    )) j = max(1        , j - 1)
      if (cq > SS(j + 1)) j = min(numSS - 1, j + 1)
      
      u = (cx - XX(i)) * dXX(i)
      v = (cq - SS(j)) * dSS(j)
      
      u1 = 1d0 - u
      v1 = 1d0 - v
      
#if 1
      if (cl < LLmin(i) .or. cl < LLmin(i + 1)) then
         res = lzero
      else
         coefs = chtab(1 + j * numL + i * numL * numSS:numL + j * numL + i * numL * numSS)
         valij = chebev(cl,coefs,numL,LLmin(i + 1),LLmax)
         
         coefs = chtab(1 + (j - 1) * numL + i * numL * numSS:numL + (j - 1) * numL + i * numL * numSS)
         valijp = chebev(cl,coefs,numL,LLmin(i + 1),LLmax)
         
         coefs = chtab(1 + (j - 1) * numL + (i - 1) * numL * numSS:numL + (j - 1) * numL + (i - 1) * numL * numSS)
         valipjp = chebev(cl,coefs,numL,LLmin(i),LLmax)
         
         coefs = chtab(1 + j * numL + (i - 1) * numL * numSS:numL + j * numL + (i - 1) * numL * numSS)
         valipj = chebev(cl,coefs,numL,LLmin(i),LLmax)
         
         res = u1 * v1 * valipjp + &
         u * v1 * valijp + &
         u1 * v * valipj + &
         u * v * valij
      end if
#endif
      
#if 0
      if (cl < LLmin(i)) then
         if (cl < LLmin(i + 1)) then
            res = lzero
         else
            coefs = chtab(1 + j * numL + i * numL * numSS:numL + j * numL + i * numL * numSS)
            valij = chebev(cl,coefs,numL,LLmin(i + 1),LLmax)
            
            coefs = chtab(1 + (j - 1) * numL + i * numL * numSS:numL + (j - 1) * numL + i * numL * numSS)
            valijp = chebev(cl,coefs,numL,LLmin(i + 1),LLmax)
            
            res = u * v * valij + &
            u * v1 * valijp
         end if
      else
         if (cl < LLmin(i + 1)) then
            coefs = chtab(1 + (j - 1) * numL + (i - 1) * numL * numSS:numL + (j - 1) * numL + (i - 1) * numL * numSS)
            valipjp = chebev(cl,coefs,numL,LLmin(i),LLmax)
            
            coefs = chtab(1 + j * numL + (i - 1) * numL * numSS:numL + j * numL + (i - 1) * numL * numSS)
            valipj = chebev(cl,coefs,numL,LLmin(i),LLmax)
            
            res = u1 * v1 * valipjp + &
            u1 * v * valipj
         else
            coefs = chtab(1 + j * numL + i * numL * numSS:numL + j * numL + i * numL * numSS)
            valij = chebev(cl,coefs,numL,LLmin(i + 1),LLmax)
            
            coefs = chtab(1 + (j - 1) * numL + i * numL * numSS:numL + (j - 1) * numL + i * numL * numSS)
            valijp = chebev(cl,coefs,numL,LLmin(i + 1),LLmax)
            
            coefs = chtab(1 + (j - 1) * numL + (i - 1) * numL * numSS:numL + (j - 1) * numL + (i - 1) * numL * numSS)
            valipjp = chebev(cl,coefs,numL,LLmin(i),LLmax)
            
            coefs = chtab(1 + j * numL + (i - 1) * numL * numSS:numL + j * numL + (i - 1) * numL * numSS)
            valipj = chebev(cl,coefs,numL,LLmin(i),LLmax)
            
            res = u1 * v1 * valipjp + &
            u * v1 * valijp + &
            u1 * v * valipj + &
            u * v * valij
         end if
      end if
#endif
      
      res = dmax1(lzero, res)
      
      return
   end function chtab_int
   
   
   !-----------------------------{  Emissivity  }-------------------------------
   
   function j_mb(nu,B,n0,gmin,gmax,qq,RMAfunc,c0,jmbtab) result(emiss)
      ! ========================================================================
      !   j_mb:
      !   -----
      !
      !   nu:      scalar   <- frequency
      !   B:       scalar   <- Magnetic field
      !   n0:      scalar   <- number densiti of gmin
      !   gmin:    scalar   <- chunk minimum Lorentz factor
      !   gmax:    scalar   <- chunk maximum Lorentz factor
      !   qq:      scalar   <- power-law index
      !   RMAfunc: function <- function to be used for analytic integration
      !   jmbtab:  array    <- table to be used for interpolation (emissivity)
      ! ========================================================================
      use misc, only: an_error
      use anaFormulae!, only: RMA_qromb
      use pwlinteg
      implicit none
      interface
         function RMAfunc(c,g) result(res)
            double precision :: res
            double precision, intent(in) :: c, g
         end function RMAfunc
      end interface
      double precision, intent(in) :: nu, B, gmin, gmax, qq, n0, c0
      double precision, intent(in), dimension(:) :: jmbtab
      integer :: i
      double precision :: emiss, chi, nu_b, I2, jmbconst2, loggmin, loggmax
      double precision :: I3max, I3min, I3diff, I3rel, I3RMAmax, I3RMAmin, I3RMAdiff, I3RMArel, I3prod
      
      jmbconst2 = 0.125d0 * jmbconst ! <- missing 1/8 factor
      nu_b = nuconst * B
      chi = nu / nu_b
      
      if (chi < dexp(XX(1))) then
         write(*,*) 'chi =', chi
         call an_error('j_nu: nu below table limits in emissivity')
      end if
      
      i = max( 1, min( idnint(dble(numXX - 1) * (dlog(chi) - XX(1)) / (XX(numXX) - XX(1))) + 1, numXX - 1 ) )
      loggmin = dlog(gmin / globgmax)
      loggmax = dlog(gmax / globgmax)
! #if 0
!       if (LLmin(i) > loggmax) then
!          emiss = 0d0
!          return
!       end if
! #else
      if (LLmin(i) >= loggmax) then
         emiss = 0d0
         return
      else if (LLmin(i) >= loggmin .and. LLmin(i) < loggmax) then
         I3min = dexp(chtab_int(chi, dexp(LLmin(i)), qq, jmbtab))
         I3max = dexp(chtab_int(chi, gmax / globgmax, qq, jmbtab))
      else
         I3min = dexp(chtab_int(chi, gmin / globgmax, qq, jmbtab))
         I3max = dexp(chtab_int(chi, gmax / globgmax, qq, jmbtab))
      end if
! #endif
      
      I3diff = I3min - I3max
      I3rel = dabs(I3diff) / dabs(I3min)
      
      !!$if (I3min <= I3max .or. I3min <= 1d-100) then
      !!$   emiss = 0d0
      !!$   return
      !!$end if
      
      !!$if (dabs(I3diff) * globgmax**(1d0 - qq) > 1d-1) then
      !!$   I2 = RMA_qromb(chi, qq, dlog(gmin), dlog(gmax), 1d0, RMA)
      !!$else
      if (I3rel < 2d-5) I3diff = 0d0 !c0*RMA_qromb(chi, qq, dlog(gmin/globgmax), dlog(gmax/globgmax), globgmax, RMAfunc)
      if (0.5d0 * LN2(globgmax, 1d-9) * (qq - 1d0)**2 < 1d-3) then
         I2 = (1d0 - LN1(globgmax, 1d-9) * (qq - 1d0)) * I3diff
      else
         I2 = globgmax**(1d0 - qq) * I3diff
         if (I2 > 1d3) I2 = 0d0!RMA_qromb(chi, qq, dlog(gmin), dlog(gmax), 1d0, RMAfunc)
         !!$   write (*,"(9ES15.6)") nu,qq,I3min,I3max,I3diff,I3rel,globgmax**(1d0 - qq),globgmax**(1d0 - qq)*I3diff,I2
      end if
      !   if (qq > 0d0 .and. I2 > 5d-2) I2 = 0d0!RMA_qromb(chi, qq, dlog(gmin), dlog(gmax), 1d0, RMAfunc)
      !!$end if
      
#if 0
      if (0.5d0 * LN2(globgmax, 1d-9) * (qq - 8d0)**2 < 1d-1) then
         if (I3diff > 1d-1) I3diff = 1.25d0 * I3max
         I2 = (1d0 - LN1(globgmax, 1d-9) * (qq - 8d0)) * I3diff / globgmax**7
      else if (0.5d0 * LN2(globgmax, 1d-9) * (qq - 1d0)**2 < 1d-1) then
         if (I3diff > 1d-1) I3diff = 1.25d0 * I3max
         I2 = (1d0 - LN1(globgmax, 1d-9) * (qq - 1d0)) * I3diff
      else
         if (I3diff > 1d-1) I3diff = 1.25d0 * I3max
         I2 = globgmax**(1d0 - qq) * I3diff
      else if (0.5d0 * LN2(globgmax, 1d-9) * (qq + 8d0)**2 < 1d-1) then
         if (I3diff > 1d-1 / globgmax**9) I3diff = 1.25d0 * I3max / globgmax**9
         I2 = globgmax**9 * (1d0 - LN1(globgmax, 1d-9) * (qq + 8d0)) * I3diff
      else
         if (I3diff > 1d-1) I3diff = 1.25d0 * I3max
         I2 = globgmax**(1d0 - qq) * I3diff
      end if
#endif
      !!$write (*,"(8ES15.6)") nu,gmin,gmax,qq,I3min,I3max,I3min-I3max,I2
      
      emiss = dmax1(1d-200, jmbconst2 * nu_b * n0 * I2 * gmin**qq)
      
      return
   end function j_mb
   
   
   function j_mb_qromb(nu,B,n0,gmin,gmax,qq,c0,RMAfunc) result(emiss)
      use anaFormulae, only: RMA_qromb
      implicit none
      interface
         function RMAfunc(c,g) result(res)
            double precision :: res
            double precision, intent(in) :: c,g
         end function RMAfunc
      end interface
      
      double precision, intent(in) :: nu,B,gmin,gmax,qq,n0,c0
      integer :: i
      double precision :: emiss,chi,nu_b,I2,jmbconst2
      
      jmbconst2 = 0.25d0 * jmbconst ! <- missing 1/4 factor
      nu_b = nuconst * B
      chi = nu / nu_b
      
      !!$I2 = globgmax**(1d0 - q) * &
      !!$     ( RMA_qromb(chi, q, dlog(gmin / globgmax), 0d0, globgmax) - &
      !!$     RMA_qromb(chi, q, dlog(gmax / globgmax), 0d0, globgmax) )
      !!$I2 = RMA_qromb(chi, q, dlog(gmin / globgmax), 0d0, globgmax) - RMA_qromb(chi, q, dlog(gmax / globgmax), 0d0, globgmax)
      !!$emiss = dmax1(1d-200, jmbconst * nu_b * n0 * gmin**q * I2)
      !!$emiss = dmax1(1d-200, jmbconst * nu_b * n0 * gmin**q * globgmax**(1d0 - q) * I2)
      
      i = max( 1, min( idnint(dble(numXX - 1) * (dlog(chi) - XX(1)) / (XX(numXX) - XX(1))) + 1, numXX - 1 ) )
      if (LLmin(i) > dlog(gmax / globgmax)) then
         emiss = 0d0
      else if (LLmin(i) > dlog(gmin / globgmax) .and. LLmin(i) <= dlog(gmax / globgmax)) then
         I2 = c0 * RMA_qromb(chi, qq, LLmin(i), dlog(gmax/globgmax), globgmax, RMAfunc)
         emiss = dmax1(1d-200, jmbconst2 * nu_b * n0 * I2 * gmin**qq)
      else
         I2 = c0 * RMA_qromb(chi, qq, dlog(gmin / globgmax), dlog(gmax / globgmax), globgmax, RMAfunc)
         emiss = dmax1(1d-200, jmbconst2 * nu_b * n0 * I2 * gmin**qq)
      end if
      
      return
   end function j_mb_qromb
   
   
   ! ----------------------------{  Absorption  }------------------------------
   
   function a_mb(nu,B,n0,gmin,gmax,qq,aRMAfunc,c0,ambtab) result(absor)
      ! =======================================================================
      !  a_mb:
      !  -----
      !
      !  nu:      scalar   <- frequency
      !  B:       scalar   <- Magnetic field
      !  n0:      scalar   <- number densiti of gmin
      !  gmin:    scalar   <- chunk minimum Lorentz factor
      !  gmax:    scalar   <- chunk maximum Lorentz factor
      !  qq:      scalar   <- power-law index
      !  RMAfunc: function <- function to be used for analytic integration
      !  ambtab:  array    <- table to be used for interpolation (absorption)
      ! =======================================================================
      use misc, only: an_error
      use anaFormulae!, only: ARMA_qromb
      use pwlinteg
      implicit none
      interface
         function aRMAfunc(c,g) result(res)
            double precision :: res
            double precision, intent(in) :: c,g
         end function aRMAfunc
      end interface
      integer :: i
      double precision, intent(in) :: nu,B,gmin,gmax,qq,n0,c0
      double precision, intent(in), dimension(:) :: ambtab
      double precision :: absor,chi,nu_b,A2,ambconst2,A3max,A3min,A3diff,A3rel
      double precision :: loggmin, loggmax
      
      ambconst2 = 2.5d-1 * 1.5625d-2 * ambconst ! <- missing 1/64/4 factor
      nu_b = nuconst * B
      chi = nu / nu_b
      
      
      if (chi < dexp(XX(1))) then
         !!$   A2 = dmax1(1d-200, ARMA_qromb(chi, qq, dlog(gmin), dlog(gmax), 1d0,  RMAfunc))
         !!$   absor = dmax1(1d-200, ambconst2 * nu_b * n0 * A2 * gmin**qq)
         write(*,*) 'chi =', chi
         call an_error('a_nu: nu below table limits in emissivity')
      end if

      i = max( 1, min( idnint(dble(numXX - 1) * (dlog(chi) - XX(1)) / (XX(numXX) - XX(1))) + 1, numXX - 1 ) )
      loggmin = dlog(gmin / globgmax)
      loggmax = dlog(gmax / globgmax)

      if (LLmin(i) >= loggmax) then
         absor = 0d0
         return
      else if (LLmin(i) >= loggmin .and. LLmin(i) < loggmax) then
         A3min = dexp(chtab_int(chi, dexp(LLmin(i)), qq, ambtab))
         A3max = dexp(chtab_int(chi, gmax / globgmax, qq, ambtab))
      else
         A3min = dexp(chtab_int(chi, gmin / globgmax, qq, ambtab))
         A3max = dexp(chtab_int(chi, gmax / globgmax, qq, ambtab))
      end if

      ! A3min = dexp(chtab_int(chi, gmin / globgmax, qq, ambtab))
      ! A3max = dexp(chtab_int(chi, gmax / globgmax, qq, ambtab))
      ! if (A3min <= A3max .or. A3min < 1d-100) then
      !    absor = c0*ARMA_qromb(chi, qq, dlog(gmin/globgmax), dlog(gmax/globgmax), globgmax, RMAfunc)!0d0!1d-200
      !    !print*,A3min,A3max
      !    return
      ! end if

      A3diff = A3min - A3max
      A3rel = dabs(A3diff) / dabs(A3min)
      
      if (A3rel < 2d-5) A3diff = 0d0 !c0*RMA_qromb(chi, qq, dlog(gmin/globgmax), dlog(gmax/globgmax), globgmax, RMAfunc)
      if (0.5d0 * LN2(globgmax, 1d-9) * (qq - 1d0)**2 < 1d-3) then
         A2 = (1d0 - LN1(globgmax, 1d-9) * (qq - 1d0)) * A3diff
      else
         A2 = globgmax**(1d0 - qq) * A3diff
         if (A2 > 1d3) A2 = 0d0!RMA_qromb(chi, qq, dlog(gmin), dlog(gmax), 1d0, RMAfunc)
         !!$   write (*,"(9ES15.6)") nu,qq,I3min,I3max,I3diff,I3rel,globgmax**(1d0 - qq),globgmax**(1d0 - qq)*I3diff,I2
      end if
      
#if 0      
      if (1d0 - A3max / A3min < 1d-6) then
         A2 = c0*ARMA_qromb(chi, qq, dlog(gmin/globgmax), dlog(gmax/globgmax), globgmax, RMAfunc)
      else
         if (0.5d0 * dlog(globgmax)**2 * (qq - 1d0)**2 < 1d-6) then
            A2 = (1d0 - dlog(globgmax) * (qq - 1d0)) * (A3min - A3max)
         else
            A2 = globgmax**(1d0 - qq) * (A3min - A3max)
         end if
         if (A2 > 1d5) A2 = c0*ARMA_qromb(chi, qq, dlog(gmin/globgmax), dlog(gmax/globgmax), globgmax, RMAfunc)
      end if
#endif      
      absor = dmax1(1d-200, ambconst2 * nu_b * n0 * A2 * gmin**qq / nu**2)
      
      return
   end function a_mb
   
   
   function a_mb_qromb(nu,B,n0,gmin,gmax,qq,c0,RMAfunc) result(absor)
      use anaFormulae, only: ARMA_qromb
      implicit none
      interface
         function RMAfunc(c,g) result(res)
            double precision :: res
            double precision, intent(in) :: c,g
         end function RMAfunc
      end interface
      double precision :: absor,chi,nu_b,A2,ambconst2
      double precision, intent(in) :: nu,B,gmin,gmax,qq,n0,c0
      
      ambconst2 = 0.25d0 * ambconst ! <- missing 1/4 factor
      nu_b = nuconst * B
      chi = nu / nu_b
      
      !!$I2 = globgmax**(1d0 - q) * &
      !!$     ( RMA_qromb(chi, q, dlog(gmin / globgmax), 0d0, globgmax) - &
      !!$     RMA_qromb(chi, q, dlog(gmax / globgmax), 0d0, globgmax) )
      !!$I2 = RMA_qromb(chi, q, dlog(gmin / globgmax), 0d0, globgmax) - RMA_qromb(chi, q, dlog(gmax / globgmax), 0d0, globgmax)
      !!$emiss = dmax1(1d-200, jmbconst * nu_b * n0 * gmin**q * I2)
      !!$emiss = dmax1(1d-200, jmbconst * nu_b * n0 * gmin**q * globgmax**(1d0 - q) * I2)
      
      A2 = c0*ARMA_qromb(chi, qq, dlog(gmin/globgmax), dlog(gmax/globgmax), globgmax, RMAfunc)
      
      absor = dmax1(1d-200, ambconst2 * nu_b * n0 * A2 * gmin**qq / nu**2)
      
      return
   end function a_mb_qromb
   
end module magnetobrem
