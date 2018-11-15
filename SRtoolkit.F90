! ************************************************************************
!  Special Relativity toolkit
!  ==========================
!
!  - bofg  : function
!  - gofb  : function
!  - gbofg : function
!  
! ************************************************************************
module SRtoolkit
   use constants
   use K2
   implicit none
   
   double precision, parameter :: sqrt2=1.414213562373095048801688d0
   double precision, parameter :: isqrt2=0.707106781186547524400844d0
   
   public

   interface bofg
      module procedure bofg_s
      module procedure bofg_v
   end interface bofg

   interface gofb
      module procedure gofb_s
      module procedure gofb_v
   end interface gofb

   interface gbofg
      module procedure gbofg_s
      module procedure gbofg_v
   end interface gbofg

   interface gofgb
      module procedure gofgb_s
      module procedure gofgb_v
   end interface gofgb

   interface RMaxwell
      module procedure RMaxwell_s
      module procedure RMaxwell_v
   end interface RMaxwell

   interface powlaw_dis
      module procedure powlaw_dis_s
      module procedure powlaw_dis_v
   end interface powlaw_dis

   interface one_over_gb2
      module procedure one_over_gb2_s
      module procedure one_over_gb2_v
   end interface one_over_gb2

   interface nu_obs_f
      module procedure nu_obs_s
      module procedure nu_obs_v
   end interface nu_obs_f

   interface nu_com_f
      module procedure nu_com_s
      module procedure nu_com_v
   end interface nu_com_f


   interface t_obs_f
      module procedure t_obs_s
      module procedure t_obs_v
   end interface t_obs_f

   interface t_com_f
      module procedure t_com_s
      module procedure t_com_v
   end interface t_com_f

contains
   
   !
   !      ::::: Converting Lorentz factor to beta :::::
   !
   function bofg_s(g) result(b)
      double precision, intent(in) :: g
      double precision :: b
      if ( g <= 1d0 ) then
         b = 0d0
      else
         b = dsqrt(1d0 - 1d0 / g**2)
      end if
   end function bofg_s
   
   function bofg_v(g) result(b)
      double precision, intent(in), dimension(:) :: g
      integer :: i
      double precision, dimension(size(g)) :: b
      do i = 1, size(g)
         if ( g(i) <= 1d0 ) then
            b(i) = 0d0
         else
            b(i) = dsqrt(1d0 - 1d0 / g(i)**2)
         end if
      end do
   end function bofg_v
   
   !
   !      ::::: Converting speed, beta, to Lorentz factor :::::
   !
   function gofb_s(b) result(g)
      double precision, intent(in) :: b
      double precision :: g
      g = 1d0 / dsqrt(1d0 - b**2)
   end function gofb_s
   
   function gofb_v(b) result(g)
      double precision, intent(in), dimension(:) :: b
      double precision, dimension(size(b)) :: g
      g = 1d0 / dsqrt(1d0 - b**2)
   end function gofb_v
   
   !
   !      ::::: Calculating gamma * beta(gamma) :::::
   !
   function gbofg_s(g) result(p)
      double precision, intent(in) :: g
      double precision :: p
      p = dsqrt(g**2 - 1d0)
   end function gbofg_s
   
   function gbofg_v(g) result(p)
      double precision, intent(in), dimension(:) :: g
      double precision, dimension(size(g)) :: p
      p = dsqrt(g**2 - 1d0)
   end function gbofg_v
   
   !
   !     :::::  1 / (gamma * beta(gamma))^2 = 1 / (gamma^2 - 1)  :::::
   !
   function one_over_gb2_s(g) result(p)
      double precision, intent(in) :: g
      double precision :: p
      if ( g < 1.5d0 ) then
         p = 0.5d0 * (g - 1d0) - 0.25d0 + 0.125d0 * (g - 1d0)
      else
         p = 1d0 / (g**2 - 1d0)
      end if
   end function one_over_gb2_s
   
   function one_over_gb2_v(g) result(p)
      integer :: i
      double precision, intent(in), dimension(:) :: g
      double precision, dimension(size(g)) :: p
      do i=1,size(g)
         if ( g(i) < 1.5d0 ) then
            p(i) = 0.5d0 * (g(i) - 1d0) - 0.25d0 + 0.125d0 * (g(i) - 1d0)
         else
            p(i) = 1d0 / (g(i)**2 - 1d0)
         end if
      end do
   end function one_over_gb2_v
   
   
   !
   !     ::::: Calculating gamma(gamma * beta) :::::
   !
   function gofgb_s(l) result(g)
      double precision, intent(in) :: l
      double precision :: g
      if (l > 1d-2) then
         g = dsqrt(l**2 + 1d0)
      else
         g = 1d0 + 0.5d0 * l**2 - 0.125d0 * l**4
      end if
   end function gofgb_s
   
   function gofgb_v(l) result(g)
      double precision, intent(in), dimension(:) :: l
      integer :: i
      double precision, dimension(size(l)) :: g
      do i=1,size(g)
         if (l(i) > 1d-2) then
            g(i) = dsqrt(l(i)**2 + 1d0)
         else
            g(i) = 1d0 + 0.5d0 * l(i)**2 - 0.125d0 * l(i)**4
         end if
      end do
   end function gofgb_v


   !
   !     ::::::  Doppler factor  ::::::
   !
   function Doppler(gamma, mu) result(D)
      implicit none
      double precision, intent(in) :: gamma, mu
      double precision :: D
      D = 1d0 / (gamma * (1d0 + bofg(gamma) * mu))
   end function Doppler


   !
   !     ::::::  Viewing angle in the observer frame   ::::::
   !
   function mu_obs_f(gamma, mu) result(muo)
      implicit none
      double precision, intent(in) :: gamma, mu
      double precision :: muo
      muo = (mu + bofg(gamma)) / (1 + bofg(gamma) * mu)
   end function mu_obs_f


   !
   !     ::::::  Viewing angle in the comoving frame   ::::::
   !
   function mu_com_f(gamma, mu) result(muc)
      implicit none
      double precision, intent(in) :: gamma, mu
      double precision :: muc
      muc = (mu - bofg(gamma)) / (1 - bofg(gamma) * mu)
   end function mu_com_f


   !
   !     ::::::  Frequency in the observer frame  ::::::
   !
   function nu_obs_s(nu, z, gamma, mu) result(nuobs)
      implicit none
      double precision, intent(in) :: nu, z, gamma, mu
      double precision :: D,nuobs
      D = Doppler(gamma, mu)
      nuobs = nu * D / (1d0 + z)
    end function nu_obs_s

    function nu_obs_v(nu, z, gamma, mu) result(nuobs)
      implicit none
      double precision, intent(in) :: z, gamma, mu
      double precision, intent(in), dimension(:) :: nu
      double precision :: D
      double precision, dimension(size(nu)) :: nuobs
      D = Doppler(gamma, mu)
      nuobs = nu * D / (1d0 + z)
    end function nu_obs_v


   !
   !     ::::::  Frequency in the comoving frame  ::::::
   !
   function nu_com_s(nu, z, gamma, mu) result(nucom)
      implicit none
      double precision, intent(in) :: nu, z, gamma, mu
      double precision :: D,nucom
      D = Doppler(gamma, mu)
      nucom = nu * (1d0 + z) / D
   end function nu_com_s

    function nu_com_v(nu, z, gamma, mu) result(nucom)
      implicit none
      double precision, intent(in) :: z, gamma, mu
      double precision, intent(in), dimension(:) :: nu
      double precision :: D
      double precision, dimension(size(nu)) :: nucom
      D = Doppler(gamma, mu)
      nucom = nu * (1d0 + z) / D
   end function nu_com_v


   !
   !     ::::::  Time in the observer frame  ::::::
   !
   ! NOTE: Argument 'x' is the position in the comoving frame
   !
   function t_obs_s(t, z, gamma, x, mu) result(tobs)
      implicit none
      double precision, intent(in) :: t, z, gamma, x, mu
      double precision :: tobs
      tobs = (1d0 + z) * ((t / Doppler(gamma, mu)) + (gamma * x * (bofg(gamma) - mu)) / cspeed)
   end function t_obs_s

   function t_obs_v(t, z, gamma, x, mu) result(tobs)
      implicit none
      double precision, intent(in) :: z, gamma, x, mu
      double precision, intent(in), dimension(:) :: t
      double precision, dimension(size(t)) :: tobs
      tobs = (1d0 + z) * ((t / Doppler(gamma, mu)) + (gamma * x * (bofg(gamma) - mu)) / cspeed)
   end function t_obs_v


   !
   !     ::::::  Time in the co-moving frame  ::::::
   !
   ! NOTE: Argument 'x' is the position in the comoving frame
   !
   function t_com_s(t, z, gamma, x, mu) result(tcom)
      implicit none
      double precision, intent(in) :: t, z, gamma, x, mu
      double precision :: tcom
      tcom = Doppler(gamma, mu) * ((t / (1d0 + z)) + (gamma * x * (mu - bofg(gamma))) / cspeed)
   end function t_com_s

   function t_com_v(t, z, gamma, x, mu) result(tcom)
      implicit none
      double precision, intent(in) :: z,gamma, x, mu
      double precision, intent(in), dimension(:) :: t
      double precision, dimension(size(t)) :: tcom
      tcom = Doppler(gamma, mu) * ((t / (1d0 + z)) + (gamma * x * (mu - bofg(gamma))) / cspeed)
   end function t_com_v


   !
   ! #####  #  ####  ##### #####  # #####  #    # ##### #  ####  #    #  ####
   ! #    # # #        #   #    # # #    # #    #   #   # #    # ##   # #
   ! #    # #  ####    #   #    # # #####  #    #   #   # #    # # #  #  ####
   ! #    # #      #   #   #####  # #    # #    #   #   # #    # #  # #      #
   ! #    # # #    #   #   #   #  # #    # #    #   #   # #    # #   ## #    #
   ! #####  #  ####    #   #    # # #####   ####    #   #  ####  #    #  ####
   !
   
   ! :::::  Relativistic Maxwell distribution  :::::
   function RMaxwell_s(g, th) result(rm)
      implicit none
      double precision, intent(in) :: g, th
      double precision :: rm
      
      rm = dmax1(1d-200, g**2 * bofg(g) * dexp(-g / th) / (th * dexp(K2_func(-dlog(th)))))
      
      return
   end function RMaxwell_s

   function RMaxwell_v(g, th) result(rm)
      implicit none
      doubleprecision, intent(in) :: th
      double precision, intent(in), dimension(:) :: g
      integer :: k
      double precision, dimension(size(g)) :: rm
      
      do k=1,size(g)
         rm(k) = dmax1( 1d-200, g(k)**2 * bofg(g(k)) * dexp(-g(k) / th) &
         / (th * dexp(K2_func(-dlog(th)))) )
      enddo
      
      return
   end function RMaxwell_v
   
   ! :::::  Power-law distribution  :::::
   function powlaw_dis_s(g,g1,g2,s) result(pwl)
      implicit none
      doubleprecision, intent(in) :: g,g1,g2,s
      doubleprecision :: pwl
      if ( g >= g1 .and. g <=g2 ) then
         pwl = g**(-s)
      else
         pwl = 0d0
      end if
   end function powlaw_dis_s

   function powlaw_dis_v(g,g1,g2,s) result(pwl)
      implicit none
      doubleprecision, intent(in) :: g1,g2,s
      doubleprecision, intent(in), dimension(:) :: g
      integer :: k
      doubleprecision, dimension(size(g)) :: pwl
      do k = 1, size(g)
         if ( g(k) >= g1 .and. g(k) <= g2 ) then
            pwl(k) = g(k)**(-s)
         else
            pwl(k) = 0d0
         end if
      end do
   end function powlaw_dis_v
   
#if 0
   !
   !     Slope of the Relativistic Maxwell distribution
   !
   function dRMaxwell(g, th) result(drm)
      use K2
      implicit none
      double precision, intent(in) :: g, th
      double precision :: drm
      
      !!$drm = -((2d0 * g**(2) - 1d0) / (g**(2) - 1d0) - g / th)
      drm = (g - g**3 - th + 2 * th * g**2) * dexp(-g/th) / ( dsqrt(g**2 - 1d0) * th**2 * dexp(K2_func(-dlog(th))) )
      
      return
   end function dRMaxwell
#endif

end module SRtoolkit
