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
   use data_types
   use constants
   implicit none

   interface bofg
      module procedure bofg_s
      module procedure bofg_v
   end interface bofg

   interface gofb
      module procedure gofb_s
      module procedure gofb_v
   end interface gofb

   interface pofg
      module procedure pofg_s
      module procedure pofg_v
   end interface pofg

   interface gofp
      module procedure gofp_s
      module procedure gofp_v
   end interface gofp

   interface ip2ofg
      module procedure ip2ofg_s
      module procedure ip2ofg_v
   end interface ip2ofg

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

   interface x_com_f
      module procedure x_com_s
      module procedure x_com_v
   end interface x_com_f

contains
   !
   !     ----------{     beta(gamma)     }----------
   !
   function bofg_s(g) result(b)
      real(dp), intent(in) :: g
      real(dp) :: b
      if ( g <= 1d0 ) then
         b = 0d0
      else
         b = dsqrt(1d0 - 1d0 / g**2)
      end if
   end function bofg_s
   
   function bofg_v(g) result(b)
      real(dp), intent(in), dimension(:) :: g
      integer :: i
      real(dp), dimension(size(g)) :: b
      do i = 1, size(g)
         if ( g(i) <= 1d0 ) then
            b(i) = 0d0
         else
            b(i) = dsqrt(1d0 - 1d0 / g(i)**2)
         end if
      end do
   end function bofg_v
   
   !
   !     ----------{     gamma(beta)     }----------
   !
   function gofb_s(b) result(g)
      real(dp), intent(in) :: b
      real(dp) :: g
      g = 1d0 / dsqrt(1d0 - b**2)
   end function gofb_s

   function gofb_v(b) result(g)
      real(dp), intent(in), dimension(:) :: b
      real(dp), dimension(size(b)) :: g
      g = 1d0 / dsqrt(1d0 - b**2)
   end function gofb_v

   !
   !     ----------{     p(gamma)     }----------
   !
   function pofg_s(g) result(p)
      real(dp), intent(in) :: g
      real(dp) :: p
      p = dsqrt(g**2 - 1d0)
   end function pofg_s

   function pofg_v(g) result(p)
      real(dp), intent(in), dimension(:) :: g
      real(dp), dimension(size(g)) :: p
      p = dsqrt(g**2 - 1d0)
   end function pofg_v

   !
   !     ----------{     gamma(p)     }----------
   !
   function gofp_s(p) result(g)
      real(dp), intent(in) :: p
      real(dp) :: g
      g = dsqrt(p**2 + 1d0)
   end function gofp_s

   function gofp_v(p) result(g)
      real(dp), intent(in), dimension(:) :: p
      real(dp), dimension(size(p)) :: g
      g = dsqrt(p**2 + 1d0)
   end function gofp_v


   !
   !     ----------{     1 / (gamma^2 - 1)     }----------
   !
   function ip2ofg_s(g) result(p)
      real(dp), intent(in) :: g
      real(dp), parameter :: eps = 1d-6
      real(dp) :: p
      if ( (g - 1d0)**3 < 32d0 * eps ) then
         p = 0.5d0 / (g - 1d0) - 0.25d0 + 0.125d0 * (g - 1d0)
      else
         p = 1d0 / (g**2 - 1d0)
      end if
   end function ip2ofg_s

   function ip2ofg_v(g) result(p)
      integer :: i
      real(dp), intent(in), dimension(:) :: g
      real(dp), parameter :: eps = 1d-6
      real(dp), dimension(size(g)) :: p
      do i=1,size(g)
         if ( (g(i) - 1d0)**3 < 32d0 * eps ) then
            p(i) = 0.5d0 / (g(i) - 1d0) - 0.25d0 + 0.125d0 * (g(i) - 1d0)
         else
            p(i) = 1d0 / (g(i)**2 - 1d0)
         end if
      end do
   end function ip2ofg_v


   !
   !     ::::::  Doppler factor  ::::::
   !
   function Doppler(gamma, mu) result(D)
      implicit none
      real(dp), intent(in) :: gamma, mu
      real(dp) :: D
      D = 1d0 / (gamma * (1d0 - bofg(gamma) * mu))
   end function Doppler


   !
   !     ::::::  Viewing angle in the observer frame   ::::::
   !
   function mu_obs_f(gamma, muc) result(muo)
      implicit none
      real(dp), intent(in) :: gamma, muc
      real(dp) :: muo
      muo = (muc + bofg(gamma)) / (1d0 + bofg(gamma) * muc)
   end function mu_obs_f


   !
   !     ::::::  Viewing angle in the comoving frame   ::::::
   !
   function mu_com_f(gamma, muo) result(muc)
      implicit none
      real(dp), intent(in) :: gamma, muo
      real(dp) :: muc
      muc = (muo - bofg(gamma)) / (1d0 - bofg(gamma) * muo)
   end function mu_com_f


   !
   !     ::::::  Frequency in the observer frame  ::::::
   !
   function nu_obs_s(nu, z, dopp) result(nuobs)
      implicit none
      real(dp), intent(in) :: nu, z, dopp
      real(dp) ::nuobs
      nuobs = nu * dopp / (1d0 + z)
    end function nu_obs_s

    function nu_obs_v(nu, z, dopp) result(nuobs)
      implicit none
      real(dp), intent(in) :: z, dopp
      real(dp), intent(in), dimension(:) :: nu
      real(dp), dimension(size(nu)) :: nuobs
      nuobs = nu * dopp / (1d0 + z)
    end function nu_obs_v


   !
   !     ::::::  Frequency in the comoving frame  ::::::
   !
   function nu_com_s(nu, z, dopp) result(nucom)
      implicit none
      real(dp), intent(in) :: nu, z, dopp
      real(dp) :: nucom
      nucom = nu * (1d0 + z) / dopp
   end function nu_com_s

    function nu_com_v(nu, z, dopp) result(nucom)
      implicit none
      real(dp), intent(in) :: z, dopp
      real(dp), intent(in), dimension(:) :: nu
      real(dp), dimension(size(nu)) :: nucom
      nucom = nu * (1d0 + z) / dopp
   end function nu_com_v


   !
   !     ::::::  Time in the observer frame  ::::::
   !
   ! NOTE: Argument 'x' is the position in the comoving frame
   !
   function t_obs_s(t, z, gamma, x, muo) result(tobs)
      implicit none
      real(dp), intent(in) :: t, z, gamma, x, muo
      real(dp) :: tobs
      tobs = (1d0 + z) * ((t / Doppler(gamma, muo)) + (gamma * x * (bofg(gamma) - muo)) / cLight)
   end function t_obs_s

   function t_obs_v(t, z, gamma, x, muo) result(tobs)
      implicit none
      real(dp), intent(in) :: z, gamma, x, muo
      real(dp), intent(in), dimension(:) :: t
      real(dp), dimension(size(t)) :: tobs
      tobs = (1d0 + z) * ((t / Doppler(gamma, muo)) + (gamma * x * (bofg(gamma) - muo)) / cLight)
   end function t_obs_v


   !
   !     ::::::  Time in the co-moving frame  ::::::
   !
   ! NOTE: Argument 'x' is the position in the comoving frame
   !
   function t_com_s(tobs, z, gamma, x, muo) result(tcom)
      implicit none
      real(dp), intent(in) :: tobs, z, gamma, x, muo
      real(dp) :: tcom
      tcom = Doppler(gamma, muo) * ((tobs / (1d0 + z)) + (gamma * x * (muo - bofg(gamma))) / cLight)
   end function t_com_s

   function t_com_v(t, z, gamma, x, muo) result(tcom)
      implicit none
      real(dp), intent(in) :: z,gamma, x, muo
      real(dp), intent(in), dimension(:) :: t
      real(dp), dimension(size(t)) :: tcom
      tcom = Doppler(gamma, muo) * ((t / (1d0 + z)) + (gamma * x * (muo - bofg(gamma))) / cLight)
   end function t_com_v


   !
   !     ::::::  Position in the co-moving frame  ::::::
   !
   function x_com_s(t, tobs, z, gamma, muo) result(xcom)
      implicit none
      real(dp), intent(in) :: t, tobs, z, gamma, muo
      real(dp) :: xcom
      xcom = cLight * ( (t / Doppler(gamma, muo)) - (tobs / (1d0 + z)) ) / ( gamma * (muo - bofg(gamma)) )
   end function x_com_s

   function x_com_v(t, tobs, z, gamma, muo) result(xcom)
      implicit none
      real(dp), intent(in) :: z, gamma, muo, tobs
      real(dp), intent(in), dimension(:) :: t
      real(dp), dimension(size(t)) :: xcom
      xcom = cLight * ( (t / Doppler(gamma, muo)) - (tobs / (1d0 + z)) ) / ( gamma * (muo - bofg(gamma)) )
   end function x_com_v

end module SRtoolkit
