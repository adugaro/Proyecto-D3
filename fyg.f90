module fyg

  use iso_fortran_env, only :wp => REAL64
  use kepler
contains

  ! ------------------------------------------------------------------

  subroutine calculo_fyg (nbig, N_dim,masa,pos,vel,deltat)

    use parametros

    ! ==================================================

    implicit none

    ! --------------------------------------------------
    ! Variables Globales
    ! -------------------------------------------------- 

    integer, intent(in) :: nbig
    integer, intent(in) :: N_dim

    real(wp), dimension(nbig), intent (in) :: masa   
    real(wp), dimension(nbig,N_dim), intent(inout) :: pos
    real(wp), dimension(nbig,N_dim), intent(inout) :: vel
    real(wp), intent(in) :: deltat

    ! --------------------------------------------------
    ! Variables Locales 
    ! --------------------------------------------------

    integer :: k
    integer :: i
    real(wp) :: mu
    real(wp) :: unosobrea
    real(wp) :: n 
    real(wp) :: v0cuad
    real(wp) :: r0mod
    real(wp) :: deltaE 
    real(wp) :: r
    real(wp) :: rsobrea 
    real(wp) :: a
    real(wp) :: asobrer0
    real(wp) :: asobrer
    real(wp) :: unosobren
    real(wp), dimension(N_dim) :: v0
    real(wp), dimension(N_dim) :: r0
    real(wp) :: funF 
    real(wp) :: funG
    real(wp) :: funFp
    real(wp) :: funGp
    real(wp) :: ecosE0
    real(wp) :: esenE0
    real(wp) :: cosdeltaE
    real(wp) :: sendeltaE
    real(wp) :: ex2
    real(wp) :: ex
    real(wp) :: se
    real(wp) :: ce
    real(wp) :: e0
    real(wp) :: m
    real(wp) :: e
    real(wp) :: r0e0

    ! ==================================================

    do k=1,nbig ! Para todos los embriones

       do i=1,N_dim    ! Para las 3 componentes
          r0(i) = pos(k,i)
          v0(i) = vel(k,i)
       end do

       ! Calculo Semieje
       r0mod = sqrt(dot_product(r0,r0))
       v0cuad = dot_product(v0,v0)
       mu = G * MSol
       unosobrea = 2.0_wp/r0mod - v0cuad/mu ! calculo 1/ a
       if (unosobrea < 0.0) then
          write(*,*) "El semieje es negativo"
          stop
       endif
       a = 1.0_wp / unosobrea ! Semieje

       !Calculo Mov. Medio
       n = sqrt(mu*unosobrea**3) ! Mov. medio
       unosobren = 1.0_wp/n 

       ! Calculo Excentricidad
       ecosE0 =  1.0_wp - r0mod * unosobrea  
       esenE0 =  dot_product(r0,v0) * unosobrea**2 * unosobren
       ex2 = ecosE0*ecosE0 + esenE0*esenE0
       ex  = sqrt(ex2) ! Excentricidad

       !   dm = deltat*n- int(deltat*n/dpi)*dpi

       e0  = atan2(esenE0,ecosE0) ! Anom. Excentrica inicial

       m   = e0-esenE0+n*deltat   ! Anom. Media

       ! ================================================= 
       !                   KEPLER  (Conociendo ex y m, te devuelve e)      
                                                         
       call solkep_pablo(ex,m,e)
       !  call solkep_simpson(ex,m,e)
       !  call solkep_newton(ex,m,e)
       ! ==================================================
       
       ! Calculamos F y G (funF y funG)

       r0e0 =  r0mod*ex  
       se   =  sin(e)  
       ce   =  cos(e)  
       funF    =  a*((ce-ex)*ecosE0+se*esenE0)/r0e0  
       funG    =  ((ecosE0-ex2)*se-(ce-ex)*esenE0)/(n*ex)  

       ! Actualizamos posiciones y calculamos r 

       do i=1,N_dim
          pos(k,i) = funF * r0(i) + funG * v0(i)
       end do
       r = sqrt(pos(k,1)**2+pos(k,2)**2+pos(k,3)**2)

       ! Calculamos F' y G'

       funFp   = -sqrt(a*mu)*(se*ecosE0-ce*esenE0)/(r*r0e0)  
       funGp   =  (1.0_wp+funG*funFp)/funF 

       ! Actualizamos velocidades

       do i=1,N_dim
          vel(k,i) = funFp * r0(i) + funGp * v0(i)
       end do

    end do ! do k=1, nbig

  end subroutine calculo_fyg

  ! ---------------------------------------------------------------------

end module fyg
