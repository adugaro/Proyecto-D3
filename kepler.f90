module kepler

  use iso_fortran_env, only :wp => REAL64

contains


  subroutine solkep_pablo(ex,m,e)

    ! ==================================================================

    implicit none

    ! --------------------------------------------------      
    ! Variables globales
    ! --------------------------------------------------

    real(wp), intent(in) :: ex
    real(wp), intent(inout) :: m
    real(wp), intent(out) :: e

    ! --------------------------------------------------      
    ! Variables locales
    ! --------------------------------------------------

    real(wp) :: ce
    real(wp) :: dex   
    real(wp) :: dpi
    real(wp) :: e0
    real(wp) :: ec
    real(wp) :: es
    real(wp) :: mk 
    real(wp) :: se
    real(wp) :: u   
    real(wp) :: xpri
    real(wp) :: xseg 
    integer :: MAXITER1       ! Número máx de iteraciones, primer  metodo
    integer :: niter
    real(wp) :: tol  ! Tolerancia para la solucion     
    real(wp) :: zero
    real(wp) :: doble
    real(wp) :: ocho      

    parameter (MAXITER1 =     20_wp)  
    parameter ( tol = 1.0E-14_wp )  
    parameter (zero     =  0.0_wp)
    parameter (doble    =  2.0_wp)
    parameter (ocho     =  8.0_wp)

    ! ==================================================================

    ! -----------------------------------------------------------------
    ! Bloque de procesamiento
    ! -----------------------------------------------------------------

    dpi   = ocho*atan(1.0_wp)  
    m = m - int(m/dpi) * dpi
    e     = m  
    do niter = 1, MAXITER1
       e0   = e  
       se   = sin(e0)  
       ce   = cos(e0)  
       es   = ex*se  
       ec   = 1.0_wp-ex*ce  
       mk   = e0-es  
       u    = (mk-m)/ec  
       xpri = e0-u  
       xseg = e0-u/(1.0_wp-u*es)  
       e    = (xpri+xseg)/doble
       dex  = abs(e-e0)
       if(dex.le.tol) return
    end do
    write(*,*) "No alcanzo las iteracion MAXITER1=20"

  end subroutine solkep_pablo

!!$  ! ---------------------------------------------------------------------

  subroutine solkep_simpson (ex,m,ea)

    implicit none

    ! --------------------------------------------------      
    ! Variables globales
    ! --------------------------------------------------

    real(wp), intent(in) :: ex
    real(wp), intent(inout) :: m
    real(wp), intent(out) :: ea
    ! --------------------------------------------------      
    ! Variables locales
    ! --------------------------------------------------
    real(wp) :: pi 
    real(wp) :: twopi
    real(wp) :: halfpi
    real(wp) :: quarterpi
    real(wp) :: m1
    real(wp) :: d

    integer :: f,j

    twopi = 8.0_wp*atan(1.0_wp)
    pi = twopi/2.0_wp
    halfpi = twopi/4.0_wp
    quarterpi = pi / 4.0_wp

    ! ==================================================================

    ! -----------------------------------------------------------------
    ! Bloque de procesamiento
    ! -----------------------------------------------------------------

    f = SGN(m)
    m = abs(m)/twopi
    m = (m-int(m))*twopi*f
    if (m < 0.0_wp) m = m + twopi
    f = 1
    if (m > pi) then
       f = -1
       m = twopi - m
    endif
    ea = halfpi
    d = quarterpi
    do j = 1, 20
       m1 = ea - ex * sin(ea)
       ea = ea + d * SGN(m-m1)
       d = 0.5_wp*d
    enddo

    ea = ea * f
  end subroutine solkep_simpson
  ! ---------------------------------------------------------------------
  function SGN(x) RESULT (y)

    implicit none
    real(wp) :: x
    integer y
    if (x .EQ. 0.0_wp) then
       y = 0
       return
    endif

    if (x .GT. 0.0_wp) then
       y = +1
    else
       y = -1
    endif
    return
  end function SGN

  ! ---------------------------------------------------------------------
!!$*
!!$*                           SUBROUTINE NEWTONM
!!$*
!!$*  This SUBROUTINE performs the Newton Rhapson iteration to find the
!!$*    Eccentric Anomaly given the Mean anomaly.  The True Anomaly is also
!!$*    calcuLated.
!!$*
!!$*  Author        : David Vallado                  303-344-6037    1 Mar 2001
!!$*
!!$*  Inputs          Description                    Range / Units
!!$*    Ecc         - Eccentricity                   0.0D0 to
!!$*    M           - Mean Anomaly                   -2Pi to 2Pi rad
!!$*
!!$*  Outputs       :
!!$*    E0          - Eccentric Anomaly              0.0D0 to 2Pi rad
!!$*    Nu          - True Anomaly                   0.0D0 to 2Pi rad
!!$*
!!$*  Locals        :
!!$*    E1          - Eccentric Anomaly, next value  rad
!!$*    Sinv        - Sine of Nu
!!$*    Cosv        - Cosine of Nu
!!$*    Ktr         - Index
!!$*    R1r         - CUBIC roots - 1 to 3
!!$*    R1i         - imaginary component
!!$*    R2r         -
!!$*    R2i         -
!!$*    R3r         -
!!$*    R3i         -
!!$*    S           - Variables for parabolic solution
!!$*    W           - Variables for parabolic solution
!!$*
!!$*  Coupling      :
!!$*    CUBIC       - Solves a CUBIC polynomial
!!$*    SINH        - Hyperbolic Sine
!!$*    COSH        - Hyperbolic Cosine
!!$*
!!$*  References    :
!!$*    Vallado       2001, 72-75, Alg 2, Ex 2-1
!!$*
!!$* ------------------------------------------------------------------------------
!!$*

  SUBROUTINE solkep_newton(Ecc,M,E0)
    IMPLICIT NONE
    REAL(wp) ::  Ecc, M, E0
    !    * -----------------------------  Locals  ------------------------------
    INTEGER :: Ktr, NumIter
    REAL(wp) :: Small, E1, Pi
    !    * -------------------------  Implementation   -------------------------
    NumIter =    100
    Small   =     0.0000000000001_wp
    Pi      = 3.14159265358979_wp

    ! -------------------- Elliptical ----------------------

    ! -----------  Initial Guess -------------
    IF ( ((M .lt. 0.0_wp) .and. (M .gt. -Pi)) .or. (M .gt. Pi) ) THEN
       E0= M - Ecc
    ELSE
       E0= M + Ecc
    ENDIF
    Ktr= 1
    E1 = E0 + ( M - E0 + Ecc*SIN(E0) )/( 1.0_wp - Ecc*COS(E0) )
    DO WHILE (( ABS(E1-E0) .gt. Small ) .and. ( Ktr .le. NumIter ))
       Ktr = Ktr + 1
       E0= E1
       E1= E0 + ( M - E0 + Ecc*SIN(E0) ) / ( 1.0_wp - Ecc*COS(E0) )
    ENDDO
    ! -------------  Find True Anomaly  ---------------
    !  Sinv= ( DSQRT( 1.0_wp-Ecc*Ecc ) * DSIN(E1) )/( 1.0_wp-Ecc*DCOS(E1) )
    !  Cosv= ( DCOS(E1)-Ecc ) / ( 1.0_wp - Ecc*DCOS(E1) )
    !  Nu  = DATAN2( Sinv,Cosv )
    RETURN
  END SUBROUTINE solkep_newton




end module kepler
