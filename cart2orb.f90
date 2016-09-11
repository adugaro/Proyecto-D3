program cartesianas_a_orbitales

  use iso_fortran_env, only :wp => REAL64
  
  implicit none

  real(wp), allocatable :: masas(:)
  real(wp), allocatable :: masa(:)
  real(wp), allocatable :: pos(:,:)
  real(wp), allocatable :: vel(:,:)
  real(wp), allocatable :: x(:)
  real(wp), allocatable :: v(:)
  real(wp), parameter :: Msol = 1.0_wp 
  real(wp), parameter ::  G = 2.959122082855911E-4_wp ! G en ua^3/(Msol*dia^2)
  real(wp) :: GGM,t
  real(wp) :: Q,E,FI,FN,PER
  integer :: nbig
  integer :: N_dim
  integer :: i
  integer :: k
  integer :: unidad
  integer :: nfilas
  integer :: planeta
  character(80)::archivo_cond_ini 
  character(80)::archivo_salida

  !------------------------------------------------------------------

  planeta = 1 
  archivo_cond_ini = "salida_venus.sal"
  archivo_salida = "venus_orb.sal"

  nbig = 1
  N_dim=3

  allocate(pos(nbig,N_dim))
  allocate(vel(nbig,N_dim))
  allocate(x(N_dim))
  allocate(v(N_dim))
  allocate (masas(9))
  allocate(masa(1))
  masas(1) = 1.66013679527193009E-07_wp ! Mercurio
  masas(2) = 2.44783833966454430E-06_wp ! Venus
  masas(3) = 3.04043264264672381E-06_wp ! Tierra
  masas(4) = 3.22715144505386530E-07_wp ! Marte
  masas(5) = 9.54791938424326609E-04_wp ! Jupiter
  masas(6) = 2.85885980666130812E-04_wp ! Saturno
  masas(7) = 4.36624404335156298E-05_wp ! Urano
  masas(8) = 5.15138902046611451E-05_wp ! Neptuno
  masas(9) = 7.39644970414201173E-09_wp ! Pluton

  masa(1) = masas(planeta)
  GGM= G*(Msol+masa(1))
  
  !------------------------------------------------------------------

  CALL calcular_numero_lineas(archivo_cond_ini, nfilas)

  open(unit=10,file=archivo_cond_ini)
  open(unit=11,file=archivo_salida)

  do i=1,nfilas

     read(10,*) t,pos(1,1),pos(1,2),pos(1,3),vel(1,1),vel(1,2),vel(1,3)

     call demo2helio(pos, vel, masa, nbig, N_dim)

     x(1) = pos(1,1) 
     x(2) = pos(1,2) 
     x(3) = pos(1,3) 
     v(1) = vel(1,1) 
     v(2) = vel(1,2) 
     v(3) = vel(1,3) 

  !   write(*,*) x(1),x(2),x(3)

     call  ELEM(GGM,x,v,Q,E,FI,FN,PER)

     write(11,*) t,Q,E,FI,FN,PER

  enddo

  close(10)
  close(11)

end program cartesianas_a_orbitales

!------------------------------------------------------------------

subroutine demo2helio(pos,vel,masa,nbig,N_dim)
  use iso_fortran_env, only :wp => REAL64
  implicit none
  !--------------------
  integer, intent(in) :: nbig
  integer, intent(in) :: N_dim
  real(wp), dimension(nbig), intent(in) :: masa
  real(wp), dimension(nbig, N_dim), intent(inout) :: pos
  real(wp), dimension(nbig, N_dim), intent(inout) :: vel
  !--------------------
  real(wp), dimension(N_dim) :: mvsuma
  real(wp) :: temp,msol
  integer :: i
  integer :: k

  msol = 1.0_wp

  mvsuma = 0.0_wp
  do k = 1, nbig
     do i=1,N_dim
        mvsuma(i) = mvsuma(i) + masa(k)*vel(k,i)          
     end do
  end do

  temp = 1.0_wp / msol
  mvsuma(1) = temp * mvsuma(1)
  mvsuma(2) = temp * mvsuma(2)
  mvsuma(3) = temp * mvsuma(3)

  do k = 1, nbig
     do i=1,N_dim
        vel(k,i) = vel(k,i) + mvsuma(i)
     end do
  end do

end subroutine demo2helio

  !--------------------------------------------------
  SUBROUTINE calcular_numero_lineas(inputfile, N)
    IMPLICIT NONE

    CHARACTER(80), INTENT(in) :: inputfile
    INTEGER, INTENT(out) :: N
    INTEGER :: io

    N=0
    OPEN(UNIT=10, FILE=inputfile, ACTION='READ')
    DO
       READ(10, *, IOSTAT=io)
       IF (io /= 0) EXIT
       N = N + 1
    END DO
    CLOSE(10)
    RETURN
  END SUBROUTINE calcular_numero_lineas

  !--------------------------------------------------

  !------------------------------------------------------------------

  SUBROUTINE ELEM(GGM,X,V,Q,E,FI,FN,PER)
    use iso_fortran_env, only :wp => REAL64

    !   CC PARAMETERS :  
    !   CC         IN :  GM    : GRAVITY-CONSTANT (GAUSSIAN CONSTANT R*8  
    !   CC                       PLANETARY APPLICATIONS)  
    !   CC                       IN  (UA**1/3)  / DAYS * SOLARMASS**1/2  
    !   CC               T     : EPOCH TO WHICH X,V (AND ELEMENTS)   R*8  
    !   CC                       REFER (MJD)  
    !   CC               X,V   : POSITION AND VELOCITY OF CELESTIAL R*8(3)  
    !   CC                       BODY AT TIME T  
    !   CC        OUT :  Q      : PERICENTRIC DISTANCE                R*8  
    !   CC               E      : NUMERICAL EXCENTRICITY              R*8  
    !   CC               FI     : INCLINATION OF ORBITAL PLANE WITH   R*8  
    !   CC                        RESPECT TO FUNDAMENTAL PLANE OF  
    !   CC                        COORDINATE SYSTEM  
    !   CC               FN     : LONGITUDE OF ASCENDING NODE         R*8  
    !   CC               PER    : ARGUMENT OF PERIGEE/HELION          R*8  
    !   CC               T0     : TIME OF PERIGEE/HELION PASSING      R*8  

    implicit none

    real(wp), dimension(3) :: H
    real(wp), dimension(3) :: v
    real(wp), dimension(3) :: x
    real(wp), dimension(3) :: xx
    real(wp) :: fn,ck,sk,fi,ci,si
    real(wp) :: u,p,r,ve,rrp,ct,alfa,e
    real(wp) :: a,cex,sex,ex,v1,per,q,ggm

    ! momento angular
    H(1)= X(2)*V(3)-X(3)*V(2)  
    H(2)=-X(3)*V(1)+X(1)*V(3)  
    H(3)= X(1)*V(2)-X(2)*V(1)

    !C-- LONG. OF NODE ,INCLINATION AND ARG OF LATITUDE  
    FN=ATAN2(H(1),H(2))  
    CK=COS(FN)  
    SK=SIN(FN)  
    FI=ATAN2(SQRT(H(1)**2+H(2)**2),H(3))  
    CI=COS(FI)  
    SI=SIN(FI)

    XX(1)=    CK*X(1)   +SK*X(2)  
    XX(2)=-CI*SK*X(1)+CI*CK*X(2)+SI*X(3)  
    XX(3)= SI*SK*X(1)-SI*CK*X(2)+CI*X(3)  
    U=ATAN2(XX(2),XX(1))  

    !C-- P/GM, R, V Y COS(ANGULO ENTRE R Y V)  
    P=(H(1)**2+H(2)**2+H(3)**2)/GGM  
    R= SQRT(X(1)**2+X(2)**2+X(3)**2)  
    VE=SQRT(V(1)**2+V(2)**2+V(3)**2)  
    RRP=X(1)*V(1)+X(2)*V(2)+X(3)*V(3)  
    CT=RRP/(R*VE)  
    !C-- ENERGIA POR UNIDAD DE MASA,Y EXCENTRICIDAD  
    ALFA= VE * VE / 2.0_wp - GGM / R  
    E= SQRT(1.0_wp + 2.0_wp*P*ALFA/GGM)  

    !C-- CASO ELIPTICO  
    A= P/(1.0_wp - E*E )  
    CEX= (A-R) / (A*E)  
    SEX= RRP/(E*SQRT(GGM*A))  
    EX= ATAN2(SEX,CEX)  

    V1= 2.0_wp * ATAN( SQRT((1.0_wp+E)/(1.0_wp-E)) * TAN(EX/2.0_wp) )  
    PER= U-V1  
    Q= A*(1.0_wp-E)  
    RETURN

  end SUBROUTINE ELEM

