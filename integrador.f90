program n_cuerpos
  use iso_fortran_env, only :wp => REAL64
  use rutinas

  ! ==================================================================

  implicit none

  real(wp), allocatable :: masa(:)
  real(wp), allocatable :: pos(:,:)
  real(wp), allocatable :: vel(:,:)
  real(wp), allocatable :: radio_hill(:)
  real(wp), allocatable :: rho(:)
  real(wp), allocatable :: velsol(:)

  real(wp) :: deltat
  real(wp) :: deltatescr
  real(wp) :: t_inicial
  real(wp) :: t_final
  real(wp) :: tescr
  real(wp) :: t
  real(wp) :: energia
  real(wp) :: energia0
  real(wp) :: enerant
  real(wp) :: mvsuma
  real(wp) :: kkk
  integer :: N_embriones
  integer :: N_dim
  integer :: i
  integer :: j
  integer :: k
  integer :: unidad

  character(80)::archivo_cond_ini  

  ! ==================================================================

  ! **************************************************
  ! Lectura de parametros del sistema
  ! **************************************************

  deltat = 6.0_wp 
  deltatescr = 365.25E3_wp
  T_inicial = 0.0_wp 
  T_final =  365.25E5_wp
  tescr = 0.0_wp

  N_embriones= 1
  archivo_cond_ini = "big.in.mercurio"
  N_dim=3

  allocate(masa(N_embriones))
  allocate(pos(N_embriones,N_dim))
  allocate(vel(N_embriones,N_dim))
  allocate(radio_hill(N_embriones))
  allocate(rho(N_embriones))
  allocate(velsol(N_dim))

  ! **************************************************
  ! Lectura de condiciones iniciales
  ! **************************************************

  open(unit=10,file=archivo_cond_ini,status='old')

  do i=1,N_embriones
     read(10,'(14x,f23.17,3x,f5.0,3x,f4.2)') masa(i),radio_hill(i),rho(i)
     read(10,*) (pos(i,j),j=1,N_dim)
     read(10,*) (vel(i,j),j=1,N_dim)
     read(10,*)
  enddo

  close(10)


  ! **************************************************
  ! Convertimos las velocidades de heliocentricas a
  ! a baricentricas, ya que el metodo utiliza coordenadas
  ! democraticas
  ! **************************************************

  call helio2demo(pos, vel, masa, N_embriones, N_dim)

  ! ******************************************************************

  open(unit=31, file='salida_mercurio.sal')
  open(unit=32, file='salida_venus.sal')
  open(unit=33, file='salida_tierra.sal')
  open(unit=34, file='salida_marte.sal')
  open(unit=35, file='salida_jupiter.sal')
  open(unit=36, file='salida_saturno.sal')
  open(unit=37, file='salida_urano.sal')
  open(unit=38, file='salida_neptuno.sal')
  open(unit=39, file='salida_plutos.sal')

  ! ******************************************************************

  ! **************************************************
  ! Bloque de Integracion
  ! **************************************************

  t = 0.0
  do while (t<t_final)

     call integracion(t,pos,vel,deltat,masa,N_embriones,N_dim)

     t = t + deltat
     tescr = tescr + deltat

     if (tescr > deltatescr) then
        do k=1,N_embriones
           unidad=30+k
           write(unidad,*) t, (pos(k,i),i=1,3), (vel(k,i),i=1,3), masa(k)
        end do
        tescr = tescr-deltatescr
     end if

  end do ! termina el do while (t < t_final)

  ! ******************************************************************

  close(31)
  close(32)
  close(33)
  close(34)
  close(35)
  close(36)
  close(37)
  close(38)
  close(39)

  ! ******************************************************************

  deallocate(masa)
  deallocate(pos)
  deallocate(vel)
  deallocate(radio_hill)
  deallocate(rho)

  ! ******************************************************************

end program n_cuerpos
