module rutinas2
  use iso_fortran_env, only :wp => REAL64

contains

  ! ------------------------------------------------------------------

  subroutine leer_par_sist(nbig,nsmall,nbody,N_dim,t_inicial,t_final,deltat,tescr,deltatescr,nbig_cond_ini,nsmall_cond_ini)

    implicit none

    integer, intent(out) :: nbig
    integer, intent(out) :: nsmall
    integer, intent(out) :: nbody
    integer, intent(out) :: N_dim

    real(wp), intent(out) :: t_inicial
    real(wp), intent(out) :: t_final
    real(wp), intent(out) :: deltat
    real(wp), intent(out) :: tescr
    real(wp), intent(out) :: deltatescr    

    character(80), intent(out) :: nbig_cond_ini 
    character(80), intent(out) :: nsmall_cond_ini

    OPEN(10,FILE='param.en',STATUS='old')

    READ(10,*) 
    READ(10,*) nbig
    READ(10,*) 
    READ(10,*) nsmall
    READ(10,*)
    READ(10,*) N_dim
    READ(10,*)
    READ(10,*) t_inicial
    READ(10,*)
    READ(10,*) t_final
    READ(10,*)
    READ(10,*) deltat
    READ(10,*)
    READ(10,*) tescr
    READ(10,*)
    READ(10,*) deltatescr
    READ(10,*)
    READ(10,*) nbig_cond_ini
    READ(10,*)
    READ(10,*) nsmall_cond_ini

    CLOSE(10)

    nbody = nbig + nsmall


  end subroutine leer_par_sist

  ! ------------------------------------------------------------

  subroutine leer_cond_ini (nbig_cond_ini,nsmall_cond_ini,nbig,nsmall,nbody,N_dim,pos,vel,masa,radio_hill,rho)

    implicit none

    ! Entrada

    integer, intent(in) :: nbig
    integer, intent(in) :: nsmall
    integer, intent(in) :: nbody
    integer, intent(in) :: N_dim

    character(80), intent(in) :: nbig_cond_ini 
    character(80), intent(in) :: nsmall_cond_ini

    ! Salida

    real(wp),dimension(nbody), intent(out) :: masa
    real(wp),dimension(nbody), intent(out) :: radio_hill
    real(wp),dimension(nbody), intent(out) :: rho
    real(wp),dimension(nbody,N_dim), intent(out) :: pos
    real(wp),dimension(nbody,N_dim), intent(out) :: vel

    ! Variables locales

    integer :: i
    integer :: j

    ! ------------------------------------------------------------

    ! C.I. para los embriones 

    open(unit=10,file=nbig_cond_ini,status='old')

    do i=1,nbig
       read(10,'(14x,f23.17,3x,f5.0,3x,f4.2)') masa(i),radio_hill(i),rho(i)
       read(10,*) (pos(i,j),j=1,N_dim)
       read(10,*) (vel(i,j),j=1,N_dim)
       read(10,*)
    enddo

    close(10)

    ! -------------

    ! C.I. para los cuerpos menores

    if (nsmall > 0) then

       open(unit=10,file=nsmall_cond_ini,status='old')

       do i=1,nsmall
          read(10,'(14x,f23.17,3x,f5.0,3x,f4.2)') masa(nbig+i),radio_hill(nbig+i),rho(nbig+i)
          read(10,*) (pos(nbig+i,j),j=1,N_dim)
          read(10,*) (vel(nbig+i,j),j=1,N_dim)
          read(10,*)
       enddo

       close(10)

    endif

    ! ------------------------------------------------------------

  end subroutine leer_cond_ini

  ! -------------

end module rutinas2
