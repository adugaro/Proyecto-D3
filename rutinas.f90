module rutinas
  use iso_fortran_env, only :wp => REAL64
  use parametros

contains

  ! ==================================================================
  
  subroutine helio2demo(pos,vel,masa,nbig,N_dim)

    implicit none

    !--------------------
    integer, intent(in) :: nbig
    integer, intent(in) :: N_dim
    real(wp), dimension(nbig), intent(in) :: masa
    real(wp), dimension(nbig, N_dim), intent(inout) :: pos
    real(wp), dimension(nbig, N_dim), intent(inout) :: vel
    !--------------------
    real(wp), dimension(N_dim) :: mvsuma
    real(wp) :: masasuma
    real(wp) :: temp
    integer :: i
    integer :: k

    mvsuma = 0.0_wp
    masasuma = 0.0_wp

    do k = 1, nbig
       do i=1,N_dim
          mvsuma(i) = mvsuma(i) + masa(k)*vel(k,i)          
       end do
       masasuma = masasuma + masa(k)
    end do

    temp = 1.0_wp / (msol + masasuma)
    mvsuma(1) = temp * mvsuma(1)
    mvsuma(2) = temp * mvsuma(2)
    mvsuma(3) = temp * mvsuma(3)

    do k = 1, nbig
       do i=1,N_dim
          vel(k,i) = vel(k,i) - mvsuma(i)
       end do
    end do

  end subroutine helio2demo

  !------------------------------------------------------------------  

  subroutine demo2helio(pos,vel,masa,nbig,N_dim)

    implicit none
    !--------------------
    integer, intent(in) :: nbig
    integer, intent(in) :: N_dim
    real(wp), dimension(nbig), intent(in) :: masa
    real(wp), dimension(nbig, N_dim), intent(inout) :: pos
    real(wp), dimension(nbig, N_dim), intent(inout) :: vel
    !--------------------
    real(wp), dimension(N_dim) :: mvsuma
    real(wp) :: temp
    integer :: i
    integer :: k

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

  ! ------------------------------------------------------------------
  real(wp) function calcular_energia(pos, vel, masa, nbig, N_dim)
    implicit none

    integer, intent(in) :: nbig
    integer, intent(in) :: N_dim
    real(wp), dimension(nbig, N_dim), intent(in) :: pos
    real(wp), dimension(nbig, N_dim), intent(in) :: vel
    real(wp), dimension(nbig), intent(in) :: masa

    !--------------------
    integer :: i
    integer :: k

    real(wp) :: mvsuma
    real(wp) :: enerk
    real(wp) :: enerp
    real(wp), dimension(N_dim) :: velsol

    do i=1,N_dim
       mvsuma = 0.0_wp
       do k=1,nbig
          mvsuma = mvsuma + masa(k)*vel(k,i)
       end do
       velsol(i) = -mvsuma/Msol
    end do

    enerk= Msol*(velsol(1)**2+velsol(2)**2+velsol(3)**2)
    !    enerk = 0._wp

    do k=1, nbig
       enerk = enerk + masa(k)*(vel(k,1)**2+vel(k,2)**2+vel(k,3)**2)
    end do
    enerk = 0.5_wp*enerk

    enerp = 0.0_wp
    do k=1,nbig        
       enerp = enerp - masa(k)/(sqrt(pos(k,1)**2+pos(k,2)**2+pos(k,3)**2)) 
    end do
    enerp = enerp*G*Msol

    calcular_energia = enerk + enerp
    return
  end function calcular_energia
  !--------------------------------------------------

end module rutinas
