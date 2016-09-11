module hybrido
  use iso_fortran_env, only :wp => REAL64
  use parametros

contains

  ! ------------------------------------------------------------------  

  SUBROUTINE integracion(t,pos,vel,deltat,masa,nbig,N_dim)
    use fyg

    ! ==================================================

    implicit none

    ! ------------------------------------------------------------------
    ! Variables Globales
    ! ------------------------------------------------------------------

    integer, intent(in) :: nbig
    integer, intent(in) :: N_dim    
    real(wp), intent(in) :: deltat
    real(wp), intent(in) :: t
    real(wp),dimension(nbig), intent(in) :: masa
    real(wp),dimension(nbig,N_dim), intent(inout) :: pos
    real(wp),dimension(nbig,N_dim), intent(inout) :: vel

    ! ------------------------------------------------------------------
    !   Variables Locales 
    ! ------------------------------------------------------------------

    real(wp) :: Qki_Qji
    real(wp) :: rij3
    real(wp) :: imp
    real(wp) :: rk3
    integer :: i
    integer :: j
    integer :: k
    real(wp) :: tau2
    real(wp), dimension(3) :: posvect, velvect

    real(wp), dimension(3) :: mvsum
    real(wp) :: temp1
    real(wp) :: zero
    !------------------------------------------------------------------

    tau2=deltat*0.5_wp
    zero= 0.0_wp
    ! ******************************************************************
    ! ******************************************************************

    ! i) Las coordenadas permanecen fijas y cada cuerpo recibe una aceleracion proveniente de los otros cuerpos (pero no del cuerpo principal) que modifica su impulso sobre un intervalo temporal t/2.

    ! k para los embriones target, j para los embriones que interactuan con k
    ! i las componentes x, y, z 

    do k=1,nbig
       do i = 1,N_dim
          Qki_Qji = 0.0_wp
          do j = 1, nbig
             if (k/=j)then
                rij3=sqrt((pos(k,1)-pos(j,1))**2 + (pos(k,2)-pos(j,2))**2 + (pos(k,3)-pos(j,3))**2)**3
                rij3=1.0_wp/rij3
                Qki_Qji = Qki_Qji + masa(j)*(pos(k,i)-pos(j,i))*rij3
             endif
          enddo
          vel(k,i) = vel(k,i) - tau2*G*Qki_Qji
       enddo
    enddo

    ! ******************************************************************
    ! ******************************************************************

    ! ii) Los impulsos permanecen fijos, y cada cuerpo sufre un desplazamiento en su posicion t/2

    mvsum(1) = 0.0_wp
    mvsum(2) = 0.0_wp
    mvsum(3) = 0.0_wp
    do j = 1, nbig
       mvsum(1) = mvsum(1)  +  masa(j) * vel(j,1)
       mvsum(2) = mvsum(2)  +  masa(j) * vel(j,2)
       mvsum(3) = mvsum(3)  +  masa(j) * vel(j,3)
    end do

    temp1 = tau2 / Msol
    mvsum(1) = temp1 * mvsum(1)
    mvsum(2) = temp1 * mvsum(2)
    mvsum(3) = temp1 * mvsum(3)
    do j = 1, nbig
       pos(j,1) = pos(j,1)  +  mvsum(1)
       pos(j,2) = pos(j,2)  +  mvsum(2)
       pos(j,3) = pos(j,3)  +  mvsum(3)
    end do


!!$    temp1 = tau2/Msol
!!$    do k = 1, nbig
!!$       do i =1,N_dim
!!$          imp=0.0_wp
!!$          do j = 1, nbig
!!$             imp = imp + vel(j,i)*masa(j)
!!$          enddo
!!$          pos(k,i) = pos(k,i) + temp1*imp 
!!$       enddo
!!$    enddo
!!$        
!!$     do k =1,3
!!$         mvsum(k) = zero
!!$      enddo     
!!$!C     -----
!!$      do i=1,nbig
!!$       ! j    = ind(i) ! Indice  global
!!$       ! masa = m(j)
!!$        do k=1,3
!!$            mvsum(k) = mvsum(k) + masa(i)*vel(i,k)
!!$        enddo
!!$      enddo    
!!$!C     -----
!!$      do k=1,3
!!$       ! p0(k)    = - summv(k)       
!!$        mvsum(k) = tau2*(mvsum(k)/Msol)
!!$      enddo
!!$!C     -----
!!$      do i=1,nbig
!!$    !    j = ind(i) ! Indice global 
!!$        do k =1,3
!!$            pos(i,k) = pos(i,k) + mvsum(k)
!!$        end do
!!$      enddo


    ! ******************************************************************
    ! ******************************************************************

    ! iii ) Cada cuerpo evoluciona en torno a una orbita kepleriana en un tiempo t
    call calculo_fyg (nbig,N_dim,masa,pos,vel,deltat)

    ! ******************************************************************
    ! ******************************************************************

    ! iv) como en el paso ii
!!$
    mvsum(1) = 0.0_wp
    mvsum(2) = 0.0_wp
    mvsum(3) = 0.0_wp
    do j = 1, nbig
       mvsum(1) = mvsum(1)  +  masa(j) * vel(j,1)
       mvsum(2) = mvsum(2)  +  masa(j) * vel(j,2)
       mvsum(3) = mvsum(3)  +  masa(j) * vel(j,3)
    end do

    temp1 = tau2 / Msol
    mvsum(1) = temp1 * mvsum(1)
    mvsum(2) = temp1 * mvsum(2)
    mvsum(3) = temp1 * mvsum(3)
    do j = 1, nbig
       pos(j,1) = pos(j,1)  +  mvsum(1)
       pos(j,2) = pos(j,2)  +  mvsum(2)
       pos(j,3) = pos(j,3)  +  mvsum(3)
    end do

 !   write(44,*) pos(1,1),pos(1,2),pos(1,3)
!  stop
!!$    
!!$    temp1 = tau2/Msol
!!$    do k = 1, nbig
!!$       do i =1,N_dim
!!$          imp=0.0_wp
!!$          do j = 1, nbig
!!$             imp = imp + vel(j,i)*masa(j)
!!$          enddo
!!$          pos(k,i) = pos(k,i) + (temp1)*imp 
!!$       enddo
!!$    enddo
!!$

!!$     do k =1,3
!!$         mvsum(k) = zero
!!$      enddo     
!!$!C     -----
!!$      do i=1,nbig
!!$       ! j    = ind(i) ! Indice  global
!!$       ! masa = m(j)
!!$        do k=1,3
!!$            mvsum(k) = mvsum(k) + masa(i)*vel(i,k)
!!$        enddo
!!$      enddo    
!!$!C     -----
!!$      do k=1,3
!!$       ! p0(k)    = - summv(k)       
!!$        mvsum(k) = tau2*(mvsum(k)/Msol)
!!$      enddo
!!$!C     -----
!!$      do i=1,nbig
!!$    !    j = ind(i) ! Indice global 
!!$        do k =1,3
!!$            pos(i,k) = pos(i,k) + mvsum(k)
!!$        end do
!!$      enddo


    ! ******************************************************************
    ! ******************************************************************

    ! v) como en el paso i) 

    do k=1,nbig
       do i = 1,N_dim
          Qki_qji = 0.0_wp
          do j = 1, nbig
             if (k/=j)then
                rij3=sqrt((pos(k,1)-pos(j,1))**2 + (pos(k,2)-pos(j,2))**2 + (pos(k,3)-pos(j,3))**2)**3
                rij3=1.0_wp/rij3
                Qki_Qji = Qki_Qji + masa(j)*(pos(k,i)-pos(j,i))*rij3
             endif
          enddo
          vel(k,i) = vel(k,i) - tau2*G*Qki_Qji
       enddo
    enddo

    ! ******************************************************************
    ! ******************************************************************

    ! ------------------------------------------------------------------

  end subroutine integracion


end module hybrido
