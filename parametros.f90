module parametros
    use iso_fortran_env, only : wp => REAL64
    implicit none
    real(wp), parameter :: Msol = 1.0_wp 
    real(wp), parameter :: G = 2.959122082855911E-4_wp ! G en ua^3/(Msol*dia^2)
end module parametros
