! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module gauss_quadrature_module

use constants_module, only: iwrite, pi
    
implicit none

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine gauss_quad_data(ngp,gp,gw)

implicit none

integer ngp, nmaxiter, niter, i, j, k
real(8) err1, err2, etol1, etol2, gp(ngp), gw(ngp)
real(8) xleft,xright,scur,snext
real(8) plegendrecur,plegendreprev1,plegendreprev2,dplegendrecur
logical foundsolution 

! Calculates Gauss points and weights for the Gauss-Legendre quadrature
! in the interval [0,1]
xleft = 0.d0
xright = 1.d0

! Maximum number of iterations in Newton's method
nmaxiter = 1000
! Error tolerances
etol1 = 1.d-12 ! for deltaX
etol2 = 1.d-12 ! for right hand side

! Since the roots are symmetric we need to find half of them
do i = 1,(ngp+1)/2
    
    ! Initial guess for Newton's method
    scur = dcos(pi*(i-1.d0/4.d0)/(ngp+1.d0/2.d0))
    
    niter = 0
    foundsolution = .false.
    do j = 1,nmaxiter
                
        plegendreprev1 = 1.d0
        plegendrecur = scur
  
        ! Evaluation of the Legendre polynomial at scur using Bonnet's recursion formula
        do k = 1,ngp-1
            plegendreprev2 = plegendreprev1
            plegendreprev1 = plegendrecur
            plegendrecur = (2.d0*k+1.d0)/(k+1.d0)*scur*plegendreprev1 - k/(k+1.d0)*plegendreprev2
            !continue
        enddo
        
        ! Evaluation of the Legendre polynomial derivative
        dplegendrecur = ngp*scur/(scur**2-1)*plegendrecur - ngp/(scur**2-1)*plegendreprev1
        
        ! Newton iteration
        snext = scur - plegendrecur/dplegendrecur
        
        err1 = dabs(snext-scur)
        err2 = plegendrecur
        !write(*,*) j,err1,err2

        if (err1 < etol1 .and. err2 < etol2) then
            foundsolution = .true.
            exit
        endif
        
        scur = snext
        
    enddo
    
    if (.not.foundsolution) then        
        write(iwrite,*) 'WARNING: Gauss point could not be computed at the required precision'
    endif
        
    ! Compute one Gauss point and its weight at the right half of the interval
    gp(ngp+1-i) = xleft + (xright-xleft)*(scur+1.d0)/2.d0
    gw(ngp+1-i) = (xright-xleft)/((1.d0-scur**2)*dplegendrecur**2)
    
    ! Compute another Gauss point and its weight at the left half of the interval
    gp(i) = xright - (xright-xleft)*(scur+1.d0)/2.d0
    gw(i) = gw(ngp+1-i)
        
enddo

!write(iwrite,*) 'Gauss points & weights found for ngp = ', ngp
!
!do i = 1,ngp
!    write(iwrite,*) gp(i), gw(i)
!enddo

return

end subroutine gauss_quad_data

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module gauss_quadrature_module