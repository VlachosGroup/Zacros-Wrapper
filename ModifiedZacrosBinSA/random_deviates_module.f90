! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module random_deviates_module

use mt19937_module

implicit none

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function uniform_int_dev(nmax)

implicit none

integer nmax

! return an integer random number (uniformly distributed from 1 to Nmax)

uniform_int_dev = int(nmax*uniform_filt_dev()) + 1

end function uniform_int_dev

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

real(8) function uniform_filt_dev()

implicit none

! Generates uniform random number but filters-discards
! the "singular" values 0 and 1.
! These values have zero Lebeque measure so they can be safely 
! discarded without changing the statistics. 
! Algorithmically we need to avoid these values because they 
! create singularity issues

do 
    uniform_filt_dev = uniform_dev()
    if (uniform_filt_dev /= 0.d0 .and. uniform_filt_dev /= 1.d0) then
        exit
    endif
enddo

return
end function uniform_filt_dev

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

real(8) function expon_dev(lamda)

implicit none

real(8) lamda

if (lamda <= tiny(1.d0)) then
    expon_dev = huge(1.d0)
else
    expon_dev = -dlog(uniform_filt_dev())/lamda
endif

return
end function expon_dev

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

real(8) function expon_dev_from_uniform(lamda,uniform_filt_val)

implicit none

real(8) lamda
real(8) uniform_filt_val

if (lamda <= tiny(1.d0)) then
    expon_dev_from_uniform = huge(1.d0)
else
    expon_dev_from_uniform = -dlog(uniform_filt_val)/lamda
endif

return
end function expon_dev_from_uniform

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module random_deviates_module
