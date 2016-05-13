! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module sampling_handle_module

use constants_module
use error_module
use parser_module
use heap_functions_module

use simulation_setup_module
use lattice_setup_module
use mechanism_setup_module

use kmc_simulation_handle_module
use lattice_handle_module
use rates_handle_module

implicit none

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine obtain_elemstep_statistics(mproc)

implicit none

integer jelemstep, mproc
integer(8) Ntmp

real(8) taunew

if (mproc == 0) then
    return
endif

jelemstep = proctypesites(mproc,0)
taunew = curtime-prevtime

! Overall number of events and average event time
elemstep_noccur(0) = elemstep_noccur(0) + 1
elemstep_avgtime(0) = curtime/dble(elemstep_noccur(0))

Ntmp = elemstep_noccur(jelemstep)

! Per elem-step number of events and average event time
elemstep_avgtime(jelemstep) = dble(Ntmp)/dble(Ntmp+1)*elemstep_avgtime(jelemstep) + taunew/dble(Ntmp+1)
elemstep_noccur(jelemstep) = Ntmp + 1

return

end subroutine obtain_elemstep_statistics

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module sampling_handle_module
