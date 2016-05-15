! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module simulation_setup_module
! Sets up simulation.
! 
! .. f:variable:: surfspecdent(:)
!    
!    Number of dentate for a given surface species


use constants_module
use gauss_quadrature_module
use error_module
use parser_module
use mt19937_module

use lattice_setup_module

implicit none

integer iniseq, ngp, kmaxnewton, ngasspecs, nsurfspecs, maxdent, walltime, ieleventoccur
integer(8) maxsteps, curstep, dkeventsnap, dkeventprocstat, dkeventspecnum,snapshnum, procstatnum, specnumnum
integer(8), allocatable :: inewtonstats(:)

integer, allocatable :: surfspecsdent(:)
integer(8), allocatable :: gasspecsnums(:)

real(8) temp, tramp, pres, maxtime, curtime, prevtime, etol1, etol2
real(8) snaptime, dtsnap, procstattime, dtprocstat, specnumtime, dtspecnum

real(8), allocatable :: gasmolfracs(:)
real(8), allocatable :: gasmolweights(:)
real(8), allocatable :: gasenergies(:)
real(8), allocatable :: gp(:)
real(8), allocatable :: gw(:)
real(8), allocatable :: rnewtonstats(:)

character(nnam0), allocatable :: gasspecsnames(:)
character(nnam0), allocatable :: surfspecsnames(:)

logical output_snapshots, snap_on_event, snap_on_logtime
logical output_procstat, procstat_on_event, procstat_on_logtime
logical output_specnum, specnum_on_event, specnum_on_logtime, specnum_on_eleventoccur
logical no_restart, report_events, readwalltime, tpdsim
logical debug_report_processes, debug_check_lattice, debug_check_processes, debug_newtons_method
logical debug_report_globenerg, debug_check_globenerg

contains

  subroutine cleanup_ssm
    ! Cleanup module globals

    if(allocated(inewtonstats))   deallocate(inewtonstats)
    if(allocated(surfspecsdent))  deallocate(surfspecsdent)
    if(allocated(gasspecsnums))   deallocate(gasspecsnums)
    if(allocated(gasmolfracs))    deallocate(gasmolfracs)
    if(allocated(gasmolweights))  deallocate(gasmolweights)
    if(allocated(gasenergies))    deallocate(gasenergies)
    if(allocated(gp))             deallocate(gp)
    if(allocated(gw))             deallocate(gw)
    if(allocated(rnewtonstats))   deallocate(rnewtonstats)
    if(allocated(gasspecsnames))  deallocate(gasspecsnames)
    if(allocated(surfspecsnames)) deallocate(surfspecsnames)

  end subroutine cleanup_ssm

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine read_simulation_setup()

implicit none

integer i, io, irec, nwords, q, setvbuf3f
integer iw(maxwords2)

real(8) sumgasmolarfracs

character(lengthrecinp) recinput

logical readtemp, readpres, readngasspecs, readnsurfspecs, readgasspnames, readgasenergies, readsurfspnames
logical readdentation, readmaxsteps, readmaxtime, readsnapsh, readprocstat, readspecnum
logical readeventreport, readiniseq, readngausspts, readkmaxnewton, readetol1newton, readetol2newton
logical readgasmolfracs, readgasmolweights, genoutop, simfex

! Parses the simulation setup file simulation_input.dat

! Level-1 keywords allowed:
!    random_seed
!    temperature
!    pressure
!    n_gas_species
!    gas_specs_names
!    gas_molec_weights
!    gas_molar_fracs
!    n_surf_species
!    surf_specs_names
!    surf_specs_dent
!    max_steps
!    max_time
!    snapshots
!    process_statistics
!    species_numbers
!    event_report
!    wall_time
!    no_restart
!    debug_report_processes
!    debug_check_processes
!    debug_check_lattice
!    finish

! Check if the lattice specification file exists
simfex = .false.
inquire(file=csimfname,exist=simfex)
if (.not.simfex) then
    moreinfo = 'Simulation specification file ' // trim(csimfname) // ' does not exist.'
    call error(1)
endif

inquire(iwrite,opened=genoutop)
if (.not.genoutop) then
    open(unit=iwrite,file=trim(cgenoutfname),status='unknown',position='append')
    q = setvbuf3f(iwrite,1,0) ! set line-buffered behavior
endif

write(iwrite,'(/,a)') 'Simulation setup:'
write(iwrite,'(a)')   '~~~~~~~~~~~~~~~~~'

open(unit=iread,file=trim(csimfname),status='old')

moreinfo = ''
readiniseq = .false.
readngausspts = .false.
readkmaxnewton = .false.
readetol1newton = .false.
readetol2newton = .false.
readtemp = .false.
readpres = .false.
readngasspecs = .false.
readgasspnames = .false.
readgasenergies = .false.
readgasmolweights = .false.
readgasmolfracs = .false.
readnsurfspecs = .false.
readsurfspnames = .false.
readdentation = .false.
readmaxsteps = .false.
readmaxtime = .false.
readsnapsh = .false.
readprocstat = .false.
readspecnum = .false.
readeventreport = .false.
readwalltime = .false.

output_snapshots = .false.
snap_on_event = .false.
snap_on_logtime = .false.
output_procstat = .false.
procstat_on_event = .false.
procstat_on_logtime = .false.
output_specnum = .false.
specnum_on_event = .false.
specnum_on_logtime = .false.
specnum_on_eleventoccur = .false.
ieleventoccur = 0 ! the elementary event after the occurrence of which the number of species molecules are to be reported
no_restart = .false.
report_events = .false.
debug_report_processes = .false.
debug_report_globenerg = .false.
debug_check_lattice = .false.
debug_check_processes = .false.
debug_check_globenerg = .false.
debug_newtons_method = .false.

maxdent = 0

curtime = 0.d0 ! current time
prevtime = 0.d0 ! previous time (of the last event before the one about to be simulated)
curstep = 0_8 ! current step

irec = 0
io = 0

! Arrays rnewtonstats and inewtonstats give statistics for Newton's method: 
! rnewtonstats(0) average number of Newton iterations
! rnewtonstats(1) max Dx error
! rnewtonstats(2) max RHS error
! inewtonstats(0) number of times Newton's method procedure was executed
! inewtonstats(1) number of times the Newton's method loop failed by exceeding the max number of iterations
allocate(rnewtonstats(0:3))
allocate(inewtonstats(0:1))
rnewtonstats = 0.d0
inewtonstats = 0_8

call getrecord(iread,recinput,irec,io)

do while (io >= 0)

    call break_words(recinput,' ',remchar,nwords,iw)

    if (nwords > 0) then
        
        if (striccompare(recinput(iw(1):iw(2)),'random_seed')) then
            
            if (readiniseq) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif
            
            if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
                iniseq = str2int(recinput(iw(3):iw(4)))

                if (iniseq == 0) then ! seed for random number generator = 0 -> invalid!
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(1102)
                endif
                
                call initialize_rng(iniseq) ! initialize random number generator
                readiniseq = .true.
                write(iwrite,'(/,a)') '    Random sequence with seed: ' // int2str(iniseq)
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1001)                
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'n_gauss_pts')) then
            
            if (readngausspts) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif
            
            if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
                ngp = str2int(recinput(iw(3):iw(4)))
                
                if (ngp <= 0) then
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(1118)
                endif
                
                readngausspts = .true.

                allocate(gp(ngp))
                allocate(gw(ngp))
                
                call gauss_quad_data(ngp,gp,gw) ! initialize Gauss points and weights in interval [0.0,1.0]
                write(iwrite,'(/,a)') '    Number of points in Gauss quadrature: ' // int2str(ngp)
                write(iwrite,'(/,a,' // int2str(ngp) // '(1x,f10.6))') '    Gauss points:  ', (gp(i),i=1,ngp)
                write(iwrite,'(a,' // int2str(ngp) // '(1x,f10.6))') '    Gauss weights: ', (gw(i),i=1,ngp)
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1117)                
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'max_newton_iter')) then
            
            if (readkmaxnewton) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif
            
            if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
                kmaxnewton = str2int(recinput(iw(3):iw(4)))
                
                if (kmaxnewton <= 0) then
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(1120)
                endif
                
                readkmaxnewton = .true.

                write(iwrite,'(/,a)') '    Number of maximum Newton iterations: ' // int2str(kmaxnewton)
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1119)                
            endif
            
        elseif (striccompare(recinput(iw(1):iw(2)),'tol_dx_newton')) then

            if (readetol1newton) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (str_check_expression(recinput,nwords,iw,'AA/R8',.true.)) then
                etol1 = str2dbl(recinput(iw(3):iw(4)))

                if (etol1 < 0) then ! etol1 < 0 -> invalid!
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(1122)
                endif
                
                readetol1newton = .true.
                write(iwrite,'(/,a)') '    Dx tolerance for Newton''s method: ' // dbl2str(etol1)
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1021)                
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'tol_rhs_newton')) then

            if (readetol2newton) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (str_check_expression(recinput,nwords,iw,'AA/R8',.true.)) then
                etol2 = str2dbl(recinput(iw(3):iw(4)))

                if (etol2 < 0) then ! etol2 < 0 -> invalid!
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(1124)
                endif
                
                readetol2newton = .true.
                write(iwrite,'(/,a)') '    RHS tolerance for Newton''s method: ' // dbl2str(etol2)
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1023)                
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'temperature')) then

            if (readtemp) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif
            
            if (str_check_expression(recinput,nwords,iw,'AA/R8',.true.)) then
                temp = str2dbl(recinput(iw(3):iw(4)))

                if (temp <= 0) then  ! temperature <= 0 -> invalid!
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(1104)
                endif
                
                tramp = 0.d0
                tpdsim = .false.
                readtemp = .true.
                write(iwrite,'(/,a)') '    Temperature: ' // trim(dbl2str(temp))
            elseif (str_check_expression(recinput,nwords,iw,'2AA/2R8',.true.) .and. &
                (striccompare(recinput(iw(3):iw(4)),'ramp')) ) then
				temp = str2dbl(recinput(iw(5):iw(6)))

                if (temp <= 0) then  ! temperature <= 0 -> invalid!
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(1104)
                endif

				tramp = str2dbl(recinput(iw(7):iw(8)))
                tpdsim = .true.
				readtemp = .true.
				write(iwrite,'(/,a)') '    Temperature ramp: Initial T: ' // trim(dbl2str(temp)) // &
				                      '. Rate of change: ' // trim(dbl2str(tramp)) // '.'                
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1002)                
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'pressure')) then

            if (readpres) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (str_check_expression(recinput,nwords,iw,'AA/R8',.true.)) then
                pres = str2dbl(recinput(iw(3):iw(4)))

                if (pres < 0) then ! pressure < 0 -> invalid!
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(1106)
                endif
                
                readpres = .true.
                write(iwrite,'(/,a)') '    Pressure: ' // trim(dbl2str(pres))
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1003)                
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'n_gas_species')) then

            if (readngasspecs) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif
            
            if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
                ngasspecs = str2int(recinput(iw(3):iw(4)))

                if (ngasspecs < 0) then ! number of gas species < 0 -> invalid!
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(1108)
                endif
                
                readngasspecs = .true.
                write(iwrite,'(/,a)') '    Number of gas species: ' // trim(int2str(ngasspecs))
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1004)                
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'gas_specs_names')) then

            if (readgasspnames) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (.not.readngasspecs) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call error(1005)
            endif
            
            if (str_check_expression(recinput,nwords,iw,trim(int2str(ngasspecs+1)) // 'AA/',.false.)) then
                if (nwords > ngasspecs + 1) then
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call warning(1001)
                endif
                allocate(gasspecsnames(ngasspecs))
                allocate(gasspecsnums(ngasspecs))
                do i = 1,ngasspecs
                    gasspecsnums(i) = 0_8
                    gasspecsnames(i) = recinput(iw(2*i+1):iw(2*i+2))
                enddo
                readgasspnames = .true.
                write(iwrite,'(/,a)',advance='no') '    Gas species names: '
                do i = 1,ngasspecs
                    write(iwrite,'(a,1x)',advance='no') trim(gasspecsnames(i))
                enddo
                write(iwrite,'(a)') ''
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1006)            
            endif
            
        elseif (striccompare(recinput(iw(1):iw(2)),'gas_molec_weights')) then

            if (readgasmolweights) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (.not.readngasspecs) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call error(1007)
            endif
            
            if (str_check_expression(recinput,nwords,iw,'AA/' // trim(int2str(ngasspecs)) // 'R8/',.false.)) then
                if (nwords > ngasspecs + 1) then
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call warning(1008)
                endif
                allocate(gasmolweights(ngasspecs))
                do i = 1,ngasspecs
                    gasmolweights(i) = str2dbl(recinput(iw(2*i+1):iw(2*i+2)))
                enddo
                
                do i = 1,ngasspecs
                    if (gasmolweights(i) <= 0) then
                        moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                        call error(1110) ! gas species molecular weight < 0 -> invalid!
                    endif
                enddo
                
                readgasmolweights = .true.
                write(iwrite,'(/,a)',advance='no') '    Gas species molecular weights: '
                do i = 1,ngasspecs
                    write(iwrite,'(a,1x)',advance='no') trim(dbl2str(gasmolweights(i)))
                enddo
                write(iwrite,'(a)') ''
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1008)
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'gas_energies')) then

            if (readgasenergies) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (.not.readngasspecs) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call error(1030)
            endif
            
            if (str_check_expression(recinput,nwords,iw,'AA/' // trim(int2str(ngasspecs)) // 'R8/',.false.)) then
                if (nwords > ngasspecs + 1) then
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call warning(1031)
                endif
                allocate(gasenergies(ngasspecs))
                do i = 1,ngasspecs
                    gasenergies(i) = str2dbl(recinput(iw(2*i+1):iw(2*i+2)))
                enddo
                
                readgasenergies = .true.
                write(iwrite,'(/,a)',advance='no') '    Gas species energies: '
                do i = 1,ngasspecs
                    write(iwrite,'(a,1x)',advance='no') trim(dbl2str(gasenergies(i)))
                enddo
                write(iwrite,'(a)') ''
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1031)
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'gas_molar_fracs')) then

            if (readgasmolfracs) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (.not.readngasspecs) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call error(1009)
            endif
            
            if (str_check_expression(recinput,nwords,iw,'AA/' // trim(int2str(ngasspecs)) // 'R8/',.false.)) then
                if (nwords > ngasspecs + 1) then
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call warning(1003)
                endif
                allocate(gasmolfracs(ngasspecs))
                do i = 1,ngasspecs
                    gasmolfracs(i) = str2dbl(recinput(iw(2*i+1):iw(2*i+2)))
                enddo

                do i = 1,ngasspecs
                    if (gasmolfracs(i) < 0.d0 .or. gasmolfracs(i) > 1.d0) then
                        moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                        call error(1112) ! gas species molar fraction < 0 -> invalid!
                    endif
                enddo                

                sumgasmolarfracs = 0.d0
                do i = 1,ngasspecs
                    sumgasmolarfracs = sumgasmolarfracs + gasmolfracs(i)
                enddo
                if (sumgasmolarfracs > 1.d0+1E-8) then !floating point Arithmetic might give a number slightly larger than one
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(1125)
                endif
                
                readgasmolfracs = .true.
                write(iwrite,'(/,a)',advance='no') '    Gas species molar fractions: '
                do i = 1,ngasspecs
                    write(iwrite,'(a,1x)',advance='no') trim(dbl2str(gasmolfracs(i)))
                enddo
                write(iwrite,'(a)') ''
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1010)
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'n_surf_species')) then

            if (readnsurfspecs) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif
            
            if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
                nsurfspecs = str2int(recinput(iw(3):iw(4)))

                if (nsurfspecs <= 0) then ! number of surface species <= 0 invalid
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(1114)
                endif

                readnsurfspecs = .true.
                write(iwrite,'(/,a)') '    Number of surface species: ' // trim(int2str(nsurfspecs))
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1011)
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'surf_specs_names')) then

            if (readsurfspnames) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (.not.readnsurfspecs) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call error(1012)
            endif
            
            if (str_check_expression(recinput,nwords,iw,trim(int2str(nsurfspecs+1)) // 'AA/',.false.)) then
                if (nwords > nsurfspecs + 1) then
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call warning(1004)
                endif
                allocate(surfspecsnames(0:nsurfspecs))
                do i = 1,nsurfspecs
                    surfspecsnames(i) = recinput(iw(2*i+1):iw(2*i+2))
                enddo

                readsurfspnames = .true.
                write(iwrite,'(/,a)',advance='no') '    Surface species names: '
                do i = 1,nsurfspecs
                    write(iwrite,'(a,1x)',advance='no') trim(surfspecsnames(i))
                enddo
                write(iwrite,'(a)') ''
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1013)
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'surf_specs_dent')) then

            if (readdentation) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (.not.readnsurfspecs) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call error(1014)
            endif
            
            if (str_check_expression(recinput,nwords,iw,'AA/' // trim(int2str(nsurfspecs)) // 'I4/',.false.)) then
                if (nwords > nsurfspecs + 1) then
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call warning(1005)
                endif            
                allocate(surfspecsdent(0:nsurfspecs))
                surfspecsdent(0) = 1 ! vacant site
                do i = 1,nsurfspecs
                    surfspecsdent(i) = str2int(recinput(iw(2*i+1):iw(2*i+2)))
                    if (maxdent < surfspecsdent(i)) maxdent = surfspecsdent(i) 
                enddo

                do i = 1,nsurfspecs
                    if (surfspecsdent(i) <= 0) then
                        moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                        call error(1116) ! surface species dentation <= 0 -> invalid!
                    endif
                enddo
                
                readdentation = .true.
                write(iwrite,'(/,a)',advance='no') '    Surface species dentation: '
                write(iwrite,'('//int2str(nsurfspecs)//'(I2,1x))') (surfspecsdent(i),i=1,nsurfspecs)
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1015)
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'max_steps')) then
            
            if (readmaxsteps) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif
            
            if (nwords /= 2) then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1016)
            elseif (striccompare(recinput(iw(3):iw(4)),'infinity')) then
                maxsteps = huge(maxsteps)
                readmaxsteps = .true.
                write(iwrite,'(/,a)') '    Maximum number of steps: ' // trim(int82str(maxsteps)) // ' (maximum allowed value)'
            else
                if (str_is_integer(recinput(iw(3):iw(4)),8)) then
                    maxsteps = str2int8(recinput(iw(3):iw(4)))
                    readmaxsteps = .true.
                    write(iwrite,'(/,a)') '    Maximum number of steps: ' // trim(int82str(maxsteps))
                else
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(1016)                
                endif
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'max_time')) then

            if (readmaxtime) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif
            
            if (nwords /= 2) then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1017)
            elseif (striccompare(recinput(iw(3):iw(4)),'infinity')) then
                maxtime = huge(maxtime)
                write(iwrite,'(/,a)') '    Max simulated time: ' // trim(dbl2str(maxtime)) // ' (maximum allowed value)'
            else
                if (str_is_real(recinput(iw(3):iw(4)),8)) then
                    maxtime = str2dbl(recinput(iw(3):iw(4)))
                    write(iwrite,'(/,a)') '    Max simulated time: ' // dbl2str(maxtime)
                else
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(1017)
                endif
            endif
            readmaxtime = .true.

        elseif (striccompare(recinput(iw(1):iw(2)),'snapshots')) then

            if (readsnapsh) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif
            
            if (nwords == 1 .or. nwords > 5) then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1018)
            endif
            
            if (nwords == 2 .and. striccompare(recinput(iw(3):iw(4)),'off')) then
                ! output_snapshots is already set to .false. upon initialization
                write(iwrite,'(/,a)') '    Snapshot saving turned off'
                
            elseif (striccompare(recinput(iw(3):iw(4)),'on')) then
                
                ! At this point output_snapshots is set to .false. from initialization
                select case (nwords)
                
                    case (3)
                        ! On-event reporting
                        if (striccompare(recinput(iw(5):iw(6)),'event')) then
                            dkeventsnap = 1_8
                            output_snapshots = .true.
                            snap_on_event = .true.
                        endif
                        
                    case (4)
                        ! On-time reporting: need to specify timestep
                        if (striccompare(recinput(iw(5):iw(6)),'time') .and. &
                            str_is_real(recinput(iw(7):iw(8)),8)) then
                            output_snapshots = .true.
                            snap_on_event = .false.
                            snaptime = 0.d0
                            dtsnap = str2dbl(recinput(iw(7):iw(8)))
                        ! On-k^th event reporting: need to specify Deltak
                        elseif (striccompare(recinput(iw(5):iw(6)),'event') .and. &
                            str_is_integer(recinput(iw(7):iw(8)),8)) then
                            dkeventsnap = str2int(recinput(iw(7):iw(8)))
                            output_snapshots = .true.
                            snap_on_event = .true.                                
                        endif
                        
                    case (5)
                        ! On-log-time reporting: need to specify initial time and multiplier
                        if (striccompare(recinput(iw(5):iw(6)),'logtime') .and. &
                            str_is_real(recinput(iw(7):iw(8)),8) .and. &
                            str_is_real(recinput(iw(9):iw(10)),8)) then
                            output_snapshots = .true.
                            snap_on_event = .false.
                            snap_on_logtime = .true.
                            snaptime = str2dbl(recinput(iw(7):iw(8)))
                            dtsnap = str2dbl(recinput(iw(9):iw(10)))
                        endif

                end select

                if (.not. output_snapshots) then
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(1018)
                else
                    open(unit=ihistory,file=trim(chistoryfname),status='unknown')
                    q = setvbuf3f(ihistory,1,0) ! set line-buffered behavior
                    snapshnum = 0_8
                    if (snap_on_event .and. dkeventsnap == 1_8) then
                        write(iwrite,'(/,a)') '    Snapshots will be saved in file ' // trim(chistoryfname) // ' at every event'
                    elseif (snap_on_event .and. dkeventsnap /=  1_8) then
                        write(iwrite,'(/,a)') '    Snapshots will be saved in file ' // trim(chistoryfname) // ' at every '  // trim(int82str(dkeventsnap)) // ' events'
                    else
                        write(iwrite,'(/,a)') '    Snapshots will be saved in file ' // trim(chistoryfname) // ' every ' // trim(dbl2str(dtsnap)) // ' time units'
                    endif
                endif
            
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1018)
                                
            endif                
                       
            readsnapsh = .true.

        elseif (striccompare(recinput(iw(1):iw(2)),'process_statistics')) then

            if (readprocstat) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (nwords == 1 .or. nwords > 5) then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1019)
            endif
            
            if (nwords == 2 .and. striccompare(recinput(iw(3):iw(4)),'off')) then
                ! output_snapshots is already set to .false. upon initialization
                write(iwrite,'(/,a)') '    Process statistics reporting turned off'
                
            elseif (striccompare(recinput(iw(3):iw(4)),'on')) then
                
                ! At this point output_procstat is set to .false. from initialization
                select case (nwords)
                
                    case (3)
                        ! On-event reporting
                        if (striccompare(recinput(iw(5):iw(6)),'event')) then
                            dkeventprocstat = 1_8
                            output_procstat = .true.
                            procstat_on_event = .true.
                        endif
                        
                    case (4)
                        ! On-time reporting: need to specify timestep
                        if (striccompare(recinput(iw(5):iw(6)),'time') .and. &
                            str_is_real(recinput(iw(7):iw(8)),8)) then
                            output_procstat = .true.
                            procstat_on_event = .false.
                            procstattime = 0.d0
                            dtprocstat = str2dbl(recinput(iw(7):iw(8)))
                        ! On-k^th event reporting: need to specify Deltak
                        elseif (striccompare(recinput(iw(5):iw(6)),'event') .and. &
                            str_is_integer(recinput(iw(7):iw(8)),8)) then
                            dkeventprocstat = str2int(recinput(iw(7):iw(8)))
                            output_procstat = .true.
                            procstat_on_event = .true.
                        endif
                        
                    case (5)
                        ! On-log-time reporting: need to specify initial time and multiplier
                        if (striccompare(recinput(iw(5):iw(6)),'logtime') .and. &
                            str_is_real(recinput(iw(7):iw(8)),8) .and. &
                            str_is_real(recinput(iw(9):iw(10)),8)) then
                            output_procstat = .true.
                            procstat_on_event = .false.
                            procstat_on_logtime = .true.
                            procstattime = str2dbl(recinput(iw(7):iw(8)))
                            dtprocstat = str2dbl(recinput(iw(9):iw(10)))
                        endif

                end select

                if (.not. output_procstat) then
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(1019)
                else
                    open(unit=iprocstat,file=trim(cprocstatfname),status='unknown')
                    q = setvbuf3f(iprocstat,1,0) ! set line-buffered behavior
                    procstatnum = 0_8
                    if (procstat_on_event .and. dkeventprocstat == 1_8) then
                        write(iwrite,'(/,a)') '    Process statistics will be reported in file ' // trim(cprocstatfname) // ' at every event'
                    elseif (procstat_on_event .and. dkeventprocstat /= 1_8) then
                        write(iwrite,'(/,a)') '    Process statistics will be reported in file ' // trim(cprocstatfname) // ' at every '  // trim(int82str(dkeventprocstat)) // ' events'
                    else
                        write(iwrite,'(/,a)') '    Process statistics will be reported in file ' // trim(cprocstatfname) // ' every ' // trim(dbl2str(dtprocstat)) // ' time units'
                    endif
                endif
            
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1019)
                                
            endif                
                       
            readprocstat = .true.

        elseif (striccompare(recinput(iw(1):iw(2)),'species_numbers')) then

            if (readspecnum) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (nwords == 1 .or. nwords > 5) then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1020)
            endif
            
            if (nwords == 2 .and. striccompare(recinput(iw(3):iw(4)),'off')) then
                ! output_snapshots is already set to .false. upon initialization
                write(iwrite,'(/,a)') '    Species numbers reporting turned off'
                
            elseif (striccompare(recinput(iw(3):iw(4)),'on')) then
                
                ! At this point output_procstat is set to .false. from initialization
                select case (nwords)
                
                    case (3)
                        ! On-event reporting
                        if (striccompare(recinput(iw(5):iw(6)),'event')) then
                            dkeventspecnum = 1_8
                            output_specnum = .true.
                            specnum_on_event = .true.
                        endif
                        
                    case (4)
                        ! On-time reporting: need to specify timestep
                        if (striccompare(recinput(iw(5):iw(6)),'time') .and. &
                            str_is_real(recinput(iw(7):iw(8)),8)) then
                            output_specnum = .true.
                            specnum_on_event = .false.
                            specnumtime = 0.d0
                            dtspecnum = str2dbl(recinput(iw(7):iw(8)))
                        ! On-k^th event reporting: need to specify Deltak
                        elseif (striccompare(recinput(iw(5):iw(6)),'event') .and. &
                            str_is_integer(recinput(iw(7):iw(8)),8)) then
                            dkeventspecnum = str2int(recinput(iw(7):iw(8)))
                            output_specnum = .true.
                            specnum_on_event = .true.
                        ! On-event occurrence reporting: need to specify elementary event
                        elseif (striccompare(recinput(iw(5):iw(6)),'elemevent') .and. &
                            str_is_integer(recinput(iw(7):iw(8)),8)) then
                            dkeventspecnum = 1_8
                            output_specnum = .true.
                            specnum_on_eleventoccur = .true.
                            specnum_on_event = .true.
                            ieleventoccur = str2int(recinput(iw(7):iw(8)))
                        endif

                    case (5)
                        ! On-log-time reporting: need to specify initial time and multiplier
                        if (striccompare(recinput(iw(5):iw(6)),'logtime') .and. &
                            str_is_real(recinput(iw(7):iw(8)),8) .and. &
                            str_is_real(recinput(iw(9):iw(10)),8)) then
                            output_specnum = .true.
                            specnum_on_event = .false.
                            specnum_on_logtime = .true.
                            specnumtime = str2dbl(recinput(iw(7):iw(8)))
                            dtspecnum = str2dbl(recinput(iw(9):iw(10)))
                        endif

                end select

                if (.not. output_specnum) then
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(1020)
                else
                    open(unit=ispecnum,file=trim(cspecnumfname),status='unknown')
                    q = setvbuf3f(ispecnum,1,0) ! set line-buffered behavior
                    specnumnum = 0_8
                    if (specnum_on_event .and. dkeventspecnum == 1_8) then
                        if (specnum_on_eleventoccur) then
                            write(iwrite,'(/,a)') '    Species number will be reported in file ' // trim(cspecnumfname) // ' after every occurrence of elementary event ' // trim(int2str(ieleventoccur))
                        else
                            write(iwrite,'(/,a)') '    Species number will be reported in file ' // trim(cspecnumfname) // ' at every event'
                        endif
                    elseif (specnum_on_event .and. dkeventspecnum /= 1_8) then
                        write(iwrite,'(/,a)') '    Species number will be reported in file ' // trim(cspecnumfname) // ' at every '  // trim(int82str(dkeventspecnum)) // ' events'
                    else
                        if (specnum_on_logtime) then
                            write(iwrite,'(/,a)') '    Species number will be reported in file ' // trim(cspecnumfname) // ' at logarithmically-spaced times, starting at ' // trim(dbl2str(specnumtime)) // ' and multiplying with factor ' // trim(dbl2str(dtspecnum))
                        else
                            write(iwrite,'(/,a)') '    Species number will be reported in file ' // trim(cspecnumfname) // ' every ' // trim(dbl2str(dtspecnum)) // ' time units'
                        endif
                    endif
                endif
            
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1019)
                                
            endif                
                       
            readspecnum = .true.

        elseif (striccompare(recinput(iw(1):iw(2)),'event_report')) then

            if (readeventreport) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif
            
            if (nwords == 2 .and. striccompare(recinput(iw(3):iw(4)),'off')) then
                report_events = .false.
                write(iwrite,'(/,a)') '    Event reporting turned off'
            elseif (nwords == 2 .and. striccompare(recinput(iw(3):iw(4)),'on')) then
                report_events = .true.
                write(iwrite,'(/,a)') '    Event reporting turned on'
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1021)
            endif
            
            readeventreport = .true.

        elseif (striccompare(recinput(iw(1):iw(2)),'wall_time')) then

            if (readwalltime) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif
            
            if (nwords /= 2) then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1022)
            else
                if (str_is_integer(recinput(iw(3):iw(4)),4)) then
                    walltime = str2int(recinput(iw(3):iw(4)))
                    readwalltime = .true.
                    write(iwrite,'(/,a)') '    Allowed walltime in seconds: ' // trim(int2str(walltime))
                else
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(1022)
                endif
            endif

        elseif (striccompare(recinput(iw(1):iw(2)),'no_restart')) then

            if (no_restart) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (nwords > 1)    then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1023)
            endif
            
            no_restart = .true.
            write(iwrite,'(/,a)') '    Keyword no_restart parsed. You will not be able to resume the simulation at a later time.'

        elseif (striccompare(recinput(iw(1):iw(2)),'debug_report_processes')) then

            if (debug_report_processes) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (nwords > 1)    then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1024)
            endif
            
            debug_report_processes = .true.
            write(iwrite,'(/,a)') '    [Debug] Process reporting enabled. Data will be written to file ' // trim(cprocdbgfname)
            open(unit=iprocdbg,file=trim(cprocdbgfname),status='unknown')
            q = setvbuf3f(iprocdbg,1,0) ! set line-buffered behavior
            
        elseif (striccompare(recinput(iw(1):iw(2)),'debug_report_global_energetics')) then

            if (debug_report_globenerg) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (nwords > 1)    then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1028)
            endif
            
            debug_report_globenerg = .true.
            write(iwrite,'(/,a)') '    [Debug] Global energetics reporting enabled. Data will be written to file ' // trim(cglbenerdbgfname)
            open(unit=iglbenergdbg,file=trim(cglbenerdbgfname),status='unknown')
            q = setvbuf3f(iglbenergdbg,1,0) ! set line-buffered behavior

        elseif (striccompare(recinput(iw(1):iw(2)),'debug_check_processes')) then

            if (debug_check_processes) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (nwords > 1)    then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1025)
            endif
            
            debug_check_processes = .true.
            write(iwrite,'(/,a)') '    [Debug] Process information structures will be checked during the course of the simulation'

        elseif (striccompare(recinput(iw(1):iw(2)),'debug_check_global_energetics')) then

            if (debug_check_globenerg) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (nwords > 1)    then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1025)
            endif
            
            debug_check_globenerg = .true.
            write(iwrite,'(/,a)') '    [Debug] Energetic cluster contribution information will be checked during the course of the simulation'

        elseif (striccompare(recinput(iw(1):iw(2)),'debug_check_lattice')) then

            if (debug_check_lattice) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (nwords > 1)    then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1026)
            endif
            
            debug_check_lattice = .true.
            write(iwrite,'(/,a)') '    [Debug] Lattice state structures will be checked during the course of the simulation'

        elseif (striccompare(recinput(iw(1):iw(2)),'debug_newtons_method')) then

            if (debug_newtons_method) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
                call error(100)
            endif

            if (nwords > 1)    then
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(1027)
            endif
            
            debug_newtons_method = .true.
            write(iwrite,'(/,a)') '    [Debug] Newton''s method output will be written in file ' // trim(cnewtndbgfname) // &
                                  ' during the course of the simulation'
            open(unit=inewtndbg,file=trim(cnewtndbgfname),status='unknown')
            q = setvbuf3f(inewtndbg,1,0) ! set line-buffered behavior

        elseif (striccompare(recinput(iw(1):iw(2)),'finish')) then

            exit
            
        else

            moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
            call error(1100)

        endif
        
    endif
    
    call getrecord(iread,recinput,irec,io)

enddo

close(iread)

if (io < 0) then ! finish keyword not encountered
    call warning(1100)
endif

! Check for validity of input

if (.not.readiniseq) then ! seed for random number generator not specified
    call system_clock(count=iniseq)
    write(iwrite,'(/,a)') '    No seed was specified for the random number generator; using value ' // trim(int2str(iniseq))
    call initialize_rng(iniseq) ! initialize random number generator
    ! call error(1101)
endif
if (.not.readtemp) then ! temperature not specified
    call error(1103)
endif
if (.not.readpres) then ! pressure not specified
    call error(1105)
endif
if (.not.readngasspecs) then ! number of gas species not specified
    call error(1107)
endif
if (.not.readgasmolweights) then ! gas species molecular weights not specified
    call error(1109)
endif
if (.not.readgasmolfracs) then ! gas species molar fractions not specified
    call error(1111)
endif
if (.not.readnsurfspecs) then ! number of surface species not specified
    call error(1113)
endif
if (.not.readdentation) then ! surface species dentations not specified
    call error(1115)
endif
    
! Automatic assignment of (some) variables that were not read
if (.not.readmaxsteps) then
    maxsteps = huge(maxsteps)
endif
if (.not.readngausspts) then
    ngp = 10
    allocate(gp(ngp))
    allocate(gw(ngp))
    gp = 0
    gw = 0
    if (tpdsim) then
        call gauss_quad_data(ngp,gp,gw) ! initialize Gauss points and weights in interval [0.0,1.0]
        write(iwrite,'(/,a)') '    Default number of points in Gauss quadrature: ' // int2str(ngp)
        write(iwrite,'(/,a,' // int2str(ngp) // '(1x,f10.6))') '    Gauss points:  ', (gp(i),i=1,ngp)
        write(iwrite,'(a,' // int2str(ngp) // '(1x,f10.6))') '    Gauss weights: ', (gw(i),i=1,ngp)
    endif
endif
if (.not.readkmaxnewton) then
    kmaxnewton = 150
    if (tpdsim) then
        write(iwrite,'(/,a)') '    Default number of maximum Newton iterations: ' // int2str(kmaxnewton)
    endif
endif
if (.not.readetol1newton) then
    etol1 = 1.d-9
    if (tpdsim) then
        write(iwrite,'(/,a)') '    Default Dx tolerance for Newton''s method: ' // dbl2str(etol1)
    endif
endif
if (.not.readetol2newton) then
    etol2 = 1.d-9
    if (tpdsim) then
        write(iwrite,'(/,a)') '    Default RHS tolerance for Newton''s method: ' // dbl2str(etol1)
    endif
endif
if (.not.readmaxtime) then
    maxtime = huge(maxtime)
endif
if (.not.readgasspnames) then
    write(iwrite,'(/,a)',advance='no') '    Automatic gas species names: '
    do i = 1,ngasspecs
        gasspecsnames(i) = 'GS' // trim(int2str(i))
        write(iwrite,'(a,1x)',advance='no') trim(gasspecsnames(i))
    enddo
    write(iwrite,'(a)') ''
endif
surfspecsnames(0) = '*'
if (.not.readsurfspnames) then
    write(iwrite,'(/,a)',advance='no') '    Automatic surface species names: '
    do i = 1,nsurfspecs
        surfspecsnames(i) = 'SS' // trim(int2str(i))
        write(iwrite,'(a,1x)',advance='no') trim(surfspecsnames(i))
    enddo
    write(iwrite,'(a)') ''
endif

!! Write part of the history file header:
!if out
!write(ihistory,'(a)',advance='no')   'Gas_Species:        '
!do i = 1,ngasspecs
!    write(ihistory,'(a)',advance='no') gasspecsnames(i)
!enddo
!write(ihistory,'(/,a)',advance='no') 'Surface_Species:    '
!do i = 1,nsurfspecs
!    write(ihistory,'(a)',advance='no') surfspecsnames(i)
!enddo

write(iwrite,'(/,a)') 'Finished reading simulation input.'

return

end subroutine read_simulation_setup

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module simulation_setup_module
