! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

program zacros_main

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Declare modules

! use ifport ! if this is uncommented, remove variable etime from the 
! definition line 36
use error_module

use simulation_setup_module
use lattice_setup_module
use energetics_setup_module
use mechanism_setup_module
use state_setup_module

use lattice_handle_module, only: prealloc_lhm 
use rates_handle_module, only: prealloc_rhm
use energetics_handle_module
use kmc_simulation_handle_module
use sampling_handle_module
use output_handle_module

use restarts_module

!$ use omp_lib, only: omp_get_num_threads

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

implicit none

integer numthreads, mproc, q, setvbuf3f

integer(8) actim1, actim2, cntrate, secs_passed_old, secs_passed_new

real(4) A(2), etime, t1, t2

logical restex, restarted, ressaved

! For the wall-time limited queues we use the following function that helps
! monitor time and stop the simulation when the wall time is approaching
call system_clock(count=actim1,count_rate=cntrate)
secs_passed_old = 0
secs_passed_new = 0

call initialize_warning_counters()

! *** Simulation setup
!$OMP PARALLEL default(shared)
!$
!$OMP SINGLE
numthreads = 1
!$ numthreads = omp_get_num_threads()
!$OMP END SINGLE
!$ 
!$OMP END PARALLEL

call load_restart_info(restarted,restex,numthreads) ! Attempt to load restart information

if (.not. restarted) then ! start a fresh simulation
    
    open(unit=iwrite,file=trim(cgenoutfname),status='unknown')
    q = setvbuf3f(iwrite,1,0) ! set line-buffered behavior

    write(iwrite,'(a)') '+---------------------------------------------------+'
    write(iwrite,'(a)') '|  ZACROS 1.02                                      |'
    write(iwrite,'(a)') '|  GRAPH-THEORETICAL KMC SIMULATION CODE            |'
    write(iwrite,'(a)') '|                                                   |'
    write(iwrite,'(a)') '|  Multiscale Computational Catalysis and           |'
    write(iwrite,'(a)') '|  Materials Science Research Group                 |'
    write(iwrite,'(a)') '|                                                   |'
    write(iwrite,'(a)') '|  Michail Stamatakis, Ph.D.                        |'
    write(iwrite,'(a)') '|  Chemical Engineering Department                  |'
    write(iwrite,'(a)') '|  University College London                        |'
    write(iwrite,'(a)') '+---------------------------------------------------+'
    
    ! *** Parse input files
    call read_simulation_setup()
    call read_lattice_setup()
    call read_energetics_setup()
    call read_mechanism_setup()
    call read_state_setup()
    
    
    call prealloc_kmc(numthreads)
    call prealloc_lhm(numthreads)
    call prealloc_rhm(numthreads)
endif

write(iwrite,'(/,a)') 'Threading information:'
write(iwrite,'(a)')   '~~~~~~~~~~~~~~~~~~~~~~'

if( .TRUE. &
!$ .and. .FALSE. &
  ) then
  write(*, *) "    NO THREADS"
  write(iwrite,*)   '    NO THREADS'
else
  write(*,*) "    WITH THREADS", numthreads
  write(iwrite,*) "    WITH THREADS", numthreads
endif 

if(.not. restarted) then

    ! *** Write headers for output files
    call write_history_header()
    call write_procstat_header()
    call write_specnums_header()

    ! *** Initialization procedures
    call initialize_lattice_state()
    call initialize_energetics()
    call catalogue_all_processes()

endif
! *** Entering the main KMC loop

write(iwrite,'(/,a)') 'Commencing simulation:'
write(iwrite,'(a)')   '~~~~~~~~~~~~~~~~~~~~~~'

call Cpu_Time(t1) ! function for calculating elapsed CPU time

! Open the output files for the sensitivity analysis data
open(unit = SAfnum, status='unknown',file=trim(SAfname),form='unformatted',ACCESS='STREAM')
open(unit = Specfnum, status='unknown',file=trim(Specfname),form='unformatted',ACCESS='STREAM')
open(unit = clusteroccwrite, status='unknown',file=trim(clusoccfname),form='unformatted',ACCESS='STREAM')
open(unit = Ewrite, status='unknown',file=trim(Efname),form='unformatted',ACCESS='STREAM')
open(unit = Histwrite, status='unknown',file=trim(Histfname),form='unformatted',ACCESS='STREAM')
open(unit = Propfnum, status='unknown',file=trim(Propfname),form='unformatted',ACCESS='STREAM')
open(unit = PropCountfnum, status='unknown',file=trim(PropCountfname),form='unformatted',ACCESS='STREAM')

do while (curtime < maxtime .and. curstep < maxsteps)
    
    call system_clock(count=actim2)
    secs_passed_new = (actim2-actim1)/cntrate

    if (readwalltime) then
    
        ! Checking if walltime has been exceeded   
        if (secs_passed_new > walltime) exit
        
        if ( secs_passed_new < secs_passed_old) then
            moreinfo = 'Emergency exit for secs_passed_new = '                 &
                       // trim(int82str(secs_passed_new))                      &
                       // char(13) // char(10)                                 &
                       // '                   secs_passed_old = '              &
                       // trim(int82str(secs_passed_old))
            call warning(1)
            exit
        endif
        
    endif
    
    mproc = event_times_labels(1)
	dtPrior = curtime - prevtime
    prevtime = curtime
    curtime = event_times_heap(1)
    if (curtime<prevtime) then
       write(*,*) prevtime, curtime
       write(*,*) 'ERROR: current time is smaller than previous time'
       stop
    endif
    curstep = curstep + 1_8
    
    call save_snaphshot()
    call save_procstatistics()
    call save_specnums(mproc)

    if (output_procstat) then
        call obtain_elemstep_statistics(mproc)
    endif

    call realize_process(mproc)

enddo

! *** Finalizing the run
call Cpu_Time(t2)
write(iwrite,'(/,a)') 'Simulation stopped:'
write(iwrite,'(a)')   '~~~~~~~~~~~~~~~~~~~'
write(iwrite,'(a)')   ' Current KMC time: ' // trim(real2str(real(curtime,4)))
write(iwrite,'(a)')   ' Events occurred:  ' // trim(int82str(curstep))
write(iwrite,'(a)')   ' Event frequency:  '                                    &
                      // trim(real2str(real(curstep/curtime,4)))

write(iwrite,'(/,a)') 'Performance facts:'
write(iwrite,'(a)')   '~~~~~~~~~~~~~~~~~~'
write(iwrite,'(/,a)') 'Elapsed CPU time:         '                             &
                      // trim(real2str(t2-t1)) // ' seconds'
write(iwrite,'(a)')   'Elapsed clock time:       '                             &
                      // trim(real2str(real(secs_passed_new,4))) // ' seconds'
write(iwrite,'(/,a)') 'Clock time per KMC event: '                             &
                      // trim( real2str( real(secs_passed_new,4)               &
                                         / real(curstep-resstep,4) ))          &
                      // ' seconds'
write(iwrite,'(a)')   'Clock time per KMC time:  '                             &
                      // trim(real2str( real(secs_passed_new,4)                &
                                         / real(curtime-restime,4) ))          &
                      // ' seconds/KMCTimeUnits'
write(iwrite,'(/,a)') 'Events per clock hour:    '                             &
                      // trim( int82str( int((curstep-resstep)                 &
                         * 3600.0/real(secs_passed_new,4),8)))
write(iwrite,'(a)')   'KMC Dt per clock hour:    '                             &
                      // trim(real2str( real(curtime-restime,4)*3600.0         &
                          /real(secs_passed_new,4))) // ' KMCTimeUnits'

if (inewtonstats(0) > 0) then
    write(iwrite,'(/,a)') 'Newton''s method statistics:'
    write(iwrite,'(a)')   '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    write(iwrite,'(/,a)') 'Total number of times run: '                        &
                          // trim(int82str(inewtonstats(0)))
    write(iwrite,'(a)') 'Number of times failed:    '                          &
                        // trim(int82str(inewtonstats(1)))
    write(iwrite,'(a)') 'Avg number of iterations:  '                          &
                        // trim(dbl2str(rnewtonstats(0)))
    write(iwrite,'(a)') 'Maximum Dx error:          '                          &
                        // trim(dbl2str(rnewtonstats(1)))
    write(iwrite,'(a)') 'Maximum RHS error:         '                          &
                        // trim(dbl2str(rnewtonstats(2)))
endif

if (.not.no_restart) then
    call save_restart_info(ressaved)
    if (.not. ressaved) then
        call warning(3)
    endif
endif

write(iwrite,'(/,a)') '> Normal termination <'

! The following two commands are used for debugging.
! For proper testing you can run either the first or the second but not both 
! at the same time.
!call dump_propensity_tables()
!call dump_cluster_contribution_tables() 
! note that the clusters reported here may differ in the non-specific sites 
! but they should be equivalent (manual review by looking at the lattice
! structure is needed)

call cleanup_kshm
call cleanup_lhm
call cleanup_lsm
call cleanup_rhm
call cleanup_ehm
call cleanup_esm
call cleanup_msm
call cleanup_ssm
call cleanup_error_module

close(iwrite)
close(ihistory)
close(iprocstat)
close(iprocdbg)
close(SAfnum)
close(Specfnum)

stop

end program zacros_main
