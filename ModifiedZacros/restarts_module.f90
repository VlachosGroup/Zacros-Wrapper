! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module restarts_module

use constants_module

use heap_functions_module
use mt19937_module

use lattice_setup_module
use mechanism_setup_module
use energetics_setup_module
use simulation_setup_module

use kmc_simulation_handle_module
use lattice_handle_module, only: prealloc_lhm
use rates_handle_module, only: prealloc_rhm

use energetics_handle_module
use sampling_handle_module
 
!$ use omp_lib, only: omp_get_num_threads

implicit none      

integer restimes
integer resstep

real(8) restime

contains
      
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine save_restart_info(ressaved)

use state_setup_module, only: nadsorb

implicit none      

integer i, j, k, m

logical ressaved

ressaved = .false.

open(unit=irestart,file=trim(crestartfname),status='unknown')

rewind(irestart)
write(irestart,'(I10)',err=100) restimes

! *** Save the state of the random number generator

write(irestart,'(a)',err=100) 'Random number generator state'

write(irestart,'(I13)',err=100) rngstate(rngn)
write(irestart,'(6I13)',err=100) (rngstate(i), i = 0,rngn-1)

! *** Save simulation setup information (used in simulation_setup_module)

write(irestart,'(a)',err=100) 'Simulation setup information'

write(irestart,'(2I20,3ES32.16E3)',err=100) (inewtonstats(i), i = 0,1), (rnewtonstats(i), i = 0,2)

write(irestart,'(5I10,I20,3I10)',err=100) iniseq, ngp, kmaxnewton, ngasspecs, nsurfspecs, maxsteps, maxdent, walltime, ieleventoccur
write(irestart,'(7I20)',err=100) curstep, snapshnum, procstatnum, specnumnum, dkeventsnap, dkeventprocstat, dkeventspecnum

write(irestart,'(' // int2str(ngp) // 'ES32.16E3)',err=100) (gp(i), i = 1,ngp)
write(irestart,'(' // int2str(ngp) // 'ES32.16E3)',err=100) (gw(i), i = 1,ngp)

write(irestart,'(' // int2str(ngasspecs) // 'I20)',err=100) (gasspecsnums(i), i = 1,ngasspecs)
write(irestart,'(' // int2str(nsurfspecs+1) // 'I10)',err=100) (surfspecsdent(i), i = 0,nsurfspecs)

write(irestart,'(' // int2str(ngasspecs) // &
     'a' // int2str(nnam0)//')',err=100) (gasspecsnames(i), i = 1,ngasspecs)
write(irestart,'(' // int2str(nsurfspecs+1) // &
     'a' // int2str(nnam0)//')',err=100) (surfspecsnames(i), i = 0,nsurfspecs)

write(irestart,'(' // int2str(ngasspecs) // &
     'ES32.16E3)',err=100) (gasmolfracs(i), i = 1,ngasspecs)
write(irestart,'(' // int2str(ngasspecs) // &
     'ES32.16E3)',err=100) (gasmolweights(i), i = 1,ngasspecs)
write(irestart,'(' // int2str(ngasspecs) // &
     'ES32.16E3)',err=100) (gasenergies(i), i = 1,ngasspecs)

write(irestart,'(14ES32.16E3)',err=100) temp, tramp, pres, etol1, etol2, maxtime, curtime, prevtime, snaptime, dtsnap, & 
                                     procstattime, dtprocstat, specnumtime, dtspecnum
                                     
write(irestart,'(20L2)',err=100) output_snapshots, snap_on_event, report_events, output_procstat, procstat_on_event, &
                        output_specnum, specnum_on_event, specnum_on_eleventoccur, readwalltime, no_restart, debug_report_processes, debug_check_lattice, &
                        snap_on_logtime, procstat_on_logtime, specnum_on_logtime, tpdsim, &
                        debug_check_processes, debug_newtons_method, debug_report_globenerg, debug_check_globenerg

! *** Save lattice setup information (used in lattice_setup_module)

write(irestart,'(a)',err=100) 'Lattice setup information'

write(irestart,'(3I10)',err=100) nsites, nsitetypes, maxcoord

write(irestart,'(' // int2str(nsitetypes) // &
     'a' // int2str(nnam0)//')',err=100) (sitetypenames(i), i = 1,nsitetypes)

do i = 1,nsites
    write(irestart,'(I10,2ES32.16E3,'//int2str(siteneighb1(i,0)+2)//'I10)',err=100) &
           i,xsite(i),ysite(i),sitetype(i),siteneighb1(i,0), &
           (siteneighb1(i,j),j=1,siteneighb1(i,0))
enddo

write(irestart,'(2ES32.16E3)',err=100) v1box(1), v1box(2)
write(irestart,'(2ES32.16E3)',err=100) v2box(1), v2box(2)

! *** Save mechanism information (used in mechanism_setup_handle_module)

write(irestart,'(a)',err=100) 'Mechanism setup information'

write(irestart,'(6I10)',err=100) nelemsteps, nreversible, elemstepnmxsites, elemstepnmxcoord, &
                                 elemstepnmxgas, elemstepnmxangles

write(irestart,'(' // int2str(nelemsteps) // &
     'a' // int2str(nnam0)//')',err=100) (elemstepnames(i), i = 1,nelemsteps)

write(irestart,'(' // int2str(nelemsteps) // 'I10)',err=100) (reverselemstep(i), i = 1,nelemsteps)

write(irestart,'(' // int2str(nelemsteps) // 'I10)',err=100) (elemstepnsites(i), i = 1,nelemsteps)

do i = 1,nelemsteps
    write(irestart,'(' // int2str(2*elemstepnsites(i)+2) // 'I10)',err=100) &
                    (elemstepreactnts(i,j), j = 0,elemstepnsites(i)), &
                    (elemstepproducts(i,j), j = 0,elemstepnsites(i))
enddo                    

do i = 1,nelemsteps
    do j = 1,elemstepnsites(i)
        write(irestart,'('//int2str(elemstepneigh(i,j,0)+3)//'I10)',err=100) &
               j,elemstepstype(i,j),elemstepneigh(i,j,0), &
               (elemstepneigh(i,j,k),k=1,elemstepneigh(i,j,0))
    enddo
enddo                    

do i = 1,nelemsteps
    do j = 1,elemstepnsites(i)
        write(irestart,'(5I10)',err=100) &
               j,(elemsteplatticestatein(i,j,k),k=1,2), &
                 (elemsteplatticestatefn(i,j,k),k=1,2)
    enddo
enddo                    

do i = 1,nelemsteps
    do j = 1,elemstepnsites(i)
        write(irestart,'('//int2str(2*maxdent+1)//'I10)',err=100) &
               j,(elemstepadsorbposin(i,j,k),k=1,maxdent), &
                 (elemstepadsorbposfn(i,j,k),k=1,maxdent)
    enddo
enddo

write(irestart,'('//int2str(nsurfspecs+1)//'I10)',err=100) &
        (nspecparticip(j), j = 0,nsurfspecs)

do i = 0,nsurfspecs
    do j = 1,nspecparticip(i)
        write(irestart,'(2I10)',err=100) (specparticip(i,j,k),k=1,2)
    enddo
enddo

do i = 1,nelemsteps
    write(irestart,'('//int2str(elemstepnmxsites+1)//'I10)',err=100) &
            (elemsteplevels(i,j),j=1,elemstepnmxsites)
enddo

do i = 1,nelemsteps
    write(irestart,'('//int2str(2*elemstepnmxgas+1)//'I10)',err=100) &
            (elemstepgases(i,j), j = 0,2*elemstepnmxgas)
enddo

if (elemstepnmxangles > 0) then
    do i = 1,nelemsteps
        write(irestart,'('//int2str(3*elemstepnmxangles+1)//'I10)',err=100) &
                (elemstepnanglessequenc(i,j), j = 0,3*elemstepnmxangles)
        write(irestart,'('//int2str(elemstepnmxangles)//'ES32.16E3)',err=100) &
                (elemstepangles(i,j),j=1,elemstepnmxangles)
    enddo
endif

do i = 1,nelemsteps
    write(irestart,'(2L2,2I10,ES32.16E3)',err=100) elemstepnomirrorimgs(i), elemstepabslorientat(i), &
             elemsteporientationedge(i,1), elemsteporientationedge(i,2), &
             elemsteporientationangles(i)
enddo

do i = 1,nelemsteps
    write(irestart,'(L2,10ES32.16E3)',err=100) preexpisconst(i), (preexp(i,j),j=0,7), acteng(i), omega(i)
enddo

! *** Save energetic cluster information (used in energetics_setup_module)

write(irestart,'(a)',err=100) 'Energetics setup information'

write(irestart,'(5I10)',err=100) nclusters, clusternmxsites, clusternmxcoord, clusternmxangles, clustermaxlevel

write(irestart,'(' // int2str(nclusters) // &
     'a' // int2str(nnam0)//')',err=100) (clusternames(i), i = 1,nclusters)

write(irestart,'(' // int2str(nclusters) // 'I10)',err=100) (clusternsites(i), i = 1,nclusters)

do i = 1,nclusters
    write(irestart,'(' // int2str(clusternsites(i)+1) // 'I10)',err=100) &
                    (clusterspecs(i,j), j = 0,clusternsites(i))
enddo                    

do i = 1,nclusters
    do j = 1,clusternsites(i)
        write(irestart,'('//int2str(clusterneigh(i,j,0)+3)//'I10)',err=100) &
               j,clusterstype(i,j),clusterneigh(i,j,0), &
               (clusterneigh(i,j,k),k=1,clusterneigh(i,j,0))
    enddo
enddo                    

do i = 1,nclusters
    do j = 1,clusternsites(i)
        write(irestart,'(3I10)',err=100) &
               j,(clusterlatticestate(i,j,k),k=1,2)
    enddo
enddo                    

do i = 1,nclusters
    do j = 1,clusternsites(i)
        write(irestart,'('//int2str(maxdent+1)//'I10)',err=100) &
               j,(clusteradsorbpos(i,j,k),k=1,maxdent)
    enddo
enddo

write(irestart,'('//int2str(nsurfspecs+1)//'I10)',err=100) &
        (nspecclusterparticip(j), j = 0,nsurfspecs)

do i = 0,nsurfspecs
    do j = 1,nspecclusterparticip(i)
        write(irestart,'(2I10)',err=100) (specclusterparticip(i,j,k),k=1,2)
    enddo
enddo

do i = 1,nclusters
    write(irestart,'('//int2str(clusternmxsites+1)//'I10)',err=100) &
            (clusterlevels(i,j),j=1,clusternmxsites)
enddo

if (clusternmxangles > 0) then
    do i = 1,nclusters
        write(irestart,'('//int2str(3*clusternmxangles+1)//'I10)',err=100) &
                (clusternanglessequenc(i,j), j = 0,3*clusternmxangles)
        write(irestart,'(' // int2str(clusternmxangles) // 'ES32.16E3)',err=100) &
                (clusterangles(i,j),j=1,clusternmxangles)
    enddo
endif

do i = 1,nclusters
    write(irestart,'(2L2,2I10,ES32.16E3)',err=100) clusternomirrorimgs(i), clusterabslorientat(i), &
             clusterorientationedge(i,1), clusterorientationedge(i,2), &
             clusterorientationangles(i)
enddo


do i = 1,nclusters
    write(irestart,'(I10,ES32.16E3)',err=100) clustergraphmultipl(i), clusterenrg(i)
enddo

! *** Save lattice state information (used in lattice_handle_module)

write(irestart,'(a)',err=100) 'Lattice state'

write(irestart,'(I10)',err=100) nadsorb

do i = 1,nsites
    write(irestart,'(3I10)',err=100) (latticestate(i,j), j = 1,3)
enddo

do i = 1,nadsorb
    write(irestart,'(' // int2str(surfspecsdent(adsorbspecposi(i,0))+1) // 'I10)',err=100) &
                adsorbspecposi(i,0), &
                (adsorbspecposi(i,j), j = 1,surfspecsdent(adsorbspecposi(i,0)))
enddo

! *** Save energetics information (used in energetics_handle_module)

write(irestart,'(a)',err=100) 'KMC energetics'

write(irestart,'(2I10)',err=100) nglobclust, globclustcapacity0
write(irestart,'(2ES32.16E3)',err=100) energconst, globalenergy

do i = 1,nglobclust
    write(irestart,'(' // int2str(clusternmxsites+1) // 'I10)',err=100) &
                (globclustypesites(i,j), j=0,clusternmxsites)
enddo

do i = 1,nadsorb
    write(irestart,'(' // int2str(2*nadsorglobclusparticip(i)+1) // 'I10)',err=100) &
                nadsorglobclusparticip(i),((adsorglobclusparticip(i,j,m), m=1,2), j=1,nadsorglobclusparticip(i))
enddo

! *** Save simulation information (used in kmc_simulation_handle_module)

write(irestart,'(a)',err=100) 'KMC simulation state'

write(irestart,'(2I10)',err=100) nprocesses, heapcapacity0

do i = 1,nprocesses
    write(irestart,'(2I10,ES32.16E3)',err=100) &
                event_times_labels(i), event_times_indexes(i), event_times_heap(i)
enddo

do i = 1,nprocesses
    write(irestart,'(' // int2str(elemstepnmxsites+1) // 'I10)',err=100) &
                (proctypesites(i,j), j=0,elemstepnmxsites)
enddo

do i = 1,nprocesses
    write(irestart,'(2ES32.16E3)',err=100) procpropenst0(i), procdeltaenrg(i)
enddo

do i = 1,nadsorb
    write(irestart,'(' // int2str(2*nadsorprocparticip(i)+1) // 'I10)',err=100) &
                nadsorprocparticip(i),((adsorprocparticip(i,j,m), m=1,2), j=1,nadsorprocparticip(i))
enddo

do i = 0,nelemsteps
    write(irestart,'(' // int2str(nelemsteps) // 'I20)',err=100) elemstep_noccur(i)
enddo

do i = 0,nelemsteps
    write(irestart,'(' // int2str(nelemsteps) // 'ES32.16E3)',err=100) elemstep_avgtime(i)
enddo

write(irestart,'(a)',err=100) 'Finished'

close(irestart)

ressaved = .true.
write(iwrite,'(/,a)') 'Restart information successfully written in file ' // trim(crestartfname) // '.'

100 return

end subroutine save_restart_info

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine load_restart_info(restarted,restex,numthreads)

use state_setup_module, only: nadsorb

implicit none

integer i, io, j, itmp, k, m, q, setvbuf3f, numthreads
logical restarted,restex

character(1) comments
character(lengthrecinp) recinput
character(256) buf256

restarted = .false.
restimes = 0

! Check if a restart file exists
restex = .false.
inquire(file=crestartfname,exist=restex)
if (.not.restex) return

open(unit=irestart,file=trim(crestartfname),status='old')
read(irestart,'(I10)') restimes
restimes = restimes + 1

open(unit=iwrite,file=trim(cgenoutfname),status='unknown',position='append')
q = setvbuf3f(iwrite,1,0) ! set line-buffered behavior
write(iwrite,'(/,a)') '--------------------------------------------------------'
write(iwrite,'(/,a)') 'Simulation will resume from previously saved state... '
write(iwrite,'(/,a)',advance='no') 'This is restart No ' // trim(int2str(restimes))

! *** Load the state of the random number generator

read(irestart,'(a)') comments

allocate(rngstate(0:rngn))
read(irestart,'(I13)') rngstate(rngn)
read(irestart,'(6I13)') (rngstate(i), i = 0,rngn-1)

! *** Load simulation information (used in simulation_setup_module)

read(irestart,'(a)') comments

allocate(inewtonstats(0:1))
allocate(rnewtonstats(0:2))
read(irestart,'(2I20,3ES32.16E3)') (inewtonstats(i), i = 0,1), (rnewtonstats(i), i = 0,2)

read(irestart,'(5I10,I20,3I10)') iniseq, ngp, kmaxnewton, ngasspecs, nsurfspecs, maxsteps, maxdent, walltime, ieleventoccur
read(irestart,'(7I20)') curstep, snapshnum, procstatnum, specnumnum, dkeventsnap, dkeventprocstat, dkeventspecnum

allocate(gp(ngp))
allocate(gw(ngp))

read(irestart,'(' // int2str(ngp) // 'ES32.16E3)') (gp(i), i = 1,ngp)
read(irestart,'(' // int2str(ngp) // 'ES32.16E3)') (gw(i), i = 1,ngp)

allocate(gasspecsnames(ngasspecs))
allocate(surfspecsnames(0:nsurfspecs))
allocate(gasspecsnums(ngasspecs))
allocate(surfspecsdent(0:nsurfspecs))
allocate(gasmolfracs(ngasspecs))
allocate(gasenergies(ngasspecs))
allocate(gasmolweights(ngasspecs))

! Initialize variables
do i = 1,ngasspecs
    gasspecsnames(i) = 'No_name'
    gasspecsnums(i) = 0
    gasmolfracs(i) = 0.d0
    gasmolweights(i) = 0.d0
enddo
do i = 0,nsurfspecs
    surfspecsnames(i) = 'No_name'
    surfspecsdent(i) = 0
enddo

read(irestart,'(' // int2str(ngasspecs) // 'I20)') (gasspecsnums(i), i = 1,ngasspecs)
read(irestart,'(' // int2str(nsurfspecs+1) // 'I10)') (surfspecsdent(i), i = 0,nsurfspecs)

read(irestart,'(' // int2str(ngasspecs) // &
     'a' // int2str(nnam0)//')') (gasspecsnames(i), i = 1,ngasspecs)
read(irestart,'(' // int2str(nsurfspecs+1) // &
     'a' // int2str(nnam0)//')') (surfspecsnames(i), i = 0,nsurfspecs)

read(irestart,'(' // int2str(ngasspecs) // &
     'ES32.16E3)') (gasmolfracs(i), i = 1,ngasspecs)
read(irestart,'(' // int2str(ngasspecs) // &
     'ES32.16E3)') (gasmolweights(i), i = 1,ngasspecs)
read(irestart,'(' // int2str(ngasspecs) // &
     'ES32.16E3)') (gasenergies(i), i = 1,ngasspecs)

read(irestart,'(14ES32.16E3)') temp, tramp, pres, etol1, etol2, maxtime, curtime, prevtime, snaptime, dtsnap, & 
                                            procstattime, dtprocstat, specnumtime, dtspecnum

read(irestart,'(20L2)') output_snapshots, snap_on_event, report_events, output_procstat, procstat_on_event, &
                        output_specnum, specnum_on_event, specnum_on_eleventoccur, readwalltime, no_restart, debug_report_processes, debug_check_lattice, &
                        snap_on_logtime, procstat_on_logtime, specnum_on_logtime, tpdsim, &
                        debug_check_processes, debug_newtons_method, debug_report_globenerg, debug_check_globenerg

if (output_snapshots) then
    open(unit=ihistory,file=trim(chistoryfname),status='unknown',position='append')
    q = setvbuf3f(ihistory,1,0) ! set line-buffered behavior
endif

if (output_procstat) then
    open(unit=iprocstat,file=trim(cprocstatfname),status='unknown',position='append')
    q = setvbuf3f(iprocstat,1,0) ! set line-buffered behavior
endif

if (output_specnum) then
    open(unit=ispecnum,file=trim(cspecnumfname),status='unknown',position='append')
    q = setvbuf3f(ispecnum,1,0) ! set line-buffered behavior
endif

if (debug_report_processes) then
    open(unit=iprocdbg,file=trim(cprocdbgfname),status='unknown',position='append')
    q = setvbuf3f(iprocdbg,1,0) ! set line-buffered behavior
endif

if (debug_report_globenerg) then
    open(unit=iglbenergdbg,file=trim(cglbenerdbgfname),status='unknown',position='append')
    q = setvbuf3f(iglbenergdbg,1,0) ! set line-buffered behavior
endif

if (debug_newtons_method) then
    open(unit=inewtndbg,file=trim(cnewtndbgfname),status='unknown',position='append')
    q = setvbuf3f(inewtndbg,1,0) ! set line-buffered behavior
endif

write(iwrite,'(a)') ' from KMC step ' // trim(int82str(curstep)) &
             // ' and t = ' // trim(dbl2str(curtime))

resstep = curstep
restime = curtime

! *** Load lattice information (used in lattice_setup_module)

read(irestart,'(a)') comments

read(irestart,'(3I10)') nsites, nsitetypes, maxcoord

allocate(sitetypenames(nsitetypes))
allocate(xsite(nsites))
allocate(ysite(nsites))
allocate(sitetype(nsites))
allocate(siteneighb1(nsites,0:maxcoord))

! Initialize variables
do i = 1,nsitetypes
    sitetypenames(i) = 'No_name'
enddo
do i = 1,nsites
    xsite(i) = 0.d0
    ysite(i) = 0.d0
    sitetype(i) = 0
    do j = 0,maxcoord
        siteneighb1(i,j) = 0
    enddo
enddo

read(irestart,'(' // int2str(nsitetypes) // &
     'a' // int2str(nnam0)//')') (sitetypenames(i), i = 1,nsitetypes)

do i = 1,nsites
    read(irestart,'(i10,2ES32.16E3,'//int2str(maxcoord+2)//'i10)') &
           itmp,xsite(i),ysite(i),sitetype(i),siteneighb1(i,0), &
           (siteneighb1(i,j),j=1,siteneighb1(i,0))
enddo

read(irestart,'(2ES32.16E3)') v1box(1), v1box(2)
read(irestart,'(2ES32.16E3)') v2box(1), v2box(2)

! *** Load mechanism information (used in mechanism_setup_handle_module)

read(irestart,'(a)') comments

read(irestart,'(6I10)') nelemsteps, nreversible, elemstepnmxsites, elemstepnmxcoord, &
                                 elemstepnmxgas, elemstepnmxangles

allocate(elemstepnames(nelemsteps))
allocate(reverselemstep(nelemsteps))
allocate(elemstepnsites(nelemsteps))
allocate(elemstepreactnts(nelemsteps,0:elemstepnmxsites))
allocate(elemstepproducts(nelemsteps,0:elemstepnmxsites))
allocate(elemstepneigh(nelemsteps,elemstepnmxsites,0:elemstepnmxcoord))
allocate(elemstepstype(nelemsteps,elemstepnmxsites))
allocate(elemsteplatticestatein(nelemsteps,elemstepnmxsites,2))
allocate(elemsteplatticestatefn(nelemsteps,elemstepnmxsites,2))
allocate(elemstepadsorbposin(nelemsteps,elemstepnmxsites,maxdent))
allocate(elemstepadsorbposfn(nelemsteps,elemstepnmxsites,maxdent))
allocate(elemstepgases(nelemsteps,0:2*elemstepnmxgas))
allocate(elemstepnanglessequenc(nelemsteps,0:3*elemstepnmxangles))
allocate(elemstepangles(nelemsteps,elemstepnmxangles))
allocate(elemstepnomirrorimgs(nelemsteps))
allocate(elemstepabslorientat(nelemsteps))
allocate(elemsteporientationedge(nelemsteps,2))
allocate(elemsteporientationangles(nelemsteps))
allocate(preexp(nelemsteps,0:7))
allocate(preexpisconst(nelemsteps))
allocate(acteng(nelemsteps))
allocate(omega(nelemsteps))
allocate(nspecparticip(0:nsurfspecs))
allocate(specparticip(0:nsurfspecs,elemstepnmxsites*nelemsteps,2))
allocate(elemsteplevels(nelemsteps,elemstepnmxsites))
allocate(elemstep_noccur(0:nelemsteps))
allocate(elemstep_avgtime(0:nelemsteps))

! Initialize variables
do i = 0,nelemsteps
    elemstep_noccur(i) = 0_8
    elemstep_avgtime(i) = 0.d0
enddo
do i = 1,nelemsteps
    elemstepnames(i) = 'No_name'
    reverselemstep(i) = 0
    elemstepnsites(i) = 0
    elemstepnomirrorimgs(i) = .false.
    elemstepabslorientat(i) = .false.
    elemsteporientationedge(i,1:2) = 0
    elemsteporientationangles(i) = 0.d0
    do j = 1,elemstepnmxsites
        do k = 0,elemstepnmxcoord
            elemstepneigh(i,j,k) = 0
        enddo
    enddo
    elemstepreactnts(i,0) = 0
    elemstepproducts(i,0) = 0
    do j = 1,elemstepnmxsites
        elemstepreactnts(i,j) = -1
        elemstepproducts(i,j) = -1
        elemstepstype(i,j) = -1
        elemsteplevels(i,j) = 0
        do k = 1,2
            elemsteplatticestatein(i,j,k) = 0
            elemsteplatticestatefn(i,j,k) = 0
        enddo
        do k = 1,maxdent
            elemstepadsorbposin(i,j,k) = 0
            elemstepadsorbposfn(i,j,k) = 0
        enddo        
    enddo
    do j = 0,2*elemstepnmxgas
        elemstepgases(i,j) = 0
    enddo
    do j = 0,3*elemstepnmxangles
        elemstepnanglessequenc(i,j) = 0
    enddo
    do j = 1,elemstepnmxangles
        elemstepangles(i,j) = 0.d0
    enddo
    preexp(i,0) = 1.d0
    preexp(i,1:7) = 0.d0
    preexpisconst(i) = .false.
    acteng(i) = 0.d0
    omega(i) = 0.5d0
enddo
do i = 0,nsurfspecs
    nspecparticip(i) = 0
    do j = 1,elemstepnmxsites*nelemsteps
        do k = 1,2
            specparticip(i,j,k) = 0
        enddo
    enddo
enddo

read(irestart,'(' // int2str(nelemsteps) // &
     'a' // int2str(nnam0)//')') (elemstepnames(i), i = 1,nelemsteps)

read(irestart,'(' // int2str(nelemsteps) // 'I10)') (reverselemstep(i), i = 1,nelemsteps)

read(irestart,'(' // int2str(nelemsteps) // 'I10)') (elemstepnsites(i), i = 1,nelemsteps)

do i = 1,nelemsteps
    read(irestart,'(' // int2str(2*elemstepnsites(i)+2) // 'I10)') &
                    (elemstepreactnts(i,j), j = 0,elemstepnsites(i)), &
                    (elemstepproducts(i,j), j = 0,elemstepnsites(i))
enddo                    

do i = 1,nelemsteps
    do j = 1,elemstepnsites(i)
        read(irestart,'('//int2str(elemstepnmxcoord+3)//'I10)') &
               itmp,elemstepstype(i,j),elemstepneigh(i,j,0), &
               (elemstepneigh(i,j,k),k=1,elemstepneigh(i,j,0))
    enddo
enddo                    

do i = 1,nelemsteps
    do j = 1,elemstepnsites(i)
        read(irestart,'(5I10)') &
               itmp,(elemsteplatticestatein(i,j,k),k=1,2), &
                 (elemsteplatticestatefn(i,j,k),k=1,2)
    enddo
enddo                    

do i = 1,nelemsteps
    do j = 1,elemstepnsites(i)
        read(irestart,'('//int2str(2*maxdent+1)//'I10)') &
               itmp,(elemstepadsorbposin(i,j,k),k=1,maxdent), &
                 (elemstepadsorbposfn(i,j,k),k=1,maxdent)
    enddo
enddo

read(irestart,'('//int2str(nsurfspecs+1)//'I10)') &
        (nspecparticip(j), j = 0,nsurfspecs)

do i = 0,nsurfspecs
    do j = 1,nspecparticip(i)
        read(irestart,'(2I10)') (specparticip(i,j,k),k=1,2)
    enddo
enddo

do i = 1,nelemsteps
    read(irestart,'('//int2str(elemstepnmxsites+1)//'I10)') &
            (elemsteplevels(i,j),j=1,elemstepnmxsites)
enddo

do i = 1,nelemsteps
    read(irestart,'('//int2str(2*elemstepnmxgas+1)//'I10)') &
            (elemstepgases(i,j), j = 0,2*elemstepnmxgas)
enddo

if (elemstepnmxangles > 0) then
    do i = 1,nelemsteps
        read(irestart,'('//int2str(3*elemstepnmxangles+1)//'I10)') &
                (elemstepnanglessequenc(i,j), j = 0,3*elemstepnmxangles)
        read(irestart,'(' // int2str(elemstepnmxangles) // 'ES32.16E3)') &
                (elemstepangles(i,j),j=1,elemstepnmxangles)
    enddo
endif

do i = 1,nelemsteps
    read(irestart,'(2L2,2I10,ES32.16E3)') elemstepnomirrorimgs(i), elemstepabslorientat(i), &
             elemsteporientationedge(i,1), elemsteporientationedge(i,2), &
             elemsteporientationangles(i)
enddo

do i = 1,nelemsteps
    read(irestart,'(L2,10ES32.16E3)') preexpisconst(i), (preexp(i,j),j=0,7), acteng(i), omega(i)
enddo

! *** Load energetic cluster information (used in energetics_setup_module)

read(irestart,'(a)') comments

read(irestart,'(5I10)') nclusters, clusternmxsites, clusternmxcoord, clusternmxangles, clustermaxlevel

allocate(clusternames(nclusters))
allocate(clusternsites(nclusters))
allocate(clusterspecs(nclusters,0:clusternmxsites))
allocate(clusterneigh(nclusters,clusternmxsites,0:clusternmxcoord))
allocate(clusterstype(nclusters,clusternmxsites))
allocate(clusterlatticestate(nclusters,clusternmxsites,2))
allocate(clusteradsorbpos(nclusters,clusternmxsites,maxdent))
allocate(clusterenrg(nclusters))
allocate(nspecclusterparticip(0:nsurfspecs))
allocate(specclusterparticip(0:nsurfspecs,clusternmxsites*nclusters,2))
allocate(clustergraphmultipl(nclusters))
allocate(clusterlevels(nclusters,clusternmxsites))
allocate(clusternanglessequenc(nclusters,0:3*clusternmxangles))
allocate(clusterangles(nclusters,clusternmxangles))
allocate(clusternomirrorimgs(nclusters))
allocate(clusterabslorientat(nclusters))
allocate(clusterorientationedge(nclusters,2))
allocate(clusterorientationangles(nclusters))

! Initialize variables
do i = 1,nclusters
    clusternames(i) = 'No_name'
    clusternsites(i) = 0
    clusternanglessequenc(i,0) = 0
    clusternomirrorimgs(i) = .false.
    clusterabslorientat(i) = .false.
    clusterorientationedge(i,1:2) = 0
    clusterorientationangles(i) = 0.d0
    do j = 1,clusternmxangles
        clusterangles(i,j) = 0.d0
        clusternanglessequenc(i,3*j-2:3*j) = 0
    enddo
    do j = 1,clusternmxsites
        do k = 0,clusternmxcoord
            clusterneigh(i,j,k) = 0
        enddo
    enddo
    clusterspecs(i,0) = 0
    do j = 1,clusternmxsites
        clusterspecs(i,j) = -1
        clusterstype(i,j) = -1
        clusterlevels(i,j) = 0
        do k = 1,2
            clusterlatticestate(i,j,k) = 0
        enddo
        do k = 1,maxdent
            clusteradsorbpos(i,j,k) = 0
        enddo        
    enddo
    clusterenrg(i) = 0.d0
    clustergraphmultipl(i) = 1
enddo
do i = 0,nsurfspecs
    nspecclusterparticip(i) = 0
    do j = 1,clusternmxsites*nclusters
        do k = 1,2
            specclusterparticip(i,j,k) = 0
        enddo
    enddo
enddo

read(irestart,'(' // int2str(nclusters) // &
     'a' // int2str(nnam0)//')') (clusternames(i), i = 1,nclusters)

read(irestart,'(' // int2str(nclusters) // 'I10)') (clusternsites(i), i = 1,nclusters)

do i = 1,nclusters
    read(irestart,'(' // int2str(clusternsites(i)+1) // 'I10)') &
                    (clusterspecs(i,j), j = 0,clusternsites(i))
enddo                    

do i = 1,nclusters
    do j = 1,clusternsites(i)
        read(irestart,'('//int2str(clusternmxcoord+3)//'I10)') &
               itmp,clusterstype(i,j),clusterneigh(i,j,0), &
               (clusterneigh(i,j,k),k=1,clusterneigh(i,j,0))
    enddo
enddo                    

do i = 1,nclusters
    do j = 1,clusternsites(i)
        read(irestart,'(3I10)') &
               itmp,(clusterlatticestate(i,j,k),k=1,2)
    enddo
enddo                    

do i = 1,nclusters
    do j = 1,clusternsites(i)
        read(irestart,'('//int2str(maxdent+1)//'I10)') &
               itmp,(clusteradsorbpos(i,j,k),k=1,maxdent)
    enddo
enddo

read(irestart,'('//int2str(nsurfspecs+1)//'I10)') &
        (nspecclusterparticip(j), j = 0,nsurfspecs)

do i = 0,nsurfspecs
    do j = 1,nspecclusterparticip(i)
        read(irestart,'(2I10)') (specclusterparticip(i,j,k),k=1,2)
    enddo
enddo

do i = 1,nclusters
    read(irestart,'('//int2str(clusternmxsites+1)//'I10)') &
            (clusterlevels(i,j),j=1,clusternmxsites)
enddo

if (clusternmxangles > 0) then
    do i = 1,nclusters
        read(irestart,'('//int2str(3*clusternmxangles+1)//'I10)') &
                (clusternanglessequenc(i,j), j = 0,3*clusternmxangles)
        read(irestart,'(' // int2str(clusternmxangles) // 'ES32.16E3)') &
                (clusterangles(i,j),j=1,clusternmxangles)
    enddo
endif

do i = 1,nclusters
    read(irestart,'(2L2,2I10,ES32.16E3)') clusternomirrorimgs(i), clusterabslorientat(i), &
             clusterorientationedge(i,1), clusterorientationedge(i,2), &
             clusterorientationangles(i)
enddo

do i = 1,nclusters
    read(irestart,'(I10,ES32.16E3)') clustergraphmultipl(i), clusterenrg(i)
enddo

! *** Load lattice state information (used in lattice_handle_module)

read(irestart,'(a)') comments

read(irestart,'(I10)') nadsorb

!allocate(latticestate(nsites,3))
call prealloc_kmc(numthreads)
call prealloc_lhm(numthreads)
call prealloc_rhm(numthreads)

latticestate => get_latticestate()

allocate(adsorbspecposi(nsites+elemstepnmxsites,0:maxdent))

do i = 1,nsites
    read(irestart,'(3I10)') (latticestate(i,j), j = 1,3)
enddo

do i = 1,nadsorb
    read(irestart,'(' // int2str(maxdent+1) // 'I10)') &
                adsorbspecposi(i,0), &
                (adsorbspecposi(i,j), j = 1,surfspecsdent(adsorbspecposi(i,0)))
enddo

! *** Load energetics information (used in energetics_handle_module)

read(irestart,'(a)') comments

read(irestart,'(2I10)') nglobclust, globclustcapacity0
read(irestart,'(2ES32.16E3)') energconst, globalenergy

allocate(globclustypesites(globclustcapacity0,0:clusternmxsites))
allocate(nadsorglobclusparticip(nsites))
allocate(adsorglobclusparticip(nsites,50*nclusters*maxcoord,2))

do i = 1,nglobclust
    read(irestart,'(' // int2str(clusternmxsites+1) // 'I10)') &
                (globclustypesites(i,j), j=0,clusternmxsites)
enddo

do i = 1,nadsorb
    read(irestart,'(' // int2str(2*50*nclusters*maxcoord+1) // 'I10)') &
                nadsorglobclusparticip(i),((adsorglobclusparticip(i,j,m), m=1,2), j=1,nadsorglobclusparticip(i))
enddo

! *** Load KMC state simulation information (used in kmc_simulation_handle_module)

read(irestart,'(a)') comments

read(irestart,'(2I10)') nprocesses, heapcapacity0

allocate(event_times_heap(heapcapacity0))
allocate(event_times_labels(heapcapacity0))
allocate(event_times_indexes(heapcapacity0))
allocate(proctypesites(heapcapacity0,0:elemstepnmxsites))
allocate(procpropenst0(heapcapacity0))
allocate(procdeltaenrg(heapcapacity0))
allocate(nadsorprocparticip(nsites))
allocate(adsorprocparticip(nsites,20*maxcoord,2))

recinput = ' '
buf256 = ' '
do i = 1,nprocesses

    read(irestart,'(a' // int2str(len(recinput)) // ')', iostat = io) recinput
    buf256 = recinput(20+1:20+32)
    call remove_leading_blanks(buf256)
    
    if (striccompare(trim(buf256),'inf') .or. striccompare(trim(buf256),'infinity')) then
        read(recinput,'(2I10)') &
                    event_times_labels(i), event_times_indexes(i)
        ! MA -- not portable. Replacing with large number.
        event_times_heap(i) = huge(1.d0)
    else
        read(recinput,'(2I10,ES32.16E3)') &
                    event_times_labels(i), event_times_indexes(i), event_times_heap(i)
    endif
enddo

!do i = 1,nprocesses
!        write(12,'(2I10,ES32.16E3)') &
!                    event_times_labels(i), event_times_indexes(i), event_times_heap(i)   
!enddo

do i = 1,nprocesses
    read(irestart,'(' // int2str(elemstepnmxsites+1) // 'I10)') &
                (proctypesites(i,j), j=0,elemstepnmxsites)
enddo

do i = 1,nprocesses
    read(irestart,'(2ES32.16E3)') procpropenst0(i), procdeltaenrg(i)
enddo

do i = 1,nadsorb
    read(irestart,'(' // int2str(2*10*maxcoord+1) // 'I10)') &
                nadsorprocparticip(i),((adsorprocparticip(i,j,m), m=1,2), j=1,nadsorprocparticip(i))
enddo

do i = 0,nelemsteps
    read(irestart,'(' // int2str(nelemsteps) // 'I20)') elemstep_noccur(i)
enddo

do i = 0,nelemsteps
    read(irestart,'(' // int2str(nelemsteps) // 'ES32.16E3)') elemstep_avgtime(i)
enddo

read(irestart,'(a)') comments

close(irestart)

restarted = .true.
write(iwrite,'(/,a)') 'Restart information successfully read from file ' // trim(crestartfname) // '.'

return

write(iwrite,'(/,a)') 'Resuming simulation:'
write(iwrite,'(a)')   '~~~~~~~~~~~~~~~~~~~~'

110 continue

! The following set of commands is executed only if something goes wrong in the parsing of the restart file

call error(2)

return

end subroutine load_restart_info

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module restarts_module
