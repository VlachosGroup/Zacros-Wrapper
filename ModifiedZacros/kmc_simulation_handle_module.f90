! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module kmc_simulation_handle_module

use constants_module
use error_module
use parser_module
use heap_functions_module
use random_deviates_module

use simulation_setup_module
use lattice_setup_module
use mechanism_setup_module

use lattice_handle_module
use energetics_handle_module
use rates_handle_module

implicit none

integer nprocesses

integer, allocatable :: event_times_labels(:)
integer, allocatable :: event_times_indexes(:)

! Array proctypesites(i1,i2) gives: for i2 = 0 the elementary step number for
! process i1 and for i2 = 1:elemstepnsites the mapping between the pattern sites
! and the actual lattice sites for this process
integer, allocatable :: proctypesites(:,:)
integer, allocatable :: nadsorprocparticip(:)
integer, allocatable :: adsorprocparticip(:,:,:)

real(8), allocatable :: event_times_heap(:)
real(8), allocatable :: procpropenst0(:)
real(8), allocatable :: procdeltaenrg(:)
real(8), allocatable :: f(:)
integer :: MPNi

! Declares a type which is singularly used to hold work arrays.
! The issue is that some functions are called over and over. We do not want to
! have to allocate and reallocate those arrays.
type shared_workarray_type
  integer, allocatable :: neighlist(:)
  integer, allocatable :: validpatrns(:,:)
  integer, allocatable :: indxs(:)
  integer, allocatable :: adsexclude(:)
  integer, allocatable :: startsites(:)
  integer, allocatable :: procsupdate(:)
  integer, allocatable :: adsorbinrange(:)
  integer, allocatable :: sitesaddto(:,:)
  real(8), allocatable :: randnumbrealize(:)
  real(8), allocatable :: randnumbadd(:)
  logical, allocatable :: adsorbookeep(:)
end type 

!$ type threaded_workarray_type
!$   real(8), allocatable, dimension(:) :: times
!$   real(8), allocatable, dimension(:) :: procpropenst0
!$   real(8), allocatable, dimension(:) :: procdeltaenrg
!$   real(8), allocatable, dimension(:) :: activenrg
!$   real(8), allocatable, dimension(:) :: activenrg0
!$   real(8), allocatable, dimension(:) :: deltaenrg
!$   real(8), allocatable, dimension(:) :: deltaenrg0
!$   integer, allocatable, dimension(:) :: processes
!$   integer :: nb
!$ end type

! A private work array.
! When/if we move to multi-processing, we may need more than one.  The
! unfortunate _kshm is because all modules are imported everywhere. E.g.
! modules are rendered completely useless...
type (shared_workarray_type), private, target:: work_shared
! Work array for procsbookkeep
! It is not clear what size this array should be, so we will have to check it is
! large enough everytime.
logical, private, target, allocatable :: prv_procsbookeep(:)
private :: get_shared_workarray, shared_workarray_type
!$ private :: threaded_workarray_type, get_threaded_workarray

!$ type(threaded_workarray_type), private, target, allocatable :: work_prv(:)
contains
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine prealloc_kmc(numthreads)
    use energetics_setup_module, only: clustermaxlevel
    use lattice_setup_module, only: maxcoord, nsites
    use energetics_setup_module, only: clusternmxsites

    integer, intent(in) :: numthreads

!$  integer i

    if( .not. allocated(work_shared%neighlist)) then
        allocate(work_shared%neighlist(100*maxcoord))
        allocate(work_shared%validpatrns(100*maxcoord, elemstepnmxsites))
        allocate(work_shared%indxs(0:max(clustermaxlevel,elemstepnmxsites)))
        allocate(work_shared%adsexclude(0:elemstepnmxsites))
        allocate(work_shared%startsites(elemstepnmxsites))
        allocate(work_shared%procsupdate(0:100*maxcoord))
        allocate(work_shared%adsorbinrange(0:100*maxcoord))
        allocate(work_shared%adsorbookeep(nsites))
        allocate(work_shared%sitesaddto(elemstepnmxsites, maxdent))
        allocate(work_shared%randnumbrealize(100*maxcoord))
        allocate(work_shared%randnumbadd(100*maxcoord))
!$      allocate(work_prv(numthreads))
!$      do i = 1, numthreads
!$        allocate(work_prv(i)%processes(maxcoord*clusternmxsites*nsites))
!$        allocate(work_prv(i)%times(maxcoord*clusternmxsites*nsites))
!$        allocate(work_prv(i)%procpropenst0(maxcoord*clusternmxsites*nsites))
!$        allocate(work_prv(i)%procdeltaenrg(maxcoord*clusternmxsites*nsites))
!$        allocate(work_prv(i)%activenrg(maxcoord*clusternmxsites*nsites))
!$        allocate(work_prv(i)%activenrg0(maxcoord*clusternmxsites*nsites))
!$        allocate(work_prv(i)%deltaenrg(maxcoord*clusternmxsites*nsites))
!$        allocate(work_prv(i)%deltaenrg0(maxcoord*clusternmxsites*nsites))
!$      enddo
    endif
  end subroutine prealloc_kmc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  function get_shared_workarray()
    ! This function returns a pointer to an allocated work array.  In
    ! multiprocessing, we may need to more than one work array. We will only
    ! need modify this function and add a locking mechanism.

    type(shared_workarray_type), pointer :: get_shared_workarray

    get_shared_workarray => work_shared
  end function get_shared_workarray

  function get_procsbookeep(length)

    integer, intent(in) :: length
    logical, dimension(:), pointer:: get_procsbookeep
    
    integer n
    n = max(length, 1)
    if(.not. allocated(prv_procsbookeep)) then
      allocate(prv_procsbookeep(n))
    else if ( n > size(prv_procsbookeep)) then
      deallocate(prv_procsbookeep)
      allocate(prv_procsbookeep(n))
    endif
    get_procsbookeep => prv_procsbookeep
  end function

!$ function get_threaded_workarray()
!$     ! This function returns a pointer to an allocated work array.
!$     ! 
!$     ! The work array is thread specific.
!$     ! This function should only be called from a parallel section.
!$     use omp_lib, only: omp_get_thread_num
!$     
!$     type(threaded_workarray_type), pointer :: get_threaded_workarray
!$     
!$     get_threaded_workarray => work_prv(omp_get_thread_num()+1)
!$ end function get_threaded_workarray

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine cleanup_kshm
      ! Deletes work arrays
      
!$    integer i
      
      if(allocated(work_shared%neighlist)) then
        deallocate(work_shared%neighlist)
        deallocate(work_shared%validpatrns)
        deallocate(work_shared%indxs)
        deallocate(work_shared%adsexclude)
        deallocate(work_shared%startsites)
        deallocate(work_shared%procsupdate)
        deallocate(work_shared%adsorbinrange)
        deallocate(work_shared%sitesaddto)
        deallocate(work_shared%adsorbookeep)
      endif 

!$    if(allocated(work_prv)) then
!$        do i = 1, size(work_prv)
!$            deallocate(work_prv(i)%processes)
!$            deallocate(work_prv(i)%times)
!$            deallocate(work_prv(i)%procpropenst0)
!$            deallocate(work_prv(i)%procdeltaenrg)
!$            deallocate(work_prv(i)%activenrg)
!$            deallocate(work_prv(i)%activenrg0)
!$            deallocate(work_prv(i)%deltaenrg)
!$            deallocate(work_prv(i)%deltaenrg0)
!$        enddo
!$        deallocate(work_prv)
!$    endif
      
      if(allocated(prv_procsbookeep))    deallocate(prv_procsbookeep)
      if(allocated(event_times_labels))  deallocate(event_times_labels)
      if(allocated(event_times_indexes)) deallocate(event_times_indexes)
      if(allocated(proctypesites))       deallocate(proctypesites)
      if(allocated(nadsorprocparticip))  deallocate(nadsorprocparticip)
      if(allocated(adsorprocparticip))   deallocate(adsorprocparticip)
      if(allocated(event_times_heap))    deallocate(event_times_heap)
      if(allocated(procpropenst0))       deallocate(procpropenst0)
      if(allocated(procdeltaenrg))       deallocate(procdeltaenrg)

  end subroutine cleanup_kshm 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine catalogue_all_processes()

  use state_setup_module, only: nadsorb

implicit none

integer i, j, k, mmolec
integer, dimension(0:1) :: adsexclude

heapcapacity0 = 5*maxcoord*nsites
nprocesses = 0

allocate(f(nSAparams))
allocate(event_times_heap(heapcapacity0))
allocate(event_times_labels(heapcapacity0))
allocate(event_times_indexes(heapcapacity0))

! Array proctypesites(i1,i2) gives: for i2 = 0 the elementary step number for
! process i1 and for i2 = 1:elemstepnsites the mapping between the pattern sites
! and the actual lattice sites for this process
allocate(proctypesites(heapcapacity0,0:elemstepnmxsites))
! Array procpropens(i1) gives the propensity of process i1 at curtime = 0 (used
! for debugging purposes)
allocate(procpropenst0(heapcapacity0))
! Array procdeltaenrg(i1) gives the lattice reaction energy of process i1 (used
! for debugging purposes) Note that this does not include the change in the
! energy of the gas species!
allocate(procdeltaenrg(heapcapacity0))
! nadsorprocparticip(i1) gives the number of processes in which adsorbate/entity
! i1 participates
allocate(nadsorprocparticip(nsites))
! adsorprocparticip(i1,i2,i3) gives: for i2 = 1:nadsorprocparticip(i1) and i3 =
! 1 the process numbers in which adsorbate/entity i1 participates. For i3 = 2
! the index of the local entity which appears in the reaction pattern.
allocate(adsorprocparticip(nsites,20*maxcoord,2))

! Initialize variables
event_times_labels = 0
event_times_indexes = 0 
proctypesites = 0
nadsorprocparticip = 0
adsorprocparticip = 0

if (debug_report_processes) then
    write(iprocdbg,'(a)') 'Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
endif

! Loop over all molecules (and empty site entities) that exist on the lattice
do i = 1,nadsorb
    
    ! We want processes in which adsorbate/entity i participates as the 1st
    ! molecule
    mmolec = 1
    adsexclude(0) = 0
    call add_species_related_processes(i,adsexclude,mmolec)
    
enddo

if (debug_check_processes) then
    call check_processes()
endif

if (nprocesses == 0) then    
    call warning(2)
    event_times_labels(1) = 0
    event_times_heap(1) = huge(curtime)
endif

return

end subroutine catalogue_all_processes

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine realize_process(mproc)

  use energetics_setup_module, only: clustermaxlevel, clusternmxsites, clusterOcc, nclusters, clustergraphmultipl
  use state_setup_module, only: nadsorb
  use random_deviates_module, only: uniform_filt_dev
  use lattice_setup_module, only: siteneighb1
!$ use lattice_handle_module, only: copy_latticestate
!$ use omp_lib, only: omp_get_num_threads, omp_get_thread_num
  implicit none

  integer i, iadsorb, istep, igas, ilattcsite, j, jelemstep, cngas, mmolec
  integer iproc, nstart, nneigh, nglobclustbefore, nglobclustafter
  integer k, kspec, m, mproc, nspat, nadsexclude
  integer procsrem(0:250), energrem(0:250), adsrbrem(0:250), mprocsitemap(250)


  real(8) activenrg, activenrg0, deltaenrg, deltaenrg0, deltalattenerg,        &
          eventpropensity0, eventtime, globenergbefore, globenergafter,        &
          procdeltaenergy

  type(shared_workarray_type), pointer:: workarray
!$ type(threaded_workarray_type), pointer:: work_threaded
!$ integer :: nbthreads, t

  logical, dimension(:), pointer:: procsbookeep

  ! gets a work array from the module.
  workarray => get_shared_workarray()

  if (debug_report_processes) then
      write(iprocdbg,'(a,I20,a)') 'KMC step ', curstep,                        &
                                  ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  endif
  if (debug_report_globenerg) then
      write(iglbenergdbg,'(a,I20,a)') 'KMC step ', curstep,                    &
                                      ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  endif

  ! *** Elementary step identification, preparations and reporting
  if (mproc == 0) then
      return
  endif    

  jelemstep = proctypesites(mproc,0)
  nspat = elemstepnsites(jelemstep)
  do i = 1,nspat ! store the sites mapping
      mprocsitemap(i) = proctypesites(mproc,i)
  enddo

  ! Compute the new value of W for each parameter
  do MPNi = 1,nSAparams  
   f(MPNi) = 0				   
  end do

	! Fill in the f values
	DO i = 1, nprocesses
		f(proctypesites(i,0)) = f(proctypesites(i,0)) + procpropenst0(i)			! add the propensity of each reaction to the appropriate place in the vector
	END DO

  ! Time and KMC step advancing is taken care of in the main program
  if (report_events) then
      write(iwrite,'(a)') 'KMC step ' // trim(int82str(curstep))
      write(iwrite,'(a)') '   Elementary step '                                &
                          // trim(elemstepnames(jelemstep))
      write(iwrite,'(a)') '   occurred at time t = ' // trim(dbl2str(curtime))
      write(iwrite,'(a)',advance='no') '   involving site(s):'
      write(iwrite,'(' // int2str(nspat) // '(I4,1x))')                        &
           (mprocsitemap(i),i = 1,nspat)
      write(iwrite,'(a)') ''
  endif

  globenergbefore = globalenergy
  procdeltaenergy = procdeltaenrg(mproc)

  ! *** Removing reactants and associated processes and energetics

  energrem(0) = 0
  procsrem(0) = 0
  adsrbrem(0) = 0

  ! First find all processes that have to be removed as a result of removing the
  ! adsorbates that participate in process mproc
  ! Loop over every participating molecule of current step
  do j = 1, elemstepreactnts(jelemstep,0)
      
      ! Find the lattice site in which the 1st dentate sits...
      ilattcsite = mprocsitemap(elemstepadsorbposin(jelemstep,j,1))
          
      ! ... then find which adsorbate this is on the lattice.
      iadsorb = latticestate(ilattcsite,1)
      adsrbrem(0) = adsrbrem(0) + 1
      adsrbrem(adsrbrem(0)) = iadsorb
      
      ! Find all the processes in which this molecule participates
      procsloop: do k = 1,nadsorprocparticip(iadsorb)
          
          if(all(procsrem(1:procsrem(0))                                       &
                    /= adsorprocparticip(iadsorb,k,1))) then
            procsrem(0) = procsrem(0) + 1
            procsrem(procsrem(0)) = adsorprocparticip(iadsorb,k,1)
          endif
          
      enddo procsloop

      ! Find all the energetic contribution clusters in which this molecule
      ! participates
      energloop: do k = 1,nadsorglobclusparticip(iadsorb)
          
          if(all(energrem(1:energrem(0))                                       &
                    /= adsorglobclusparticip(iadsorb,k,1))) then
              energrem(0) = energrem(0) + 1
              energrem(energrem(0)) = adsorglobclusparticip(iadsorb,k,1)
          endif
          
      enddo energloop

  enddo             
                  
  continue

  ! Remove the processes identified
  call remove_energetics(energrem)
  call remove_processes(procsrem)
  call remove_non_participating_entities(adsrbrem)

  ! Check lattice structure and processes (for debugging)
  if (debug_check_lattice) then
      call check_lattice_state(1,nadsorb)
  endif
  if (debug_check_processes) then
      call check_processes()
  endif
  if (debug_check_globenerg) then
      call check_energetics()
  endif

  ! *** Adding products and associated processes and energetics

  ! Array adsexclude contains the numbers of the newly added entities
  ! (adsorbates)
  nadsexclude = 0

  ! for every product molecule of that step
  do j = 1,elemstepproducts(jelemstep,0) 
      
      nadsorprocparticip(nadsorb+j) = 0
      nadsorglobclusparticip(nadsorb+j) = 0
      
      nadsexclude = nadsexclude + 1
      workarray%adsexclude(nadsexclude) = nadsorb+j
      
      do k = 1,surfspecsdent(elemstepproducts(jelemstep,j))
      
          ilattcsite = mprocsitemap(elemstepadsorbposfn(jelemstep,j,k))
          
          workarray%sitesaddto(j,k) = ilattcsite

      enddo
  enddo

  call add_entities_to_lattice(elemstepproducts(jelemstep,0), &
                               elemstepproducts(jelemstep,1:),&
                               workarray%sitesaddto)
  continue


  !if (.false.) then

  ! Now that all products have been added, find and add the energetic
  ! contribution clusters in which they can participate
  nglobclustbefore = nglobclust

  ! for every product molecule of that step
  do j = 1,elemstepproducts(jelemstep,0) 

      workarray%adsexclude(0) = 0
      ! Find the lattice site in which the 1st dentate sits...
      ilattcsite = mprocsitemap(elemstepadsorbposfn(jelemstep,j,1))
          
      ! ... then find which adsorbate this is on the lattice.
      iadsorb = latticestate(ilattcsite,1)
      
      ! We first want to add clusters in which entity iadsorb participates as
      ! the first one in the pattern
      call add_species_related_energy_contributions(iadsorb,workarray%adsexclude,1)

      workarray%adsexclude(0) = nadsexclude
      do mmolec = 2,clusternmxsites
          ! We now want to add clusters in which entity iadsorb participates as
          ! the second, third etc. -one in the pattern, noting that we have
          ! registered the excluded adsorbates
          call add_species_related_energy_contributions(iadsorb,workarray%adsexclude,mmolec)
      enddo

  enddo


  nglobclustafter = nglobclust

  globenergafter = globalenergy

  if (debug_check_globenerg) then

      if (dabs(globenergafter - globenergbefore - procdeltaenergy) > 1e-10) then
      
          moreinfo = ' KMC step ' // trim(int82str(curstep))                   &
                     // char(13) // char(10) //                                &
                     ' Elementary step ' // trim(elemstepnames(jelemstep))     &
                     // char(13) // char(10) //                                &
                     ' Occurred at time t = ' // trim(dbl2str(curtime))        &
                     // char(13) // char(10) // &
                     ' Involved site(s): '
                     
          do i = 1,nspat
              moreinfo = trim(moreinfo) // ' ' // trim(int2str(mprocsitemap(i)))
          enddo
          
          moreinfo = trim(moreinfo) // char(13) // char(10) //                 & 
                     ' Global lattice energy before  = '                       &
                     // dbl2str(globenergbefore) // char(13) // char(10) //    &
                     ' Lattice delta energy reaction = '                       &
                     // dbl2str(procdeltaenergy) // char(13) // char(10) //    &
                     ' Sum of the above              = '                       &
                     // dbl2str(globenergbefore + procdeltaenergy)             &
                     // char(13) // char(10) //                                &
                     ' Global lattice energy after   = '                       &
                     // dbl2str(globenergafter) // char(13) // char(10)
          
          call error(5020)
      
      endif
      
  endif

  ! Subsequently find the processes which will need updating due to the
  ! energetic interactions. To do so, we first find all adsorbates in the
  ! "range" (level) of the newly added cluster contributions
  ! (nglobclustbefore+1:nglobclustafter) and then the processes in which they
  ! participate. Since we have not yet added any new processes, these are
  ! existing processes that will need updating

  workarray%adsorbinrange(0) = 0
  workarray%procsupdate(0) = 0

  nstart = 0
  ! for every product molecule of that step
  do j = 1,elemstepproducts(jelemstep,0) 

      kspec = elemstepproducts(jelemstep,j)
      
      do k = 1,surfspecsdent(kspec)

          ! Find the lattice site in which the kth dentate sits...
          ilattcsite = mprocsitemap(elemstepadsorbposfn(jelemstep,j,k))

          nstart = nstart + 1
          workarray%startsites(nstart) = ilattcsite

      enddo

      call find_dlevel_neighbors( siteneighb1, clustermaxlevel, nstart,        &
                                  workarray%startsites, nneigh,                &
                                  workarray%indxs, workarray%neighlist)

  enddo


  workarray%adsorbookeep = .false.
  procsbookeep => get_procsbookeep(nprocesses)
  procsbookeep = .false.

  ! Now find the adsorbates in range...
  do j = workarray%indxs(0)+1,workarray%indxs(clustermaxlevel)

      ilattcsite = workarray%neighlist(j)
      iadsorb = latticestate(ilattcsite,1)
      
      ! if this is a new adsorbate add it to the list
      if (.not.workarray%adsorbookeep(iadsorb)) then
      
          workarray%adsorbookeep(iadsorb) = .true.
      
          workarray%adsorbinrange(0) = workarray%adsorbinrange(0) + 1
          workarray%adsorbinrange(workarray%adsorbinrange(0)) = iadsorb
      
      endif
          
  enddo

  ! ... and the processes in which they participate
  do j = 1,workarray%adsorbinrange(0)

      iadsorb = workarray%adsorbinrange(j)

      updprocloop: do m = 1,nadsorprocparticip(iadsorb)

          iproc = adsorprocparticip(iadsorb,m,1)
          
          ! if this is a new process add it to the list
          if (.not. procsbookeep(iproc)) then
              
              procsbookeep(iproc) = .true.
              
              workarray%procsupdate(0) = workarray%procsupdate(0) + 1
              workarray%procsupdate(workarray%procsupdate(0)) = iproc
              
          endif
                  
      enddo updprocloop

  enddo
  do i=1,workarray%procsupdate(0)
     workarray%randnumbrealize(i) = uniform_filt_dev()
  enddo
  ! mstam: Debug - update all processes
  !workarray%procsupdate(0) = nprocesses
  !do i = 1,nprocesses
  !    workarray%procsupdate(i) = i
  !enddo
!$OMP PARALLEL default(shared),                                                &
!$OMP          private(i, iproc, istep, activenrg, activenrg0, deltaenrg,      &
!$OMP                  deltaenrg0, deltalattenerg, eventpropensity0,           &
!$OMP                  eventtime, work_threaded, moreinfo),                    &
!$OMP          shared(nbthreads)
!$  call copy_latticestate()
!$  work_threaded => get_threaded_workarray()
!$  work_threaded%nb = 0
!$OMP SINGLE
!$  nbthreads = omp_get_num_threads()
!$OMP END SINGLE

    ! Now calculate the new rates for the processes that need to be updated
!$OMP DO 
    do i = 1,workarray%procsupdate(0)
        
        iproc = workarray%procsupdate(i)
        istep = proctypesites(iproc,0)
        
        call calculate_elemstep_rate( istep,                                   &
                                      proctypesites(iproc,1:elemstepnmxsites), &
                                      activenrg, activenrg0, deltaenrg,        &
                                      deltaenrg0, deltalattenerg,              &
                                      eventpropensity0, eventtime,             &
                                      workarray%randnumbrealize(i))

!$      work_threaded%nb = work_threaded%nb + 1
!$      work_threaded%processes(work_threaded%nb) = iproc
!$      work_threaded%procpropenst0(work_threaded%nb) = eventpropensity0
!$      work_threaded%procdeltaenrg(work_threaded%nb) = deltalattenerg
!$      work_threaded%times(work_threaded%nb) = eventtime
!$      work_threaded%activenrg(work_threaded%nb) = activenrg
!$      work_threaded%activenrg0(work_threaded%nb) = activenrg0
!$      work_threaded%deltaenrg(work_threaded%nb) = deltaenrg
!$      work_threaded%deltaenrg0(work_threaded%nb) = deltaenrg0
!$  enddo
!$OMP END DO
!$  
!$OMP END PARALLEL
!$ 
!$  do t = 1, nbthreads
!$    do i = 1, work_prv(t)%nb
!$
!$      iproc = work_prv(t)%processes(i)
!$      istep = proctypesites(iproc,0)
!$      eventpropensity0 = work_prv(t)%procpropenst0(i)
!$      deltalattenerg = work_prv(t)%procdeltaenrg(i)
!$      eventtime = work_prv(t)%times(i)
!$      activenrg = work_prv(t)%activenrg(i)
!$      activenrg0 = work_prv(t)%activenrg0(i)
!$      deltaenrg = work_prv(t)%deltaenrg(i)
!$      deltaenrg0 = work_prv(t)%deltaenrg0(i)

        procpropenst0(iproc) = eventpropensity0
        procdeltaenrg(iproc) = deltalattenerg
        call heap_update_element(event_times_heap,event_times_labels, &
                           event_times_indexes,nprocesses,iproc,eventtime)
    
        if (debug_report_processes) then
            write(iprocdbg,'(a)') 'The rate of process '                       &
                                  // trim(int2str(iproc)) // ' was updated.'
            write(iprocdbg,'(a)') '    Elementary step number: '               &
                                  // trim(int2str(istep))
            write(iprocdbg,'(a)') '    Elementary step description: '          &
                                  // trim(elemstepnames(istep))
            write(iprocdbg,'(a)',advance='no')                                 &
              '    Mapping of lattice to pattern sites: '
            write(iprocdbg,'(' // int2str(elemstepnsites(istep)) //            &
                           '(I4,1x))')                                         &
              (proctypesites(iproc,m),m = 1,elemstepnsites(istep))
            write(iprocdbg,'(a)')                                              &
              '    Its activation energy at the zero-coverage limit is '       &
              // trim(dbl2str(activenrg0))
            write(iprocdbg,'(a)')                                              &
              '    Its activation energy for the given configuration is '      &
              // trim(dbl2str(activenrg))
            write(iprocdbg,'(a)')                                              &
              '    Its energy of reaction at the zero-coverage limit is '      &
              // trim(dbl2str(deltaenrg0))
            write(iprocdbg,'(a)')                                              & 
              '    Its energy of reaction for the given configuration is '     &
              // trim(dbl2str(deltaenrg))
            write(iprocdbg,'(a)') '    Its propensity at T0 is '               &
                                  // trim(dbl2str(eventpropensity0))
            write(iprocdbg,'(a)') '    It will occur at t = '                  &
                                  // trim(dbl2str(eventtime))                  &
                                  // ' after Dt = '                            &
                                  // trim(dbl2str(eventtime-curtime))
        endif
        
!$    enddo
    enddo

  continue

  ! Next, find and add the new processes in which the products can participate
  ! for every product molecule of that step
  do j = 1,elemstepproducts(jelemstep,0) 

      ! Find the lattice site in which the 1st dentate sits...
      ilattcsite = mprocsitemap(elemstepadsorbposfn(jelemstep,j,1))
          
      ! ... then find which adsorbate this is on the lattice.
      iadsorb = latticestate(ilattcsite,1)
      
      ! We first want to add processes in which entity iadsorb participates as
      ! the first one in the elementary step
      workarray%adsexclude(0) = 0
      call add_species_related_processes(iadsorb,workarray%adsexclude,1)
  !    call add_species_related_energy_contributions(iadsorb,workarray%adsexclude,1)

  enddo

  ! We subsequently add processes in which the entities participate as the
  ! second, third etc.  CAUTION: In order to avoid double-counting we will have
  ! to be careful not to include processes that involve newly added entities. To
  ! this end, if the order with which the entity participates in the elementary
  ! step is mmolec > 1 one has to check that the newly added entities do not
  ! appear as the mmolec-1, mmolec-2 etc species

  workarray%adsexclude(0) = nadsexclude
  do mmolec = 2,elemstepnmxsites
      ! for every product molecule of that step
      do j = 1,elemstepproducts(jelemstep,0) 
      
          ! Find the lattice site in which the 1st dentate sits...
          ilattcsite = mprocsitemap(elemstepadsorbposfn(jelemstep,j,1))
                  
          ! ... then find which adsorbate this is on the lattice.
          iadsorb = latticestate(ilattcsite,1)
          call add_species_related_processes(iadsorb,workarray%adsexclude,mmolec)

      enddo
      continue
  enddo

  !else ! mstam: Debug - Reinitialize everything
  !
  !    deallocate(globclustypesites)
  !    deallocate(nadsorglobclusparticip)
  !    deallocate(adsorglobclusparticip)
  !    deallocate(event_times_heap)
  !    deallocate(event_times_labels)
  !    deallocate(event_times_indexes)
  !    deallocate(proctypesites)
  !    deallocate(procpropens0)
  !    deallocate(procdeltaenrg0)
  !    deallocate(nadsorprocparticip)
  !    deallocate(adsorprocparticip)
  !    
  !    call initialize_energetics()
  !    call catalogue_all_processes()
  !
  !    globenergafter = globalenergy
  !
  !    if (debug_check_globenerg) then
  !        if (dabs(globenergafter - globenergbefore - procesdeltaenrg) > 1e-10) then
  !            moreinfo = ' KMC step ' // trim(int82str(curstep)) // char(13) // char(10) // &
  !                       ' Elementary step ' // trim(elemstepnames(jelemstep)) // char(13) // char(10) // &
  !                       ' Occurred at time t = ' // trim(dbl2str(curtime)) // char(13) // char(10) // &
  !                       ' Involved site(s): '
  !                       
  !            do i = 1,nspat
  !                moreinfo = trim(moreinfo) // ' ' // trim(int2str(mprocsitemap(i)))
  !            enddo
  !            
  !            moreinfo = trim(moreinfo) // char(13) // char(10) // & 
  !                       ' Global energy before  = ' // dbl2str(globenergbefore) // char(13) // char(10) // &
  !                       ' Delta energy reaction = ' // dbl2str(procesdeltaenrg) // char(13) // char(10) // &
  !                       ' Sum of the above      = ' // dbl2str(globenergbefore + procesdeltaenrg) // char(13) // char(10) // &
  !                       ' Global energy after   = ' // dbl2str(globenergafter) // char(13) // char(10)
  !            
  !            call error(5020)
  !        
  !        endif
  !    endif
  !
  !endif

  continue

  !deallocate(procsupdate)

  do i = 1,elemstepgases(jelemstep,0)
      igas = elemstepgases(jelemstep,2*i-1)
      cngas = elemstepgases(jelemstep,2*i)
      gasspecsnums(igas) = gasspecsnums(igas) + cngas
  enddo

  ! Report global lattice energy (for debugging)
  if (debug_report_globenerg) then
      write(iglbenergdbg,'(a)') 'Current total lattice energy is ' // dbl2str(globalenergy)
  endif

  ! Check lattice structure and processes (for debugging)
  if (debug_check_lattice) then
      call check_lattice_state(1,nadsorb)
  endif
  if (debug_check_processes) then
      call check_processes()
  endif
  if (debug_check_globenerg) then
      call check_energetics()
  endif

  if (nprocesses == 0) then    
      call warning(2)
      event_times_labels(1) = 0
      event_times_heap(1) = huge(curtime)
  endif

  return

end subroutine realize_process

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine add_species_related_processes(i,adsexclude,mmolec)

implicit none

integer, intent(in):: adsexclude(0:)

integer i, ipatrnsite, ilattcsite, j, jelemstep, k, kspec
integer m, mmolec, nspat, dmax, nstart, nneigh, cntrpatrn

type(shared_workarray_type), pointer:: workarray

real(8) activenrg, activenrg0, deltaenrg, deltaenrg0, deltalattenerg,          &
        eventpropensity0, eventtime


! gets a work array from the module.
workarray => get_shared_workarray()

! Add all valid processes in which adsorbate i participates
kspec = adsorbspecposi(i,0)
   
! Loop over all elementary steps to which species kspec participates
elemsteps: do j = 1,nspecparticip(kspec) 
    
    jelemstep = specparticip(kspec,j,1)
    nspat = elemstepnsites(jelemstep)
    
    ! *** INITIAL CHECKS ***
       
    if (specparticip(kspec,j,2) /= mmolec) cycle
    
    ! We will immediately check to see whether molecule i occupies on the
    ! lattice the site types that the pattern jelemstep requires        
    do k = 1,surfspecsdent(kspec)
        ipatrnsite = elemstepadsorbposin(jelemstep,mmolec,k) 
        ilattcsite = adsorbspecposi(i,k)
        if( elemstepstype(jelemstep,ipatrnsite) > 0 &
            .and. elemstepstype(jelemstep,ipatrnsite) /= sitetype(ilattcsite)) &
            cycle elemsteps
    enddo
    
    ! *** CANDIDATE SITES ***
    
    ! Make a list of candidate sites of the lattice. These sites could be
    ! participating in a process pertaining to elementary step jelemstep. 
    ! dmax = max level of elementary step i starting from sites occupied by
    ! molecule i
    dmax = elemsteplevels(jelemstep,mmolec) 
    
    nstart = surfspecsdent(kspec)
    call find_dlevel_neighbors( siteneighb1,dmax,nstart,                       &
                                adsorbspecposi(i, 1:nstart),nneigh,            &
                                workarray%indxs,workarray%neighlist )
        
    ! write(*,('(a)')) 'Entity ' // trim(int2str(i)) // ' (species ' // trim(surfspecsnames(kspec)) // &
    ! ')  could participate in elementary step ' // trim(int2str(specparticip(kspec,j,1))) // &
    ! ' involving some permutation of ' // trim(int2str(nspat)) // & 
    ! ' site(s) out of: '
    ! write(*,('(' // int2str(nneigh) // '(I4,1x))')) (workarray%neighlist(k),k=1,nneigh)
                        
    ! *** FINDING VALID PATTERNS ***
                      
!    call find_elemstep_match( &
!                jelemstep,mmolec, &
!                nneigh,neighlist, & ! the number of candidate sites and the vector containing them
!                adsexclude, &
!                cntrpatrn,validpatrns)

    call find_pattern_match( &
                jelemstep,mmolec,elemstepnmxcoord,elemstepnmxangles, &
                elemstepnsites,elemstepreactnts,elemstepstype, &
                elemstepneigh,elemsteplatticestatein,elemstepadsorbposin, &
                elemstepnanglessequenc,elemstepangles, &
                elemstepnomirrorimgs, &
                elemstepabslorientat,elemsteporientationedge,elemsteporientationangles, &
                ! the number of candidate sites and the vector containing them
                nneigh,workarray%neighlist, & 
                adsexclude, &
                cntrpatrn,workarray%validpatrns)


            
    ! write(*,*) 'Valid permutations:'
    ! do k = 1,cntrpatrn
    !     write(*,*) (validpatrns(k,m), m = 1,nspat)
    ! enddo
    ! continue
        
    ! *** ADDING IDENTIFIED PROCESSES ***
    
    do k=1,cntrpatrn
       workarray%randnumbadd(k) = uniform_filt_dev()
    enddo
    patternloop: do k = 1,cntrpatrn ! For every pattern...
        
        ! At this point we have identified a valid process (with non-zero
        ! propensity) which we need to add to the list. First, calculate the
        ! propensity and the random time that it will take place
        
        call calculate_elemstep_rate( jelemstep, workarray%validpatrns(k,:),   &
                                      activenrg, activenrg0, deltaenrg,        &
                                      deltaenrg0, deltalattenerg,              &
                                      eventpropensity0, eventtime,             &
                                      workarray%randnumbadd(k) )

        call add_process( jelemstep, nspat, workarray%validpatrns(k,:),        &
                          eventpropensity0, deltalattenerg, eventtime )
        
        if (debug_report_processes) then
            write(iprocdbg,'(a)') 'Process ' // trim(int2str(nprocesses))      &
                                  // ' identified:'
            write(iprocdbg,'(a)') '    Elementary step number: '               &
                                  // trim(int2str(jelemstep))
            write(iprocdbg,'(a)') '    Elementary step description: '          &
                                  // trim(elemstepnames(jelemstep))
            write(iprocdbg,'(a)',advance='no')                                 &
              '    Mapping of lattice to pattern sites: '
            write(iprocdbg,'(' // int2str(nspat) // '(I4,1x))')                &
              (workarray%validpatrns(k,m),m = 1,nspat)
            write(iprocdbg,'(a)')                                              &
              '    Its activation energy at the zero-coverage limit is '       &
              // trim(dbl2str(activenrg0))
            write(iprocdbg,'(a)')                                              &
              '    Its activation energy for the given configuration is '      &
              // trim(dbl2str(activenrg))
            write(iprocdbg,'(a)')                                              &
              '    Its energy of reaction at the zero-coverage limit is '      &
              // trim(dbl2str(deltaenrg0))
            write(iprocdbg,'(a)')                                              &
              '    Its energy of reaction for the given configuration is '     &
              // trim(dbl2str(deltaenrg))
            write(iprocdbg,'(a)') '    Its propensity at T0 is '               &
                                  // trim(dbl2str(eventpropensity0))
            write(iprocdbg,'(a)') '    It will occur at t = '                  &
                                  // trim(dbl2str(eventtime)) // ' after Dt = '&
                                  // trim(dbl2str(eventtime-curtime))
        endif

    enddo patternloop
                
enddo elemsteps ! loop over all elementary steps for species kspec

return

end subroutine add_species_related_processes

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine add_process( jelemstep, nspat, validpattern, eventpropensity0,      &
                        deltalattenerg, eventtime)

implicit none

integer jelemstep, nspat, ipatrnsite, ilattcsite, m, p
integer validpattern(:)

real(8) eventpropensity0, deltalattenerg, eventtime

! Add the new time of the process to the heap
call heap_insert_element(event_times_heap,event_times_labels, &
                   event_times_indexes,nprocesses,eventtime,nprocesses + 1)

! In the array that gives the sites mapping, add the new process in the list 
proctypesites(nprocesses,0) = jelemstep
do m = 1,nspat
    proctypesites(nprocesses,m) = validpattern(m)
enddo
procpropenst0(nprocesses) = eventpropensity0
procdeltaenrg(nprocesses) = deltalattenerg

! In the array that gives the processes in which a molecule/entity participates,
! add the new process in the list of the corresponding entity.
do m = 1,elemstepreactnts(jelemstep,0) ! For all reactant molecules...
    ipatrnsite = elemstepadsorbposin(jelemstep,m,1)
    ! find the lattice site where the first dentate sits...
    ilattcsite = validpattern(ipatrnsite) 
    ! and set p equal to the molecule/entity number
    p = latticestate(ilattcsite,1) 
    nadsorprocparticip(p) = nadsorprocparticip(p) + 1
    ! process in which molecule participates
    adsorprocparticip(p,nadsorprocparticip(p),1) = nprocesses
    ! it participates as the mth molecule
    adsorprocparticip(p,nadsorprocparticip(p),2) = m 
enddo            

return

end subroutine add_process

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine remove_processes(procsrem)

implicit none

integer i, iprocrem, iadsorb, ilattcsite, j, jelemstep, nspat, nprocs, m, k
integer nlastprocess
integer procsrem(0:)

do i = 1,procsrem(0)
    
    nlastprocess = nprocesses

    iprocrem = procsrem(i)
    
    jelemstep = proctypesites(iprocrem,0)
    nspat = elemstepnsites(jelemstep)

    ! Remove the processes from the list
    call heap_remove_element(event_times_heap,event_times_labels, &
                       event_times_indexes,nprocesses,iprocrem)
    
    ! Remove that process from the array which gives the processes in which a
    ! molecule/entity participates
    do m = 1,elemstepreactnts(jelemstep,0) ! For all reactant molecules...

        ! find the lattice site in which the 1st dentate sits
        ilattcsite = proctypesites(iprocrem,elemstepadsorbposin(jelemstep,m,1))
        
        ! then find which adsorbate this is on the lattice
        iadsorb = latticestate(ilattcsite,1)
        nprocs = nadsorprocparticip(iadsorb)
        
        do k = 1,nprocs
            if (adsorprocparticip(iadsorb,k,1) == iprocrem) then
                if (k < nprocs) then
                    adsorprocparticip(iadsorb,k,1)                             &
                       = adsorprocparticip(iadsorb,nprocs,1)
                    adsorprocparticip(iadsorb,k,2)                             &
                       = adsorprocparticip(iadsorb,nprocs,2)
                endif
                nadsorprocparticip(iadsorb) = nadsorprocparticip(iadsorb) - 1
            endif
        enddo

    enddo            
    
    ! In the array that gives the sites mapping, remove that process from the
    ! list 
    proctypesites(iprocrem,0:nspat) = -1
    procpropenst0(iprocrem) = -1.d0
    procdeltaenrg(iprocrem) = -1.d0
    
    if (debug_report_processes) then
        write(iprocdbg,'(a)') 'Process ' // trim(int2str(iprocrem))            &
                              // ' was removed.'
    endif
    
    ! Relabel if necessary
    if (iprocrem < nlastprocess) then
        
        jelemstep = proctypesites(nlastprocess,0)
        nspat = elemstepnsites(jelemstep)

        ! The last process nprocessesprev will take the place of process
        ! with id iprocrem
        do m = 0,nspat
            proctypesites(iprocrem,m) = proctypesites(nlastprocess,m)
        enddo
        procpropenst0(iprocrem) = procpropenst0(nlastprocess)
        procdeltaenrg(iprocrem) = procdeltaenrg(nlastprocess)
        
        do m = 1,elemstepreactnts(jelemstep,0) ! For all reactant molecules...
            
            ! find the lattice site in which the 1st dentate sits
            ilattcsite = proctypesites( nlastprocess,                          &
                                        elemstepadsorbposin(jelemstep,m,1) )
            
            ! then find which adsorbate this is on the lattice
            iadsorb = latticestate(ilattcsite,1)
            nprocs = nadsorprocparticip(iadsorb)
            
            do k = 1,nprocs             
                
                ! Change the number of process nprocessesprev in the adsorbate
                ! participation array adsorprocparticip
                if (adsorprocparticip(iadsorb,k,1) == nlastprocess) then
                    adsorprocparticip(iadsorb,k,1) = iprocrem
                endif
                
            enddo
        
        enddo
        
        ! Relabel the process in the heap
        call heap_relabel( event_times_heap, event_times_labels,               &
                           event_times_indexes, nprocesses, nlastprocess,      &
                           iprocrem)
        
        if (debug_report_processes) then
            write(iprocdbg,'(a)') 'Process ' // trim(int2str(nlastprocess))    & 
                                  // ' was relabeled to '                      &
                                  // trim(int2str(iprocrem)) // '.'
        endif

        ! If the relabelling affects the list of subsequent entities to be
        ! removed, take care of it 
        do j = i+1,procsrem(0)
            if (procsrem(j) == nlastprocess) procsrem(j) = iprocrem
        enddo
        
        ! Decrease the number of the last process stored by one
        nlastprocess = nlastprocess - 1
            
    endif
        
enddo                       

return

end subroutine remove_processes

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine remove_non_participating_entities(adsrbrem)

implicit none 

integer i, iprev, ipost, k, m, nprocs, nclusts, nrelabel
integer adsrbrem(0:)
! set in remove_entities... 
! should exist for the duration of the call
integer, pointer :: oldnewlabels(:,:)

! This pointer is set in remove_entities_from_lattice
do i = 1,adsrbrem(0)
    
    k = adsrbrem(i)
    if (nadsorprocparticip(k) > 0) then
        moreinfo = 'Adsorbate ' // trim(int2str(k)) // ' cannot be removed'    &
                   // ' because it participates in '                           &
                   // trim(int2str(nadsorprocparticip(k))) // ' processes.'
        call error(5002)
    endif

    if (nadsorglobclusparticip(k) > 0) then
        moreinfo = 'Adsorbate ' // trim(int2str(k)) // ' cannot be removed'    &
                   // ' because it participates in '                           &
                   // trim(int2str(nadsorglobclusparticip(k)))                 &
                   // ' energetic clusters.'
        call error(5002)
    endif

enddo

nrelabel = 0
call remove_entities_from_lattice( adsrbrem(0), adsrbrem(1:),                  &
                                   nrelabel, oldnewlabels )

! Take care of the relabelling of the adsorbates in the nadsorprocparticip and
! adsorprocparticip as well as the nadsorglobclusparticip and
! adsorglobclusparticip arrays
do i = 1,nrelabel
    
    iprev = oldnewlabels(i,1)
    ipost = oldnewlabels(i,2)

    nprocs = nadsorprocparticip(iprev)
    nadsorprocparticip(ipost) = nadsorprocparticip(iprev)

    do m = 1,nprocs
        adsorprocparticip(ipost,m,1) = adsorprocparticip(iprev,m,1)
        adsorprocparticip(ipost,m,2) = adsorprocparticip(iprev,m,2)
    enddo

    nclusts = nadsorglobclusparticip(iprev)
    nadsorglobclusparticip(ipost) = nadsorglobclusparticip(iprev)

    do m = 1,nclusts
        adsorglobclusparticip(ipost,m,1) = adsorglobclusparticip(iprev,m,1)
        adsorglobclusparticip(ipost,m,2) = adsorglobclusparticip(iprev,m,2)
    enddo

enddo

return

end subroutine remove_non_participating_entities

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_processes()

  use state_setup_module, only: nadsorb

implicit none

integer i, ispec, ilatsite, iadsorb, j, k, m, p, jelemstep

logical foundproc

! First, check the following:
! for each process, the adsorbates involved must have that process listed in
! their participation list
do i = 1,nprocesses ! For each process
    
    ! find which elementary step (type) it corresponds to
    jelemstep = proctypesites(i,0) 
    
    if (jelemstep <= 0) cycle

    ! for every participating molecule of that step
    do j = 1,elemstepreactnts(jelemstep,0) 
    
        ! find the lattice site in which the 1st dentate sits
        ilatsite = proctypesites(i,elemstepadsorbposin(jelemstep,j,1))
        
        ! then find which adsorbate this is on the lattice
        iadsorb = latticestate(ilatsite,1)
        
        ! Check that the mapping between elementary step sites and lattice sites
        ! of that adsorbate is correct
        do m = 1,surfspecsdent(elemstepreactnts(jelemstep,j))
            k = proctypesites(i,elemstepadsorbposin(jelemstep,j,m)) ! lattice site in which 
                                                ! the mth dentate of adsorbate iadsorb sits
            if (iadsorb /= latticestate(k,1)) then
                
                moreinfo = 'Adsorbate ' // trim(int2str(iadsorb))              &
                           // ' is referenced as occupying different sites '   &
                           // 'in elementary step '                            &
                           // trim(elemstepnames(jelemstep)) // ' (process '   &
                           // trim(int2str(i))                                 &
                           // ') than those it really occupies on the lattice.'
                call error(5003)                
                                
            endif
        enddo
    
        ! Check that this process appears in the process-participation list of
        ! that adsorbate
        foundproc = .false.
        do k = 1,nadsorprocparticip(iadsorb)
            if (adsorprocparticip(iadsorb,k,1) == i .and. &
                adsorprocparticip(iadsorb,k,2) == j) then
                foundproc = .true.
                exit
            endif
        enddo
        
        if (.not. foundproc) then

            moreinfo = 'Process ' // trim(int2str(i)) // ' involves adsorbate '&
                       // trim(int2str(iadsorb))                               &
                       // ' which does not list that process in '              &
                       // 'its participation list.'
            call error(5004)                

        endif
                
    enddo
    
enddo

! Second, check the following:
! for each adsorbate on the lattice, each of the processes listed must contain
! that molecule as a reactant
do iadsorb = 1,nadsorb ! For each adsorbate
    ispec = adsorbspecposi(iadsorb,0)
    
    ! for each process in which it participates
    do j = 1,nadsorprocparticip(iadsorb) 
    
        ! the process...
        i = adsorprocparticip(iadsorb,j,1) 
        ! ... in which iadsorb participates as the pth molecule
        p = adsorprocparticip(iadsorb,j,2) 
    
        ! find which elementary step (type) it corresponds to
        jelemstep = proctypesites(i,0) 
        
        do k = 1,surfspecsdent(ispec)
            ! lattice site
            ilatsite = proctypesites(i,elemstepadsorbposin(jelemstep,p,k)) 
           
            if (adsorbspecposi(iadsorb,k) /= ilatsite) then
            
                moreinfo = 'Process ' // trim(int2str(i))                      &
                           // ' lists adsorbate ' //  trim(int2str(iadsorb))   &
                           // ' at different sites than those it actually '    &
                           // 'occupies.'
                call error(5005)                

            endif                
           
            if (latticestate(ilatsite,1) /= iadsorb) then

                moreinfo = 'Process ' // trim(int2str(i))                      &
                           // ' involves adsorbate ' // trim(int2str(iadsorb)) &
                           // ' as adsorbate '                                 &
                           // trim(int2str(latticestate(ilatsite,1))) // '.'
                call error(5006)                

            endif
           
            if (latticestate(ilatsite,2) /= ispec) then

                moreinfo = 'Process ' // trim(int2str(i))                      &
                           // ' lists adsorbate ' // trim(int2str(iadsorb))    &
                           // ' as species '                                   &
                           // trim(int2str(latticestate(ilatsite,2))) // '.'
                call error(5007)                

            endif
            
            if (latticestate(ilatsite,3) /= k) then

                moreinfo = 'Process ' // trim(int2str(i))                      &
                           // ' improperly lists adsorbate '                   &
                           // trim(int2str(iadsorb)) // '.'
                call error(5008)                

            endif
            
        enddo
        
    enddo
    
enddo

return

end subroutine check_processes

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module kmc_simulation_handle_module
