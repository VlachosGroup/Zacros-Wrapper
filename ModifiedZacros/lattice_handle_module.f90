! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module lattice_handle_module

implicit none

integer, pointer :: latticestate(:,:) ! lattice state variable
integer, allocatable :: adsorbspecposi(:,:) ! adsorbate molecule/entity position

! Declares a type which is singularly used to hold work arrays.
! The issue is that some functions are called over and over. We do not want to
! have to allocate and reallocate those arrays. 
! Currently, only holds arrays previously 
type work_array_type
  ! used in find_pattern_match
  integer, allocatable :: fixedsites(:)
  integer, allocatable :: patrnstate(:,:)

  ! used in find_valid_patterns
  integer, allocatable :: speclevs(:)
  integer, allocatable :: nonspeclevs(:)
  integer, allocatable :: sitemap(:)
  integer, allocatable :: sitestocheck(:)
  integer, allocatable :: entitiesmap(:)
  logical, allocatable :: isfixed(:)
  logical, allocatable :: isspecific(:)
  logical, allocatable :: isused(:)
end type 


! A private work array.
! When/if we move to multi-processing, we may need more than one.
type (work_array_type), private, target, allocatable:: prv_work_array(:)
integer, private, allocatable, target:: prv_work_array_labels(:, :)
integer, allocatable, target :: prv_latticestate(:, :, :)
private :: get_work_array, get_work_array_labels, work_array_type

contains
  
  subroutine prealloc_lhm(numthreads)
      ! Allocates memory once and for all.
      use lattice_setup_module, only: maxcoord, nsites
      use mechanism_setup_module, only: elemstepnsites
      use energetics_setup_module, only: clusternmxsites
      !$ use omp_lib, only: omp_in_parallel
      
      integer, intent(in) :: numthreads
      integer :: npatrnsites, i
      
!$    if(omp_in_parallel()) then
!$      stop "This bit of code should not be called in parallel section."
!$    endif
      if( .not. allocated(prv_work_array)) then
          allocate(prv_work_array(numthreads))
          do i = 1, numthreads
              npatrnsites = maxval(elemstepnsites)
              npatrnsites = max(clusternmxsites,npatrnsites)
              allocate(prv_work_array(i)%fixedsites(npatrnsites))
              allocate(prv_work_array(i)%patrnstate(npatrnsites, 3))
              allocate(prv_work_array(i)%speclevs(0:npatrnsites))
              allocate(prv_work_array(i)%nonspeclevs(0:npatrnsites))
              allocate(prv_work_array(i)%sitemap(npatrnsites))
              allocate(prv_work_array(i)%entitiesmap(npatrnsites))
              allocate(prv_work_array(i)%sitestocheck(0:npatrnsites))
              allocate(prv_work_array(i)%isfixed(npatrnsites))
              allocate(prv_work_array(i)%isspecific(npatrnsites))
              allocate(prv_work_array(i)%isused(maxcoord*100))
          enddo
      endif
      if(.not. allocated(prv_latticestate)) then
          allocate(prv_latticestate(nsites, 3, numthreads))
      endif
  end subroutine prealloc_lhm

  function get_work_array()
      ! This function returns a pointer to an allocated work array.
      ! In multiprocessing, we may need to more than one work array. We will only
      ! need modify this function and add a locking mechanism.
      ! Should be called in parallel section only
      
!$    use omp_lib, only: omp_in_parallel, omp_get_thread_num
      
      type(work_array_type), pointer :: get_work_array
      
!$    if(omp_in_parallel()) then
!$      get_work_array => prv_work_array(omp_get_thread_num()+1)
!$    else
        get_work_array => prv_work_array(1)
!$    endif
  end function get_work_array

  function get_work_array_labels(requested_size)
      ! This function returns a pointer to an allocated work array.  In
      ! multiprocessing, we may need to more than one work array. We will only
      ! need modify this function and add a locking mechanism.
      ! :param int noldnew:
      !   Optional parameter giving the size of oldnewlabels
      !   oldnewlabels is only allocated when this routine is called with the
      !   appropriate argument.
      
!$    use omp_lib, only: omp_in_parallel
      
      integer, intent(in) :: requested_size
      integer, pointer :: get_work_array_labels(:,:)
      
      integer :: n
      integer, allocatable :: temporary(:, :)
      
!$    if(omp_in_parallel()) then
!$      stop "This bit of code should only be called in single section."
!$    endif
      n = max(requested_size, 10)
      if(.not. allocated(prv_work_array_labels)) then
        allocate(prv_work_array_labels(n, 2))
      else if(size(prv_work_array_labels, 1) < n) then
        allocate(temporary(size(prv_work_array_labels, 1), 2))
        temporary = prv_work_array_labels
        deallocate(prv_work_array_labels)
        allocate(prv_work_array_labels(n+10, 2))
        prv_work_array_labels = temporary
        deallocate(temporary)
      endif
      
      get_work_array_labels => prv_work_array_labels
  end function get_work_array_labels

!$ subroutine copy_latticestate()
!$     ! Sets lattice state to thread private
!$
!$     use omp_lib, only: omp_get_thread_num
!$     integer :: threadid 
!$     threadid = omp_get_thread_num() + 1
!$  
!$     if(threadid > 1) then
!$       prv_latticestate(:, :, threadid) = prv_latticestate(:, :, 1)
!$     endif
!$OMP  BARRIER
!$ 
!$ end subroutine copy_latticestate

  function get_latticestate()
      ! Gets a thread private lattice state.
!$    use omp_lib, only: omp_in_parallel, omp_get_thread_num
      integer, pointer, dimension(:, :) :: get_latticestate
!$ 
!$    if(omp_in_parallel()) then
!$        get_latticestate => prv_latticestate(:, :, omp_get_thread_num()+1)
!$    else
          get_latticestate => prv_latticestate(:, :, 1)
!$    endif
  end function

  subroutine cleanup_lhm
      ! Deletes work arrays
      
!$    use omp_lib, only: omp_in_parallel
      
      integer i
      
!$    if(omp_in_parallel()) then
!$      stop "This bit of code should not be called in parallele section."
!$    endif
      if(allocated(prv_work_array)) then
        do i = 1, size(prv_work_array)
          deallocate(prv_work_array(i)%fixedsites)
          deallocate(prv_work_array(i)%patrnstate)
          deallocate(prv_work_array(i)%speclevs)
          deallocate(prv_work_array(i)%nonspeclevs)
          deallocate(prv_work_array(i)%sitemap)
          deallocate(prv_work_array(i)%entitiesmap)
          deallocate(prv_work_array(i)%sitestocheck)
          deallocate(prv_work_array(i)%isfixed)
          deallocate(prv_work_array(i)%isspecific)
          deallocate(prv_work_array(i)%isused)
        enddo
        deallocate(prv_work_array)
      endif 

      if(allocated(prv_work_array_labels)) then
        deallocate(prv_work_array_labels)
      endif
      
      if(allocated(prv_latticestate)) deallocate(prv_latticestate)
      if(allocated(adsorbspecposi)) deallocate(adsorbspecposi)

  end subroutine cleanup_lhm 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine initialize_lattice_state()

use state_setup_module, only: specseedmultneigh, speciesonsites, specseedmult, &
                              nseedmult, specseedonsites, nentronsites,        &
                              nentrmult, nadsorb, specseedmultstypes
use lattice_setup_module, only: nsites, maxcoord, siteneighb1
use simulation_setup_module, only: maxdent, surfspecsdent, surfspecsnames
use parser_module, only: int2str
use error_module, only: moreinfo, error
use random_deviates_module, only: uniform_int_dev
use graph_functions_module, only: find_max_level, find_dlevel_neighbors

implicit none

integer i, ident, isite, j, k, m, p, q, nadsrbrem, nrelabel
integer adsexclude(0:1)
integer, allocatable :: adsrbrem(:)
integer, pointer :: oldnewlabels(:,:)
integer, allocatable :: specsadd(:)
integer, allocatable :: sitesspecsadd(:,:)

integer nobjsites, dmax, cntrpatrn1, cntrpatrn
integer, allocatable :: objsites(:), indxs(:)
integer, allocatable :: patrnstate(:,:)
integer, allocatable :: fixedsites(:)
integer, allocatable :: validpatrns(:,:)
integer, allocatable :: patrnsforsite(:,:)
integer, allocatable :: availablepatrns(:)

logical checkneigh, checkstate, checkstypes

character(255) tmpstr

! Initialize the lattice with empty sites everywhere
call initialize_empty_lattice()

! Start seeding the adsorbates as requested in the initial state specification file
allocate(adsrbrem(maxdent))
allocate(specsadd(1))
allocate(sitesspecsadd(1,maxdent))

! We first deal with the explicit seeding instructions, since they are the most 
! specific/restrictive.

! For every explicit seeding instruction, check if the sites referenced are empty
! and if so, seed the species requested on those sites

do i = 1, nentronsites
    
    nadsrbrem = 0
    specsadd(1) = speciesonsites(i)
    
    do j = 1,surfspecsdent(speciesonsites(i))
        
        isite = specseedonsites(i,j)
        sitesspecsadd(1,j) = isite
        
        if (latticestate(isite,2) /= 0) then
        
            tmpstr = ''
            do k = 1,speciesonsites(i)
                tmpstr = tmpstr // trim(int2str(specseedonsites(i,k))) // ' '
            enddo

            moreinfo = 'Initial state seeding instruction ' // trim(int2str(i)) // &
                       ' for species ' // trim(surfspecsnames(speciesonsites(i))) // ' on site(s) ' // &
                       trim(tmpstr) // 'cannot be executed, since site ' // trim(int2str(isite)) // &
                       ' is not empty.'
            call error(5010)            

        endif
        
        nadsrbrem = nadsrbrem + 1
        adsrbrem(nadsrbrem) = latticestate(isite,1)
        
    enddo

    call remove_entities_from_lattice(nadsrbrem, adsrbrem, nrelabel,           &
                                      oldnewlabels)
    call add_entities_to_lattice(1,specsadd,sitesspecsadd)
        
enddo

! Now we deal with the multiple adsorbates seeding instructions, in which the requested number of 
! molecules is seeded randomly on the surface

adsexclude(0) = 0
do i = 1,nentrmult
    
    ident = surfspecsdent(specseedmult(i))
    allocate(patrnstate(ident,3))
    allocate(fixedsites(ident))
    do k = 1,ident
        patrnstate(k,1) = k
        patrnstate(k,2) = 0
        patrnstate(k,3) = 1
        fixedsites(k) = 0
    enddo

    checkneigh = .true.
    checkstate = .true.
    checkstypes = .true.

    call find_max_level(ident,specseedmultneigh(i,:,:),dmax,1,(/1/))

    allocate(indxs(0:dmax))
    allocate(objsites(1+maxcoord**(ident-1)))
    allocate(validpatrns(nsites*maxcoord**(ident-1),ident))
    ! Array patrnsforsite(i1,i2) gives the number of patterns in which site i1 participates for i2 = 0,
    ! and the indexes of these patterns in array validpatrns. Thus, while searching for patterns one
    ! can eliminate overlapping patterns that use the same site.
    allocate(patrnsforsite(nsites,0:ident*maxcoord**(ident-1)))
    allocate(availablepatrns(0:nsites*maxcoord**(ident-1)))

    ! Initialize the patrnsforsite(0)
    patrnsforsite(1:nsites,0) = 0
    
    cntrpatrn = 0
    do k = 1,nsites
        
        call find_dlevel_neighbors(siteneighb1,dmax,1,(/k/), &
                        nobjsites,indxs,objsites)
        
        fixedsites(1) = 1
        
        call find_valid_patterns( &
                nobjsites,objsites, & ! the number of candidate sites and the vector containing them
                ident,specseedmultneigh(i,:,:), & ! number of sites in the pattern, neighboring structure
                specseedmultstypes(i,:),patrnstate, & ! site types and state in the pattern
                (/0/),(/0.d0/), & ! number of angle specifications and pertinent info
                .false., &        ! .false. -> do not prevent mirror images from being detected (i.e. detect them!)
                .false., &        ! .false. -> do not restrict to absolute orientation
                (/0,0/),0.d0, &
                fixedsites,adsexclude,1, &
                checkneigh,.false.,checkstate,checkstypes, &
                cntrpatrn1,validpatrns(cntrpatrn+1:,:))
        
        do j = 1,cntrpatrn1
            do p = 1,ident
                q = patrnsforsite(validpatrns(cntrpatrn+j,p),0)
                patrnsforsite(validpatrns(cntrpatrn+j,p),q+1) = cntrpatrn + j
                patrnsforsite(validpatrns(cntrpatrn+j,p),0) = q + 1
            enddo                
        enddo
        
        cntrpatrn = cntrpatrn + cntrpatrn1
        
    enddo

    availablepatrns(0) = cntrpatrn
    do j = 1,cntrpatrn
        availablepatrns(j) = 1
    enddo

    do j = 1,nseedmult(i)

        !write(*,*) '________________________________________________'        
        !write(*,*) 'Lattice state'
        !do m = 1,nsites
        !    write(*,*) m, ' - ', latticestate(m,:)
        !enddo

    
        !write(*,*) ''
        !write(*,*) 'Valid patterns'
        !do m = 1,cntrpatrn
        !    write(*,*) m, ' - ', validpatrns(m,:)
        !enddo
        !write(*,*) ''
        
        if (availablepatrns(0) == 0) then
            moreinfo = 'Seeding instruction ' // trim(int2str(i)) // ' for ' // trim(int2str(nseedmult(i))) // &
                       ' molecules of species ' // trim(surfspecsnames(specseedmult(i))) // ' can no longer proceed.'
            if (j > 1) moreinfo = trim(moreinfo) // ' ' // trim(int2str(j-1)) // ' molecule(s) has(have) been seeded already. '
            moreinfo = trim(moreinfo) // ' No more molecules can be seeded.'
            call error(4501)
        endif
        
        ! We now prepare the input for the lattice handling subroutines:
        
        ! We have all the empty site patterns available for seeding and randomly choose one of them
        k = 0
        p = uniform_int_dev(availablepatrns(0))
        do m = 1,cntrpatrn
            k = k + availablepatrns(m)
            if (k == p) exit
        enddo
        
        do k = 1,ident
            do p = 1,patrnsforsite(validpatrns(m,k),0)
                if (availablepatrns(patrnsforsite(validpatrns(m,k),p)) == 1) then
                    availablepatrns(patrnsforsite(validpatrns(m,k),p)) = 0
                    availablepatrns(0) = availablepatrns(0) - 1
                endif
            enddo            
        enddo        
        
        ! Remove empty site entities from lattice:
        nadsrbrem = ident
        do k = 1,nadsrbrem
            adsrbrem(k) = latticestate(validpatrns(m,k),1) ! these are empty sites -> single dentate entities
        enddo
        call remove_entities_from_lattice(nadsrbrem, adsrbrem, nrelabel,       &
                                          oldnewlabels)
        
        ! Add the adsorbate:
        specsadd(1) = specseedmult(i)
        do k = 1,ident
            sitesspecsadd(1,k) = validpatrns(m,k)
        enddo
        call add_entities_to_lattice(1,specsadd,sitesspecsadd)

        continue

    enddo

    deallocate(indxs)
    deallocate(objsites)
    deallocate(validpatrns)
    deallocate(patrnsforsite)
    deallocate(availablepatrns)
    
    deallocate(patrnstate)
    deallocate(fixedsites)

enddo    

deallocate(adsrbrem)
deallocate(specsadd)
deallocate(sitesspecsadd)

call check_lattice_state(1,nadsorb)

return

end subroutine initialize_lattice_state

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine initialize_empty_lattice()

use lattice_setup_module, only: nsites
use simulation_setup_module, only: maxdent, debug_check_lattice
use mechanism_setup_module, only: elemstepnmxsites
use state_setup_module, only: nadsorb

implicit none

integer i

! latticestate(i1,i2) gives for site i1: 
!   the molecule number for i2 = 1
!   the species number for i2 = 2
!   the dentate number for i2 = 3
latticestate => get_latticestate()

! adsorbspecposi(i1,i2) gives for molecule number i1, 
!   the species for i2 = 0
!   the site occupied by dentate number 1 < i2 < surfspecsdent(surfspecsdent(i1,0))
allocate(adsorbspecposi(nsites+elemstepnmxsites,0:maxdent))

nadsorb = nsites
do i = 1,nsites
    latticestate(i,1) = i
    latticestate(i,2) = 0
    latticestate(i,3) = 1
    adsorbspecposi(i,0) = 0
    adsorbspecposi(i,1) = i
enddo

! Check lattice structure and processes (for debugging)
if (debug_check_lattice) then
    call check_lattice_state(1,nadsorb)
endif

end subroutine initialize_empty_lattice

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               nobjsites,objsites, &     ! the number of candidate sites and the vector containing them
!               npatrnsites,patrnneigh, & ! number of sites in the pattern, neighboring structure
!               patrnstypes,patrnstate, & ! site types and state in the pattern
!               patrnnanglessequenc,patrnangles, & ! number of and sites in angle specifications and angle values

subroutine find_valid_patterns( &
                nobjsites,objsites, &     
                npatrnsites,patrnneigh, & 
                patrnstypes,patrnstate, & 
                patrnnanglessequenc,patrnangles, & 
                patrnnomirrors, &
                patrnabslorientat, & ! require absolute orientation or not
                patrnorientationedge, &
                partnorientationangle, &
                fixedsites,adsexclude,mmolec, &
                checkneigh,checkgeom,checkstate,checkstypes, &
                cntrpatrn,validpatrns)
use graph_functions_module, only: create_next_permutation
! Subroutine that performs pattern matching.
!
! This works by finding the permutations of nneigh that match the graph pattern as
! well as the site-type requirements and coverage
! Example:
! 
!   |                           1   2   3   4   5   6   7
!   | objsites(1:nobjsites) = [ 14  12  13  15  16  21  23]
!
! the ith element of sitemap gives the neighlist index of the site that corresponds/matches 
! site i of the pattern, for instance if:
!
!   |    sitemap(1:npatrnsites)  = [ 1 6 7 ]
!
! Then: 
!   - site 1 of the pattern corresponds to site 14 of the lattice,
!   - site 2 of the pattern corresponds to site 21 of the lattice,
!   - site 3 of the pattern corresponds to site 23 of the lattice
!
! Say now that we want to fix sites 1 and 3 of the pattern to 14 and 23
! that is objsites(1) and neighlist(7), respectively. Then we would set:
!
!   |   fixedsites = [ 1 0 7 ]
!
! Furthermore, there are cases where we may not be interested in the exact lattice site match for
! some of the pattern sites (when their state is unknown/unspecified). For these we want one solution
! but do not want to iterate in order to find all valis permutations. Then we would set the isspecific
! array to true for these sites. For instance if in the above pattern we had no fixed sites but site 2 was
! non-specific, we would set: 
!
!   |   fixedsites = [ 1 0 7 ]
!   |   isspecific = [ T F T ]
! 
! The subroutine works by generating permutations of npatrnsites out of the nobjsites
! and testing each unfinished permutation for compliance with the pattern. Based on
! the result of the test the next permutation is generated in an efficient way, without
! having to go through all possible permutations. For testing each unfinished permutation
! the user can test neighboring patterns, coverage (state), and site types by setting
! the input variables checkneigh, checkstate, and checkstypes respectively to .true.
! Finds a valid pattern

integer, intent(in)  :: nobjsites ! Number of sites in pattern?
integer, intent(in)  :: npatrnsites, mmolec
integer, intent(in)  :: adsexclude(0:)
integer, intent(in)  :: objsites(:), patrnneigh(:,0:)
integer, intent(in)  :: patrnstypes(:), patrnstate(:,:)
integer, intent(in)  :: patrnnanglessequenc(0:)
integer, intent(in)  :: fixedsites(:)
real(8), intent(in)  :: patrnangles(:)
logical, intent(in)  :: checkneigh, checkgeom, checkstate, checkstypes
logical, intent(in)  :: patrnnomirrors
integer, intent(out) :: cntrpatrn
integer, intent(out) :: validpatrns(:,:)

logical, intent(in)  :: patrnabslorientat
real(8), intent(in)  :: partnorientationangle
integer, intent(in)  :: patrnorientationedge(:)

integer i, nlevel, nindx, lastindx, nperms

logical checkflag, statusflag, completeflag

! A work array to avoid allocation/deallocation in much called routine.
type(work_array_type), pointer:: workarray

cntrpatrn = 0

! First check if we can 
if (nobjsites < npatrnsites) return

! The work array is allocated the first time this function is called.
! It returns a pointer  to a valid array. 
workarray => get_work_array()

! Populate sitemap and sitestocheck vectors with the fixed sites
! The sitestocheck array contains the sites that will be checked for 
! the isomorphism. This saves time, since every time we add a new
! site in the candidate pattern, we do not need to check all
! previously added sites.
workarray%sitestocheck(0) = 0

workarray%sitemap(:npatrnsites) = fixedsites(:npatrnsites);
do i = 1,npatrnsites
    if (fixedsites(i) /= 0) then
        workarray%sitestocheck(0) = workarray%sitestocheck(0) + 1
        ! sitestocheck numbering refers to the pattern
        workarray%sitestocheck(workarray%sitestocheck(0)) = i 
    endif
enddo

workarray%isspecific(:npatrnsites) = patrnstate(:npatrnsites, 1) .ne. -1

! Check validity for the fixed sites of this pattern (if more or equal to one)
if (workarray%sitestocheck(0) == 1) then ! if single fixed site do not check neighboring or geometry (angles)
    call check_pattern(npatrnsites,patrnneigh,patrnstypes,patrnstate, &
                       patrnnanglessequenc,patrnangles,patrnnomirrors, &
                       patrnabslorientat,partnorientationangle,patrnorientationedge, &
                       workarray%sitemap(:npatrnsites),&
                       objsites,workarray%sitestocheck(0:npatrnsites),adsexclude,mmolec, &
                       .false.,.false.,checkstate,checkstypes, &
                       checkflag)
elseif (workarray%sitestocheck(0) > 1) then
    call check_pattern(npatrnsites,patrnneigh,patrnstypes,patrnstate, &
                       patrnnanglessequenc,patrnangles,patrnnomirrors, &
                       patrnabslorientat,partnorientationangle,patrnorientationedge, &
                       workarray%sitemap(:npatrnsites), &
                       objsites,workarray%sitestocheck(0:npatrnsites),adsexclude,mmolec, &
                       checkneigh,checkgeom,checkstate,checkstypes, &
                       checkflag)
else
    checkflag = .true.
endif

if (.not.checkflag) return

! Enter the loop in which permutations of sites are generated and the pattern is checked
nindx = 0
statusflag = .true.
completeflag = .false.
nperms = 1

call create_next_permutation(nobjsites,npatrnsites, &
                             workarray%isfixed(:npatrnsites), &
                             workarray%isspecific(:npatrnsites), &
                             workarray%isused(:nobjsites),    &
                             workarray%sitemap(:npatrnsites), &
                             nindx,lastindx,nlevel, &
                             workarray%speclevs,    &
                             workarray%nonspeclevs, &
                             checkflag,statusflag,completeflag)

if (lastindx == 0) then ! In this case all sites are fixed
    
    ! If the pattern checks, add it to the valid patterns and return
    if (checkflag) then 
        cntrpatrn = cntrpatrn + 1
        do i = 1,npatrnsites
            validpatrns(cntrpatrn,i) = objsites(workarray%sitemap(i))
        enddo
    endif
    
    return ! (return since there is nothing more to do)

endif

! There are more permutations to be generated
do while (statusflag)
    
    workarray%sitestocheck(0) = 1
    workarray%sitestocheck(1) = nlevel

    call check_pattern(npatrnsites,patrnneigh,patrnstypes,patrnstate, &
                       patrnnanglessequenc,patrnangles,patrnnomirrors, &
                       patrnabslorientat,partnorientationangle,patrnorientationedge, &
                       workarray%sitemap(:npatrnsites), &
                       objsites,workarray%sitestocheck(0:npatrnsites), &
                       adsexclude,mmolec, &
                       checkneigh,checkgeom,checkstate,checkstypes, &
                       checkflag)
    
    if (checkflag .and. completeflag) then
        cntrpatrn = cntrpatrn + 1
        do i = 1,npatrnsites
            validpatrns(cntrpatrn,i) = objsites(workarray%sitemap(i))
        enddo
    endif

    call create_next_permutation(nobjsites,npatrnsites, &
                                 workarray%isfixed(:npatrnsites), &
                                 workarray%isspecific(:npatrnsites), &
                                 workarray%isused(:nobjsites),    &
                                 workarray%sitemap(:npatrnsites), &
                                 nindx,lastindx,nlevel, &
                                 workarray%speclevs,    &
                                 workarray%nonspeclevs, &
                                 checkflag,statusflag,completeflag)
    nperms = nperms + 1

enddo

return

end subroutine find_valid_patterns

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_pattern(npatrnsites,patrnneigh,patrnstypes,patrnstate, &
                patrnnanglessequenc,patrnangles,patrnnomirrors, &
                patrnabslorientat,partnorientationangle,patrnorientationedge, &
                sitemap,objsites,sitestocheck,adsexclude,mmolec, &
                checkneigh,checkgeom,checkstate,checkstypes, &
                checkflag )

use lattice_setup_module, only: sitetype, nsites, siteneighb1, v1box, v2box,   &
                                xsite, ysite
use graph_functions_module, only: measureangle, check_isomorph
use constants_module, only: angletol

integer, intent(in)  :: npatrnsites, mmolec
integer, intent(in)  :: objsites(:), patrnneigh(:,0:)
integer, intent(in)  :: patrnstypes(:), patrnstate(:,:)
integer, intent(in)  :: sitemap(npatrnsites), sitestocheck(0:npatrnsites)
integer, intent(in)  :: adsexclude(0:)
integer, intent(in)  :: patrnnanglessequenc(0:)
real(8), intent(in)  :: patrnangles(:)
logical, intent(in)  :: checkneigh, checkgeom, checkstate, checkstypes
logical, intent(in)  :: patrnnomirrors
logical, intent(out) :: checkflag

logical, intent(in)  :: patrnabslorientat
real(8), intent(in)  :: partnorientationangle
integer, intent(in)  :: patrnorientationedge(:)

integer i, j, ipatrnsite, ilattcsite, ilattcsite1, ilattcsite2, ilattcsite3
integer mirrormult

real(8) latticeangle, desiredangle
real(8) x1, x2, x3, y1, y2, y3, c11, c12, c21, c22, c31, c32

logical firstanglecheck

! A work array to avoid allocation/deallocation in much called routine.
type(work_array_type), pointer:: workarray
integer, pointer :: lattstate(:, :)

! The work array is allocated the first time this function is called.
! It returns a pointer  to a valid array. 
workarray => get_work_array()
! Awkward way of making latticestate less of a global...
lattstate => get_latticestate()

checkflag = .true.

if (checkneigh) then
    call check_isomorph(nsites,siteneighb1,npatrnsites,patrnneigh, &
                        sitemap,objsites,sitestocheck,checkflag)
    if (.not.checkflag) return
endif    

if (checkstypes) then
    do i = 1,sitestocheck(0) ! For each of the sites to be checked
        ipatrnsite = sitestocheck(i)
        ilattcsite = objsites(sitemap(ipatrnsite))

        if (patrnstypes(ipatrnsite) <= 0) cycle

        ! check site types in pattern and lattice
        if (sitetype(ilattcsite) /= patrnstypes(ipatrnsite)) then
            checkflag = .false.
            return
        endif
    enddo
endif

if (checkstate) then
    do i = 1,sitestocheck(0) ! For each of the sites to be checked
        
        ipatrnsite = sitestocheck(i)
        ilattcsite = objsites(sitemap(ipatrnsite))

        if (patrnstate(ipatrnsite,1) == -1) cycle ! do not check the unknown/unspecific site(s)
        
        ! check species number in pattern and lattice
        if (lattstate(ilattcsite,2) /= patrnstate(ipatrnsite,2)) then
            checkflag = .false.
            return
        endif
        
        ! check dentate number in pattern and lattice
        if (lattstate(ilattcsite,3) /= patrnstate(ipatrnsite,3)) then
            checkflag = .false.
            return
        endif
        
    enddo
    
    workarray%entitiesmap(:npatrnsites) = 0 ! vector(1:npatrnsites) assignment 
    ! Check entity mapping: for this we need to consider all pattern sites that 
    ! have been assigned to a lattice site (not just the checksites):

    do ipatrnsite = 1,npatrnsites ! For each of the assigned sites
    
        if (sitemap(ipatrnsite) == 0) cycle
        if (patrnstate(ipatrnsite,1) == -1) cycle ! do not check the unknown/unspecific sites

        ilattcsite = objsites(sitemap(ipatrnsite))
        
        if (workarray%entitiesmap(patrnstate(ipatrnsite,1)) == 0) then
        
            workarray%entitiesmap(patrnstate(ipatrnsite,1)) = lattstate(ilattcsite,1)
            
            if (patrnstate(ipatrnsite,1) < mmolec .and. adsexclude(0) > 0) then ! See if any of the entities is in the excluded list
                do i = 1,adsexclude(0)
                    if (workarray%entitiesmap(patrnstate(ipatrnsite,1)) == adsexclude(i)) then
                        checkflag = .false.
                        return
                    endif
                enddo
            endif
            
        else if ( workarray%entitiesmap(patrnstate(ipatrnsite,1)) &
                  .ne. lattstate(ilattcsite,1) ) then
            checkflag = .false.
            return
        endif
        
    enddo

endif

if (checkgeom) then
    ! Check the geometry of the pattern, namely if the angles specified match those of the candidate pattern
    
    ! The firstanglecheck and mirrormult take care of the fact that a mirrored pattern will
    ! have all the angle values equal to minus those specified. Since we have to count that 
    ! mirrored pattern as well, we check to see if the first angle has the negative of the desired value
    ! and then we set the mirrormult = -1 so that in subsequent checks we verify that the subsequent angles have
    ! the negative of the values desired.
    firstanglecheck = .true.
    mirrormult = 1

    anglesloop: do i = 1,patrnnanglessequenc(0) ! For each of the angle specifications
        
        do j = 1,3
            ipatrnsite = patrnnanglessequenc(3*(i-1)+j)
            if (sitemap(ipatrnsite) == 0) then
                cycle anglesloop
            endif
        enddo

        ilattcsite1 = objsites(sitemap(patrnnanglessequenc(3*(i-1)+1)))
        ilattcsite2 = objsites(sitemap(patrnnanglessequenc(3*(i-1)+2)))
        ilattcsite3 = objsites(sitemap(patrnnanglessequenc(3*(i-1)+3)))
        
        x1 = xsite(ilattcsite1)
        y1 = ysite(ilattcsite1)
        call maptoboxfraccoords(v1box,v2box,x1,y1,c11,c12)

        x2 = xsite(ilattcsite2)
        y2 = ysite(ilattcsite2)
        
        call maptoboxfraccoords(v1box,v2box,x2,y2,c21,c22)
        if (c21-c11 < -0.5d0) then
            x2 = x2 + v1box(1)
            y2 = y2 + v1box(2)
            c21 = c21 + 1
        elseif (c21-c11 > 0.5d0) then
            x2 = x2 - v1box(1)
            y2 = y2 - v1box(2)
            c21 = c21 - 1
        endif
        if (c22-c12 < -0.5d0) then
            x2 = x2 + v2box(1)
            y2 = y2 + v2box(2)
            c22 = c22 + 1
        elseif (c22-c12 > 0.5d0) then
            x2 = x2 - v2box(1)
            y2 = y2 - v2box(2)
            c22 = c22 - 1
        endif
        
        x3 = xsite(ilattcsite3)
        y3 = ysite(ilattcsite3)
        
        call maptoboxfraccoords(v1box,v2box,x3,y3,c31,c32)
        if (c31-c21 < -0.5d0) then
            x3 = x3 + v1box(1)
            y3 = y3 + v1box(2)
            c31 = c31 + 1
        elseif (c31-c21 > 0.5d0) then
            x3 = x3 - v1box(1)
            y3 = y3 - v1box(2)
            c31 = c31 - 1
        endif
        if (c32-c22 < -0.5d0) then
            x3 = x3 + v2box(1)
            y3 = y3 + v2box(2)
            c32 = c32 - 1
        elseif (c32-c22 > 0.5d0) then
            x3 = x3 - v2box(1)
            y3 = y3 - v2box(2)
            c32 = c32 + 1
        endif

        latticeangle = measureangle(x1,x2,x3,y1,y2,y3)
        desiredangle = patrnangles(i)

        if (dabs(dabs(latticeangle)-180.d0) < angletol .and. &
            dabs(dabs(desiredangle)-180.d0) < angletol) then
            cycle anglesloop ! treat 180 and -180 as the same
        endif
        
        if (firstanglecheck) then ! in the following line, if patrnnomirrors is .true. mirror image pattern detection is suppressed
            if (.not.(patrnnomirrors) .and. (dabs(latticeangle + desiredangle) < angletol)) then
            ! if (dabs(latticeangle + desiredangle) < angletol) then
                mirrormult = -1
            endif
            firstanglecheck = .false.
        endif
        
        if (dabs(latticeangle - mirrormult*desiredangle) > angletol) then
            checkflag = .false.
            return
        endif

!        if (dabs(latticeangle - desiredangle) > angletol) then
!            checkflag = .false.
!            return
!        endif
                       
    enddo anglesloop
    
    if (patrnabslorientat) then ! if absolute orientation has been specified
        
        if (sitemap(patrnorientationedge(1)) /= 0 .and. &
            sitemap(patrnorientationedge(2)) /= 0 ) then ! if absolute orientation has been specified
                
            ilattcsite1 = objsites(sitemap(patrnorientationedge(1)))
            ilattcsite2 = objsites(sitemap(patrnorientationedge(2)))
        
            x1 = xsite(ilattcsite1)+1.d0
            y1 = ysite(ilattcsite1)
            call maptoboxfraccoords(v1box,v2box,x1,y1,c11,c12)     
        
            x2 = xsite(ilattcsite1)
            y2 = ysite(ilattcsite1)
            call maptoboxfraccoords(v1box,v2box,x2,y2,c21,c22)

            x3 = xsite(ilattcsite2)
            y3 = ysite(ilattcsite2)
            call maptoboxfraccoords(v1box,v2box,x3,y3,c31,c32)
            if (c31-c21 < -0.5d0) then
                x3 = x3 + v1box(1)
                y3 = y3 + v1box(2)
                c31 = c31 + 1
            elseif (c31-c21 > 0.5d0) then
                x3 = x3 - v1box(1)
                y3 = y3 - v1box(2)
                c31 = c31 - 1
            endif
            if (c32-c22 < -0.5d0) then
                x3 = x3 + v2box(1)
                y3 = y3 + v2box(2)
                c32 = c32 - 1
            elseif (c32-c22 > 0.5d0) then
                x3 = x3 - v2box(1)
                y3 = y3 - v2box(2)
                c32 = c32 + 1
            endif
        
            latticeangle = measureangle(x1,x2,x3,y1,y2,y3)
            desiredangle = partnorientationangle
                    
            if (dabs(latticeangle - desiredangle) > angletol) then ! test for exactly the same angle
                checkflag = .false. ! don't return yet though: need to check if angle = 180 or -180
            endif

            if (dabs(dabs(latticeangle)-180.d0) < angletol .and. &
                dabs(dabs(desiredangle)-180.d0) < angletol) then
                checkflag = .true. ! treat 180 and -180 as the same
            endif
        
        endif  
            
    endif
    
endif

return

end subroutine check_pattern

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine maptoboxfraccoords(v1box,v2box,xst,yst,c1,c2)

implicit none

real(8) v1box(2), v2box(2), xst, yst, c1, c2, determ

determ = v1box(1)*v2box(2) - v2box(1)*v1box(2)

if (determ < tiny(1.d0)) then
    c1 = 0.d0
    c2 = 0.d0
    return
endif

c1 = 1.d0/determ*( v2box(2)*xst - v2box(1)*yst)
c2 = 1.d0/determ*(-v1box(2)*xst + v1box(1)*yst)

return

end subroutine maptoboxfraccoords

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine remove_entities_from_lattice(nadsrbrem, adsrbrem, nrelabel,         &
                                        oldnewlabels)

use simulation_setup_module, only: surfspecsdent
use state_setup_module, only: nadsorb
use error_module, only: moreinfo, error
use parser_module, only: int2str

implicit none 

integer, intent(in)    :: nadsrbrem
integer, intent(out)   :: nrelabel
! Pointer to array with old and new labels is internal.
! This arguments returns a pointer to that array.
integer, pointer, intent(out)   :: oldnewlabels(:,:)
integer, intent(inout) :: adsrbrem(:)


integer i, ilattcsite, j, k, kspec

nrelabel = 0 ! counter for the re-labeled entities
oldnewlabels => get_work_array_labels(1)

do i = 1,nadsrbrem

!    continue
!    write(*,*) '______'
!    do j = 1,nsites
!        write(*,*) j,(latticestate(j,k),k=1,3)
!    enddo
    
    k = adsrbrem(i)
    if (k < 1 .or. k > nadsorb) then
        moreinfo = 'Adsorbate ' // trim(int2str(k)) // ' does not exist and cannot be removed.'
        call error(5009)
    endif
    
    kspec = adsorbspecposi(k,0)
    
    ! Remove the entity from the latticestate and the adsorbspecposi arrays
    do j = 1,surfspecsdent(kspec)
        ilattcsite = adsorbspecposi(k,j)
        latticestate(ilattcsite,1) = -1
        latticestate(ilattcsite,2) = -1
        latticestate(ilattcsite,3) = -1
        adsorbspecposi(k,j) = -1
    enddo
    
    ! Relabel the entities if needed
    if (k < nadsorb) then
    
        nrelabel = nrelabel + 1
        oldnewlabels => get_work_array_labels(nrelabel)
        oldnewlabels(nrelabel,1) = nadsorb ! old label
        oldnewlabels(nrelabel,2) = k       ! new label
    
        adsorbspecposi(k,0) = adsorbspecposi(nadsorb,0)
        
        do j = 1,surfspecsdent(adsorbspecposi(nadsorb,0))
            
            adsorbspecposi(k,j) = adsorbspecposi(nadsorb,j)
            
            ilattcsite = adsorbspecposi(nadsorb,j)
            latticestate(ilattcsite,1) = k
                       
        enddo

        ! If the relabelling affects the list of subsequent entities to be removed, 
        ! take care of it 
        do j = i+1,nadsrbrem
            if (adsrbrem(j) == nadsorb) adsrbrem(j) = k
        enddo

    endif

    ! Decrease the number of entities on the lattice by one
    nadsorb = nadsorb - 1
    
enddo

return

end subroutine remove_entities_from_lattice

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine add_entities_to_lattice(nadsrbadd,specsadd,sitesaddto)

use simulation_setup_module, only: surfspecsdent
use parser_module, only: int2str
use error_module, only: moreinfo, error
use lattice_setup_module, only: nsites
use state_setup_module, only: nadsorb

implicit none 

integer i, ilattcsite, j, k, nadsrbadd
integer specsadd(:)
integer sitesaddto(:,:)

do i = 1,nadsrbadd

    nadsorb = nadsorb + 1
    adsorbspecposi(nadsorb,0) = specsadd(i)

    do j = 1,surfspecsdent(specsadd(i))
        
        ilattcsite = sitesaddto(i,j)
    
        if (ilattcsite < 1 .or. ilattcsite > nsites) then
            moreinfo = 'Site ' // trim(int2str(k)) // ' does not exist.'
            call error(5010)
        else
            if (latticestate(ilattcsite,1) > 0) then
                moreinfo = 'Site ' // trim(int2str(k)) // ' is already occupied.'
                call error(5010)            
            endif
        endif
        
        latticestate(ilattcsite,1) = nadsorb
        latticestate(ilattcsite,2) = specsadd(i)
        latticestate(ilattcsite,3) = j
        
        adsorbspecposi(nadsorb,j) = ilattcsite

    enddo
enddo

return

end subroutine add_entities_to_lattice

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_lattice_state(iadsorb1, iadsorb2, myadsorbspecposi, mylatticestate)

use simulation_setup_module, only: surfspecsdent
use parser_module, only: int2str
use error_module, only: moreinfo, error

implicit none

integer, intent(in) :: iadsorb1, iadsorb2
integer, optional, intent(in) :: myadsorbspecposi(:,0:), mylatticestate(:,:)

integer i, ispec, j, jsite

if(present(mylatticestate))then
   do i = iadsorb1,iadsorb2
      ispec = myadsorbspecposi(i,0)
      !write(*,*) surfspecsdent(ispec)
      do j = 1,surfspecsdent(ispec)
         !write(*,*) "test", i, j
         jsite = myadsorbspecposi(i,j)
         !write(*,*) jsite
         !write(*,*) mylatticestate(jsite,1)
         if (mylatticestate(jsite,1) /= i) then
            moreinfo = 'Molecule numbering for lattice site ' // trim(int2str(jsite)) // &
                 ' and molecule/entity ' // trim(int2str(i)) // ' is inconsistent.'
            call error(5001)
         endif
         if (mylatticestate(jsite,2) /= ispec) then
            moreinfo = 'Species numbering for lattice site ' // trim(int2str(jsite)) // &
                 ' and molecule/entity ' // trim(int2str(i)) // ' is inconsistent.'
            call error(5001)
         endif
         if (mylatticestate(jsite,3) /= j) then
            moreinfo = 'Dentate numbering for lattice site ' // trim(int2str(jsite)) // &
                 ' and molecule/entity ' // trim(int2str(i)) // ' is inconsistent.'
            call error(5001)
         endif
      enddo
   enddo
else
! Check for validity of molecule species and position inverse mapping
   do i = iadsorb1,iadsorb2
      ispec = adsorbspecposi(i,0)
      do j = 1,surfspecsdent(ispec)
         jsite = adsorbspecposi(i,j)
         !write(*,*) latticestate(jsite,1)
         if (latticestate(jsite,1) /= i) then
            moreinfo = 'Molecule numbering for lattice site ' // trim(int2str(jsite)) // &
                 ' and molecule/entity ' // trim(int2str(i)) // ' is inconsistent.'
            call error(5001)
         endif
         if (latticestate(jsite,2) /= ispec) then
            moreinfo = 'Species numbering for lattice site ' // trim(int2str(jsite)) // &
                 ' and molecule/entity ' // trim(int2str(i)) // ' is inconsistent.'
            call error(5001)
         endif
         if (latticestate(jsite,3) /= j) then
            moreinfo = 'Dentate numbering for lattice site ' // trim(int2str(jsite)) // &
                 ' and molecule/entity ' // trim(int2str(i)) // ' is inconsistent.'
            call error(5001)
         endif
      enddo
   enddo
endif


return

end subroutine check_lattice_state

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine find_pattern_match( &
                jpattern,mmolec,patrnnmxcoord,patrnnmxangles, &
                patrnlistnsites,patrnlistspecies,patrnliststype, &
                patrnlistneigh,patrnlistlatticestate,patrnlistadsorbpos, &
                patrnlistnanglessequenc,patrnlistangles, &
                patrnlistnomirrors, &
                patrnlistabslorientat,patrnlistorientationedge,partnlistorientationangles, &
                nneigh,neighlist, & ! the number of candidate sites and the vector containing them
                adsexclude, &
                cntrpatrn,validpatrns)

use simulation_setup_module, only: surfspecsdent
implicit none

integer, intent(in)  :: jpattern, mmolec, patrnnmxcoord, patrnnmxangles, nneigh
integer, intent(in)  :: patrnlistnsites(:)
integer, intent(in)  :: patrnlistspecies(:,0:)
integer, intent(in)  :: patrnliststype(:,:)
integer, intent(in)  :: patrnlistneigh(:,:,0:)
integer, intent(in)  :: patrnlistlatticestate(:,:,:)
integer, intent(in)  :: patrnlistadsorbpos(:,:,:)
integer, intent(in)  :: patrnlistnanglessequenc(:,0:)
integer, intent(in)  :: patrnlistorientationedge(:,:)
integer, intent(in)  :: adsexclude(0:)
integer, intent(in)  :: neighlist(:)
real(8), intent(in)  :: patrnlistangles(:,:)
real(8), intent(in)  :: partnlistorientationangles(:)
logical, intent(in)  :: patrnlistnomirrors(:)
logical, intent(in)  :: patrnlistabslorientat(:)
integer, intent(out) :: cntrpatrn
integer, intent(out) :: validpatrns(:, :)


integer i, ident, j, k
integer kspec, npatrnsites

logical checkneigh,checkgeom,checkstate,checkstypes
type(work_array_type), pointer:: workarray

! This subroutine acts as an interface to the find_valid_patterns subroutine. 
! The latter finds elementary step pattern matches. 

! Explanation of the input:
!  jpattern points to an elementary step or cluster pattern to be searched for.
!  mmolec is the entity number in that elementary step which remains fixed during 
!         the pattern matching.
!  patrnnmxcoord is the max coordination number encountered in any pattern
!  patrnlistnsites(i1) gives the number of sites for pattern i1
!  patrnlistspecies(i1,i2) gives the i2^th molecule species number involved
!         in pattern i1. If i2 = 0 it gives the total number of molecules.
!  patrnliststype(i1,i2) gives the site type of site i2 in pattern i1
!  patrnlistneigh(i1,i2,i3) gives the i3^th neighbor of the i2^th site
!         involved in pattern i1. If i3 = 0 it gives the total number of
!         neighbors for site i2.
!  patrnlistlatticestate(i1,i2,i3) gives the lattice state (coverage) for the 
!         pattern i1. For i3 = 1 it gives the entity number occupying site i2, 
!         and for i3 it gives the dentate number
!  patrnlistadsorbpos(i1,i2,i3) gives the site occupied by dentate i3 of 
!         molecule i2 for the pattern i1. It is the inverse mapping of 
!         patrnlistlatticestate
!  nneigh is the number of sites that mmolec occupies as well as those that are 
!         the former neighboring.
!  neighlist is the list of the sites just mentioned. Permutations of these sites 
!         will be considered during the pattern matching procedure.
!  adsexclude(1:) is a list of adsorbates that must be excluded as the mmolec-1, etc entities
!         while searching for patterns. adsexclude(0) gives the number of these adsorbates. This 
!         is done to avoid double counting processing upon adding them.
!  cntrpatrn is an output variable giving the number of valid patterns found.
!  validpatrns is a 2-dimensional array, each row of which gives a valid pattern as
!         a list of lattice sites which the pattern sites are mapped to.
!
! The subroutine find_valid_patterns that is subsequently invoked works by finding the permutations 
! of nneigh sites that match the graph pattern as well as the site-type requirements and coverage 
! (namely entity, species, dentate)
!
! Example: suppose that an adsorabate occupying site 14 (the first ident sites mentioned in
! neighlist are always the ones that the fixed adsorbate occupies). Thus:
!
!                      1   2   3   4   5   6   7
!       neighlist = [ 14  12  13  15  16  21  23]
!
! Say that jelemstep is a 2-site pattern. The subroutine goes through the pattern macthing
! procedures and determines that there are two valid patterns cntrpatrn = 2). The ith row 
! of validpatrns gives the ith permutation of the sites contained in neighlist, for instance 
! validpatrns could be:
!
!       validpatrns = [ 14  12
!                       14  23 ]

npatrnsites = patrnlistnsites(jpattern)
kspec = patrnlistspecies(jpattern,mmolec)
ident = surfspecsdent(kspec)

! *** Prepare the input for the pattern matching subroutine

! We want pattern matches in which molecule mmolec is fixed on the sites it currently
! occupies on the lattice, thus set fixedsites

workarray => get_work_array()
workarray%fixedsites = 0
do k = 1,ident
    workarray%fixedsites(patrnlistadsorbpos(jpattern,mmolec,k)) = k
enddo

! The neighboring, site-type and state pattern to be searched is the one corresponding to
! the elementary step jelemstep

workarray%patrnstate(:npatrnsites, 1:3:2)                                      &
   = patrnlistlatticestate(jpattern,:npatrnsites,:)
do i = 1,npatrnsites
    if (patrnlistlatticestate(jpattern,i,1) /= -1) then
        workarray%patrnstate(i,2) = patrnlistspecies(jpattern,patrnlistlatticestate(jpattern,i,1))
    else
        workarray%patrnstate(i,1) = -1
        workarray%patrnstate(i,3) = -1
    endif
enddo
!write(*,*) workarray%patrnstate

checkneigh = .true.
checkgeom = .true.
checkstate = .true.
checkstypes = .true.

! *** Finally find the valid patterns

call find_valid_patterns(                                  &
                nneigh, neighlist(1:nneigh),               & ! the number of candidate sites and the vector containing them
                npatrnsites,                               &
                patrnlistneigh(jpattern, 1:npatrnsites,0:patrnnmxcoord), & ! number of sites in the pattern, neighboring structure
                patrnliststype(jpattern,:),                &
                workarray%patrnstate(:npatrnsites, :),     & ! site types and state in the pattern
                patrnlistnanglessequenc(jpattern, :),      &
                patrnlistangles(jpattern,:),               & ! number of and sites in angle specifications and angle values
                patrnlistnomirrors(jpattern),              & ! allow for mirror images or not
                patrnlistabslorientat(jpattern),           & ! require absolute orientation or not
                patrnlistorientationedge(jpattern,:),      &
                partnlistorientationangles(jpattern),      &
                workarray%fixedsites,adsexclude,mmolec,    &
                checkneigh,checkgeom,checkstate,checkstypes, &
                cntrpatrn,validpatrns)
return

end subroutine find_pattern_match

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module lattice_handle_module
