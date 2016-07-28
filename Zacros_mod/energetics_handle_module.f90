! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module energetics_handle_module

implicit none

integer nglobclust
integer globclustcapacity0

!integer, allocatable :: event_times_labels(:)
!integer, allocatable :: event_times_indexes(:)

real(8) globalenergy
real(8) energconst

integer, allocatable :: globclustypesites(:,:)
integer, allocatable :: nadsorglobclusparticip(:)
integer, allocatable :: adsorglobclusparticip(:,:,:)

!real(8), allocatable :: event_times_heap(:)
!real(8), allocatable :: procpropens0(:)

! Declares a type which is singularly used to hold work arrays.
! The issue is that some functions are called over and over. We do not want to
! have to allocate and reallocate those arrays.  Currently, only holds arrays
! previously 
type work_array_type
  integer, allocatable :: startsites(:)
  integer, allocatable :: neighlist(:)
  integer, allocatable :: validpatrns(:,:)
  integer, allocatable :: indxs(:)
end type 

! A private work array.
! When/if we move to multi-processing, we may need more than one. 
type (work_array_type), private, target:: prv_work_array
private :: work_array_type, get_work_array

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  function get_work_array()
    ! This function returns a pointer to an allocated work array.  In
    ! multiprocessing, we may need to more than one work array. We will only
    ! need modify this function and add a locking mechanism.
    use energetics_setup_module, only: clusterlevels, clusternmxsites
    use lattice_setup_module, only: maxcoord

    type(work_array_type), pointer :: get_work_array

    if( .not. allocated(prv_work_array%startsites)) then
      allocate(prv_work_array%startsites(clusternmxsites))
      allocate(prv_work_array%neighlist(100*maxcoord))
      allocate(prv_work_array%validpatrns(100*maxcoord,clusternmxsites))
      allocate(prv_work_array%indxs(0:maxval(clusterlevels)))
    endif

    get_work_array => prv_work_array
  end function get_work_array

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine cleanup_ehm
    ! Deletes work arrays

    if(allocated(prv_work_array%startsites)) then
      deallocate(prv_work_array%startsites)
      deallocate(prv_work_array%neighlist)
      deallocate(prv_work_array%validpatrns)
      deallocate(prv_work_array%indxs)
    endif 
    if(allocated(globclustypesites))      deallocate(globclustypesites)
    if(allocated(nadsorglobclusparticip)) deallocate(nadsorglobclusparticip)
    if(allocated(adsorglobclusparticip))  deallocate(adsorglobclusparticip)

  end subroutine cleanup_ehm 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine initialize_energetics()

use simulation_setup_module, only: debug_report_globenerg
use energetics_setup_module, only: clusterenrg, clusterlatticestate,           &
                                   clusternsites, clusternmxsites, nclusters, clusterOcc, clustergraphmultipl
use simulation_setup_module, only: debug_check_globenerg
use constants_module, only: iglbenergdbg
use lattice_setup_module, only: nsites, maxcoord
use state_setup_module, only: nadsorb
use parser_module, only: dbl2str, int2str

implicit none

integer i, j, k, mmolec
integer, allocatable :: adsexclude(:)

!call find_cluster_graphmultiplicities()

globclustcapacity0 = 50*nclusters*maxcoord*nsites
nglobclust = 0
globalenergy = 0.d0

! Array globclustypesites(i1,i2) gives: for i2 = 0 the cluster number for
! global-cluster i1 and for i2 = 1:clusternsites the mapping between the pattern
! sites and the actual lattice sites for this global-cluster
allocate(globclustypesites(globclustcapacity0,0:clusternmxsites))
! nadsorglobclusparticip(i1) gives the number of global-clusters in which
! adsorbate/entity i1 participates
allocate(nadsorglobclusparticip(nsites))
! adsorglobclusparticip(i1,i2,i3) gives: for i2 = 1:nadsorglobclusparticip(i1)
! and i3 = 1 the indexes of the global-clusters in which adsorbate/entity i1 participates.
! For i3 = 2 the number with which that adsorbate appears in the cluster
allocate(adsorglobclusparticip(nsites,50*nclusters*maxcoord,2))

! Initialize variables
do i = 1,globclustcapacity0
    do j = 0,clusternmxsites
        globclustypesites(i,j) = 0
    enddo
enddo
do i = 1,nsites
    nadsorglobclusparticip(i) = 0
enddo
do i = 1,nsites
    do j = 1,10*maxcoord
        do k = 1,2
            adsorglobclusparticip(i,j,k) = 0
        enddo
    enddo
enddo

if (debug_report_globenerg) then
    write(iglbenergdbg,'(a)') 'Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
endif

allocate(adsexclude(0:1))

! If an empty cluster has been specified (careful: this is not an unoccupied
! site) then find the value of the constant energy:
energconst = 0.d0
do i = 1,nclusters
    if (clusternsites(i) == 1 .and. clusterlatticestate(i,1,1) == -1) then
        energconst = nsites*clusterenrg(i)
    endif
enddo
globalenergy = globalenergy + energconst

if (debug_report_globenerg) then
    write(iglbenergdbg,'(a)') 'Total empty-cluster energy constant = '         &
                              // dbl2str(energconst)
endif

! Loop over all molecules (and empty site entities) that exist on the lattice
do i = 1,nadsorb
    
    ! We want processes in which adsorbate/entity i participates as the 1st
    ! molecule
    mmolec = 1
    
    call add_species_related_energy_contributions(i,adsexclude,mmolec)
    
enddo
deallocate(adsexclude)

if (debug_check_globenerg) then
    call check_energetics()
endif

! Report global lattice energy (for debugging)
if (debug_report_globenerg) then
    write(iglbenergdbg,'(a)') 'Current total lattice energy is ' // dbl2str(globalenergy)
endif

return

end subroutine initialize_energetics
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine add_species_related_energy_contributions(i,adsexclude,mmolec)

use simulation_setup_module, only: debug_report_globenerg, surfspecsdent,      &
                                   surfspecsnames
use lattice_handle_module, only: adsorbspecposi, find_pattern_match
use lattice_setup_module, only: sitetype, siteneighb1
use energetics_setup_module, only: clusterangles, specclusterparticip,         &
                                   clusternsites, clusteradsorbpos,            &
                                   clusternomirrorimgs,                        &
                                     clusterabslorientat,                      &
                                       clusterorientationedge,                 &
                                         clusterorientationangles,             &
                                   clusterstype, clusterlevels,                &
                                   clustergraphmultipl, clusternames,          &
                                   clusterspecs, nspecclusterparticip,         &
                                   clusterenrg, clusterlatticestate,           &
                                   clusternanglessequenc, clusterneigh,        &
                                   clusternmxangles, clusternmxcoord,          &
                                   clusternames, clusterOcc
use parser_module, only: int2str, dbl2str
use constants_module, only: iglbenergdbg
use graph_functions_module, only: find_dlevel_neighbors

implicit none
integer, intent(in) :: i, mmolec
integer, intent(in) :: adsexclude(0:)

integer ipatrnsite, ilattcsite, j, jcluster, k, kspec
integer m, nspat, dmax, nstart, nneigh, cntrpatrn


! A work array to avoid allocation/deallocation in much called routine.
type(work_array_type), pointer:: workarray
workarray => get_work_array()

! Add all valid energetic contributions in which adsorbate i participates

! Note mstam 09-Nov-2011: this algorithm is not appropriate for small lattices
! since in such cases the same site appears more than one in the neighboring
! list of a site. This algorithm, however, works under the assumption that there
! is no such multiplicity, and thus, will miss patterns that involve the image
! of a site on the periodic structure. A resolution to this problem would be to
! keep a larger lattice in the memory in which the images are accounted for
! explicitly. The KMC operations would be performed on the small lattice and the
! sites that correspond to images of sites therein, will be updated
! automatically. The energetics will be calculated correctly with the aid of the
! explicit images.

kspec = adsorbspecposi(i,0)

! Loop over all clusters to which species kspec participates
clusters: do j = 1, nspecclusterparticip(kspec)
    
    jcluster = specclusterparticip(kspec,j,1)
    nspat = clusternsites(jcluster)
    
!    tmpstr = clusternames(jcluster)
!    if (striccompare(trim(tmpstr),'o_pair_4nn')) then
!        continue
!    endif
    
    ! *** INITIAL CHECKS ***
       
    if (specclusterparticip(kspec, j, 2) /= mmolec) cycle clusters
    
    ! We will immediately check to see whether molecule i occupies on the
    ! lattice the site types that the pattern jcluster requires        
    do k = 1,surfspecsdent(kspec)
        ipatrnsite = clusteradsorbpos(jcluster, mmolec, k) 
        ilattcsite = adsorbspecposi(i, k)
        m = clusterstype(jcluster, ipatrnsite)
        if ( m > 0 .and. m /= sitetype(ilattcsite) ) cycle clusters
    enddo
    
    
    ! *** CANDIDATE SITES ***
    
    ! Make a list of candidate sites of the lattice. These sites could be
    ! participating in a global cluster pertaining to cluster jcluster.
    ! dmax = max level of cluster i starting from sites occupied by molecule i
    dmax = clusterlevels(jcluster,mmolec) 
    
    nstart = surfspecsdent(kspec)
    workarray%startsites(1:nstart) = adsorbspecposi(i,1:nstart)
    
    call find_dlevel_neighbors( siteneighb1, dmax, nstart,                     &
                                workarray%startsites, nneigh,                  &
                                workarray%indxs(0:dmax),                       &
                                workarray%neighlist )
        
     ! write(*,('(a)')) 'Entity ' // trim(int2str(i)) // ' (species ' // trim(surfspecsnames(kspec)) // &
     ! ')  could participate in cluster ' // trim(int2str(specclusterparticip(kspec,j,1))) // &
     ! ' involving some permutation of ' // trim(int2str(nspat)) // & 
     ! ' site(s) out of: '
     ! write(*,('(' // int2str(nneigh) // '(I4,1x))')) (workarray%neighlist(k),k=1,nneigh)
                        
    ! *** FINDING VALID PATTERNS ***
                      
    call find_pattern_match( &
                jcluster,mmolec,clusternmxcoord,clusternmxangles, &
                clusternsites,clusterspecs,clusterstype, &
                clusterneigh,clusterlatticestate,clusteradsorbpos, &
                clusternanglessequenc,clusterangles, &
                clusternomirrorimgs, &
                clusterabslorientat, &
                    clusterorientationedge, &
                        clusterorientationangles, &
                ! the number of candidate sites and the vector containing them
                nneigh, workarray%neighlist, & 
                adsexclude, &
                cntrpatrn, workarray%validpatrns)
    
    ! write(*,*) 'Valid permutations:'
    ! do k = 1,cntrpatrn
    !     write(*,*) (validpatrns(k,m), m = 1,nspat)
    ! enddo
    ! continue
        
    ! *** ADDING IDENTIFIED ENERGETIC CONTRIBUTIONS ***
    
    patternloop: do k = 1,cntrpatrn ! For every pattern...
        
        ! At this point we have identified a valid cluster that will have a
        ! contribution on the overall energy of the system.
        
        call add_global_cluster(jcluster,nspat, workarray%validpatrns(k,:))
        
        if (debug_report_globenerg) then
            write(iglbenergdbg,'(a)') 'Global-cluster '                        &
                                      // trim(int2str(nglobclust))             &
                                      // ' identified:'
            write(iglbenergdbg,'(a)') '    Cluster number: '                   &
                                      // trim(int2str(jcluster))
            write(iglbenergdbg,'(a)') '    Cluster description: '              &
                                      // trim(clusternames(jcluster))
            write(iglbenergdbg,'(a)',advance='no')                             &
              '    Mapping of lattice to pattern sites: '
            write(iglbenergdbg,'(' // int2str(nspat) // '(I4,1x))')            &
              ( workarray%validpatrns(k,m),m = 1,nspat)
            write(iglbenergdbg,'(a)') '    Cluster graph-multiplicity: '       &

                              // trim(int2str(clustergraphmultipl(jcluster)))
            write(iglbenergdbg,'(a)') '    Its energy contribution is '        &
                   // trim( dbl2str( clusterenrg(jcluster)                     &
                                     / clustergraphmultipl(jcluster) ) )
        endif

		clusterOcc(jcluster) = clusterOcc(jcluster) + 1
		
    enddo patternloop
    
enddo clusters ! loop over all clusters for species kspec

return

end subroutine add_species_related_energy_contributions

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine add_global_cluster(jcluster,nspat,validpattern)

use energetics_setup_module, only: clusteradsorbpos, clusterspecs,             &
                                   clustergraphmultipl, clusterenrg
use lattice_handle_module, only: latticestate
implicit none

integer, intent(in) :: jcluster, nspat
integer, intent(in) :: validpattern(:)

integer ipatrnsite, ilattcsite, m, p

nglobclust = nglobclust + 1

! Add the new contribution to the global energy
globalenergy = globalenergy + clusterenrg(jcluster)                            &
                              / clustergraphmultipl(jcluster)

! In the array that gives the sites mapping, add the new process in the list 
globclustypesites(nglobclust,0) = jcluster
globclustypesites(nglobclust,1:nspat) = validpattern(1:nspat)

! In the array that gives the global clusters in which a molecule/entity
! participates, add the new process in the list of the corresponding entity.
do m = 1,clusterspecs(jcluster,0) ! For all reactant molecules...
    ipatrnsite = clusteradsorbpos(jcluster,m,1)
    ! find the lattice site where the first dentate sits...
    ilattcsite = validpattern(ipatrnsite) 
    ! and set p equal to the molecule/entity number
    p = latticestate(ilattcsite,1) 
    nadsorglobclusparticip(p) = nadsorglobclusparticip(p) + 1
    ! global cluster in which molecule participates
    adsorglobclusparticip(p,nadsorglobclusparticip(p),1) = nglobclust 
    ! it participates as the mth molecule
    adsorglobclusparticip(p,nadsorglobclusparticip(p),2) = m 
enddo            

return

end subroutine add_global_cluster

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine remove_energetics(energrem)

use energetics_setup_module, only: clustergraphmultipl, clusterenrg,          &
                                   clusteradsorbpos, clusternsites,           &
                                   clusterspecs, clusterOcc
use simulation_setup_module, only: debug_report_globenerg
use constants_module, only: iglbenergdbg
use lattice_handle_module, only: latticestate
use parser_module, only: int2str

implicit none
integer, intent(inout) :: energrem(0:)

integer i, iclustrem, iadsorb, ilattcsite, j, jcluster, nspat, nclusts, m, k
integer nlastcluster

continue

do i = 1,energrem(0)

    nlastcluster = nglobclust

    iclustrem = energrem(i)
    
    jcluster = globclustypesites(iclustrem,0)
    nspat = clusternsites(jcluster)

    ! Remove that energy contribution thereby updating the total energy value
    globalenergy = globalenergy                                                &
                   - clusterenrg(jcluster)/clustergraphmultipl(jcluster)
    ! Decrease the number of global clusters by one
    nglobclust = nglobclust - 1
    
    ! Remove that cluster from the array which gives the clusters in which a
    ! molecule/entity participates
    do m = 1,clusterspecs(jcluster,0) ! For all participating molecules...

        ! find the lattice site in which the 1st dentate sits
        ilattcsite = globclustypesites( iclustrem,                             &
                                        clusteradsorbpos(jcluster,m,1) )
        
        ! then find which adsorbate this is on the lattice
        iadsorb = latticestate(ilattcsite,1)
        nclusts = nadsorglobclusparticip(iadsorb)
        
        do k = 1,nclusts
            if (adsorglobclusparticip(iadsorb,k,1) == iclustrem) then
                if (k < nclusts) then
                    adsorglobclusparticip(iadsorb,k,1)                         &
                       = adsorglobclusparticip(iadsorb,nclusts,1)
                    adsorglobclusparticip(iadsorb,k,2)                         &
                       = adsorglobclusparticip(iadsorb,nclusts,2)
                endif
                nadsorglobclusparticip(iadsorb)                                &
                  = nadsorglobclusparticip(iadsorb) - 1
            endif
        enddo

    enddo            
    
    ! In the array that gives the sites mapping, remove that process from the
    ! list 
    globclustypesites(iclustrem,0:nspat) = -1
    
    if (debug_report_globenerg) then
        write(iglbenergdbg,'(a)') 'Cluster ' // trim(int2str(iclustrem))       &
                                  // ' was removed.'
    endif
    clusterOcc(jcluster) = clusterOcc(jcluster) - 1
	
    ! Relabel if necessary
    if (iclustrem < nlastcluster) then
        
        jcluster = globclustypesites(nlastcluster,0)
        nspat = clusternsites(jcluster)

        ! The last cluster will take the place of cluster with id iclustrem
        do m = 0,nspat
            globclustypesites(iclustrem,m) = globclustypesites(nlastcluster,m)
        enddo
        
        do m = 1,clusterspecs(jcluster,0) ! For all participating molecules...
            
            ! find the lattice site in which the 1st dentate sits
            ilattcsite = globclustypesites( nlastcluster,                      &
                                            clusteradsorbpos(jcluster,m,1) )
             
            ! then find which adsorbate this is on the lattice
            iadsorb = latticestate(ilattcsite,1)
            nclusts = nadsorglobclusparticip(iadsorb)
            
            do k = 1,nclusts             
                
                ! Change the number of cluster nclusterprev in the adsorbate
                ! participation array adsorglobclusparticip
                if (adsorglobclusparticip(iadsorb,k,1) == nlastcluster) then
                    adsorglobclusparticip(iadsorb,k,1) = iclustrem
                endif
                
            enddo
        
        enddo
                
        if (debug_report_globenerg) then
            write(iglbenergdbg,'(a)') 'Cluster '                               &
                                      // trim(int2str(nlastcluster))           &
                                      // ' was relabeled to '                  &
                                      // trim(int2str(iclustrem)) // '.'
        endif

        ! If the relabelling affects the list of subsequent entities to be
        ! removed, take care of it 
        do j = i+1,energrem(0)
            if (energrem(j) == nlastcluster) energrem(j) = iclustrem
        enddo
        
        ! Decrease the number of the last cluster stored by one
        nlastcluster = nlastcluster - 1
            
    endif
        
enddo                       

!call check_energetics()

return

end subroutine remove_energetics

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_energetics()

use parser_module, only: int2str
use error_module, only: moreinfo, error
use state_setup_module, only: nadsorb
use energetics_setup_module, only: clusteradsorbpos, clusterspecs, clusternames
use lattice_handle_module, only: adsorbspecposi, latticestate
use simulation_setup_module, only: surfspecsdent
implicit none

integer i, ispec, ilatsite, iadsorb, j, k, m, p, jcluster

logical foundproc

! First, check the following:
! for each global cluster, the adsorbates involved must have that cluster listed
! in their participation list
do i = 1,nglobclust ! For each cluster
    
    ! find which cluster (type) it corresponds to
    jcluster = globclustypesites(i,0) 
    
    if (jcluster <= 0) cycle
    
    ! for every participating molecule of that cluster
    do j = 1,clusterspecs(jcluster,0) 
    
        ! find the lattice site in which the 1st dentate sits
        ilatsite = globclustypesites(i,clusteradsorbpos(jcluster,j,1))
        
        ! then find which adsorbate this is on the lattice
        iadsorb = latticestate(ilatsite,1)
        
        ! Check that the mapping between cluster sites and lattice sites
        ! of that adsorbate is correct
        do m = 1,surfspecsdent(clusterspecs(jcluster,j))
            ! lattice site in which 
            k = globclustypesites(i,clusteradsorbpos(jcluster,j,m)) 
            ! the mth dentate of adsorbate iadsorb sits
            if (iadsorb /= latticestate(k,1)) then
                
                moreinfo = 'Adsorbate ' // trim(int2str(iadsorb))              &
                            // ' is referenced as occupying different sites '  &
                            // 'in cluster ' // trim(clusternames(jcluster))  &
                            // ' (global cluster '// trim(int2str(i))          &
                            // ') than those it really occupies on the lattice.'
                call error(5013)                
                                
            endif
        enddo
    
        ! Check that this process appears in the process-participation list of
        ! that adsorbate
        foundproc = .false.
        do k = 1,nadsorglobclusparticip(iadsorb)
            if (adsorglobclusparticip(iadsorb,k,1) == i .and. &
                adsorglobclusparticip(iadsorb,k,2) == j) then
                foundproc = .true.
                exit
            endif
        enddo
        
        if (.not. foundproc) then

            moreinfo = 'Cluster ' // trim(int2str(i))                          &
                       // ' involves adsorbate ' // trim(int2str(iadsorb))     &
                       // ' which does not list that process in its '          &
                       // ' participation list.'
            call error(5014)                

        endif
                
    enddo
    
enddo

! Second, check the following:
! for each adsorbate on the lattice, each of the clusters listed must contain
! that molecule as a participating entity
do iadsorb = 1,nadsorb ! For each adsorbate
    ispec = adsorbspecposi(iadsorb,0)
    
    ! for each cluster in which it participates
    do j = 1,nadsorglobclusparticip(iadsorb) 
    
        ! the cluster...
        i = adsorglobclusparticip(iadsorb,j,1) 
        ! ... in which iadsorb participates as the pth molecule
        p = adsorglobclusparticip(iadsorb,j,2) 
    
        ! find which cluster (type) it corresponds to
        jcluster = globclustypesites(i,0) 
        
        do k = 1,surfspecsdent(ispec)
            ! lattice site
            ilatsite = globclustypesites(i,clusteradsorbpos(jcluster,p,k)) 
           
            if (adsorbspecposi(iadsorb,k) /= ilatsite) then
            
                moreinfo = 'Cluster ' // trim(int2str(i))                      &
                           // ' lists adsorbate ' // trim(int2str(iadsorb))    &
                           // ' at different sites than those it actually '    &
                           // 'occupies.'
                call error(5015)                

            endif                
           
            if (latticestate(ilatsite,1) /= iadsorb) then

                moreinfo = 'Cluster ' // trim(int2str(i))                      &
                           // ' involves adsorbate ' // trim(int2str(iadsorb)) &
                           // ' as adsorbate '                                 &
                           // trim(int2str(latticestate(ilatsite,1))) // '.'
                call error(5016)                

            endif
           
            if (latticestate(ilatsite,2) /= ispec) then

                moreinfo = 'Cluster ' // trim(int2str(i))                     &
                           // ' lists adsorbate ' // trim(int2str(iadsorb))   &
                           // ' as species '                                  &
                           // trim(int2str(latticestate(ilatsite,2))) // '.'
                call error(5017)                

            endif
            
            if (latticestate(ilatsite,3) /= k) then

                moreinfo = 'Cluster ' // trim(int2str(i))                     &
                           // ' improperly lists adsorbate '                  &
                           // trim(int2str(iadsorb)) // '.'
                call error(5018)                

            endif
            
        enddo
        
    enddo
    
enddo

return

end subroutine check_energetics

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module energetics_handle_module


