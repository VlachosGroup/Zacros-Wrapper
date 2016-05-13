! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module energetics_setup_module

use constants_module, only: nnam0

implicit none

integer nclusters, clusternmxsites, clusternmxcoord, clusternmxangles, clustermaxlevel

integer, allocatable :: clusternsites(:)
integer, allocatable :: clusterspecs(:,:)
integer, allocatable :: clusterneigh(:,:,:)
integer, allocatable :: clusterstype(:,:)
integer, allocatable :: clusternanglessequenc(:,:)
integer, allocatable :: clusterorientationedge(:,:)
integer, allocatable :: clusterlatticestate(:,:,:)
integer, allocatable :: clusteradsorbpos(:,:,:)

integer, allocatable :: nspecclusterparticip(:)
integer, allocatable :: clustergraphmultipl(:)
integer, allocatable :: specclusterparticip(:,:,:)

integer, allocatable :: clusterlevels(:,:)

real(8), allocatable :: clusterenrg(:)
real(8), allocatable :: clusterangles(:,:)
real(8), allocatable :: clusterorientationangles(:)

character(nnam0), allocatable :: clusternames(:)

logical, allocatable :: clusternomirrorimgs(:)
logical, allocatable :: clusterabslorientat(:)

integer, allocatable :: clusterOcc(:)
contains

  subroutine cleanup_esm
    ! Cleanup module globals

    if(allocated(clusternsites))         deallocate(clusternsites)
    if(allocated(clusterspecs))          deallocate(clusterspecs)
    if(allocated(clusterneigh))          deallocate(clusterneigh)
    if(allocated(clusterstype))          deallocate(clusterstype)
    if(allocated(clusternanglessequenc)) deallocate(clusternanglessequenc)
    if(allocated(clusterlatticestate))   deallocate(clusterlatticestate)
    if(allocated(clusteradsorbpos))      deallocate(clusteradsorbpos)
    if(allocated(nspecclusterparticip))  deallocate(nspecclusterparticip)
    if(allocated(clustergraphmultipl))   deallocate(clustergraphmultipl)
    if(allocated(specclusterparticip))   deallocate(specclusterparticip)
    if(allocated(clusterlevels))         deallocate(clusterlevels)
    if(allocated(clusterenrg))           deallocate(clusterenrg)
    if(allocated(clusterangles))         deallocate(clusterangles)
    if(allocated(clusternames))          deallocate(clusternames)
    if(allocated(clusternomirrorimgs))   deallocate(clusternomirrorimgs)

  end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine read_energetics_setup()

use parser_module, only: int2str, striccompare, dbl2str, iwrite, break_words,  &
                         getrecord
use constants_module, only: maxwords2, lengthrecinp, remchar, csimfname,       &
                            cgenoutfname, cenergfname, ienergread
use error_module, only: moreinfo, error, warning
use simulation_setup_module, only: surfspecsdent, surfspecsnames, nsurfspecs,  &
                                   maxdent
use graph_functions_module, only: find_max_level

implicit none

integer i, iadsorb, identat, io, irec, icluster, j, m, mmolec
integer k, nwords, nstart, kspec, kdent, dmax, q, setvbuf3f
integer maxcharsname, iw(maxwords2)

integer, allocatable :: startsites(:)
!integer, allocatable :: uniquespecbindingconfig(:,:)

character(lengthrecinp) recinput

logical genoutop, clusteropen, energopen, defaultnaming, energfex

! Parses the energetics setup file energetics_input.dat

! Level-1 keywords allowed in this subroutine:
!    energetics
!    end_energetics
!    cluster

! The program first performs a parsing of the entire input to check whether the 
! blocks are opened and closed properly and collect information that will be used 
! for determining the bounds of the arrays used for energetics specification.

! Check if the energetics specification file exists
energfex = .false.
inquire(file=cenergfname,exist=energfex)
if (.not.energfex) then
    moreinfo = 'Energetics specification file ' // trim(cenergfname) // ' does not exist.'
    call error(1)
endif

inquire(iwrite,opened=genoutop)
if (.not.genoutop) then
    open(unit=iwrite,file=trim(cgenoutfname),status='unknown',position='append')
    q = setvbuf3f(iwrite,1,0) ! set line-buffered behavior
endif

write(iwrite,'(/,a)') 'Energetics setup:'
write(iwrite,'(a)')   '~~~~~~~~~~~~~~~~~'

open(unit=ienergread,file=trim(cenergfname),status='old')

! Perform an initial pass to check for proper syntax and opening and closing
! of blocks and find the following:
!    number of clusters: nclusters
!    the max number of sites involved in a cluster: clusternmxsites
!    the max coordination number of sites involved in a cluster: clusternmxcoord

icluster = 0
energopen = .false.
clusteropen = .false.
nclusters = 0
clusternmxsites = 0
clusternmxcoord = 0
clusternmxangles = 0

irec = 0
io = 0

call getrecord(ienergread,recinput,irec,io)

do while (io >= 0)
    
    call break_words(recinput,' ',remchar,nwords,iw)
    
    if (nwords == 0) then
        call getrecord(ienergread,recinput,irec,io)
        cycle
    endif        
        
    if (striccompare(recinput(iw(1):iw(2)),'energetics')) then
    
        if (energopen) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
            call error(100)
        endif
    
        energopen = .true.

    elseif (striccompare(recinput(iw(1):iw(2)),'end_energetics')) then
    
        energopen = .false.
        exit

    elseif ( striccompare(recinput(iw(1):iw(2)),'cluster') ) then
        
        if (.not. energopen) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(8001)
        endif
        
        call parse_cluster(.true.,recinput,nwords,irec,io,icluster,iw,clusteropen)            
                        
    else
        
        moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
        call error(8100)

    endif

    call getrecord(ienergread,recinput,irec,io)

enddo

if (clusteropen) then
    moreinfo = 'Cluster information not completed in line ' // trim(int2str(irec)) // '.'
    call error(8024)
endif
if (energopen) then
    call warning(8100)
endif
if (nclusters == 0) then
    moreinfo = 'File does not appear to contain any cluster information.'
    call error(8001)
endif

! clusternames(i1) gives the name of cluster i1
allocate(clusternames(nclusters))

! clusternsites(i1) gives the number of sites involved in cluster i1
allocate(clusternsites(nclusters))

! clusterspecs(i1,i2) gives the i2^th molecule species number involved
!   in cluster i1. If i2 = 0 it gives the total number of molecules.
allocate(clusterspecs(nclusters,0:clusternmxsites))

! clusterneigh(i1,i2,i3) gives the i3^th neighbor of the i2^th site
!   involved in cluster i1. If i3 = 0 it gives the total number of
!   neighbors for site i2. Note that this structure is the subgraph (cluster pattern)
!   that will be searched for in the lattice, thereby identifying an energetic contribution.
allocate(clusterneigh(nclusters,clusternmxsites,0:clusternmxcoord))
allocate(clusterstype(nclusters,clusternmxsites))

! clusterlatticestate(i1,i2,i3) gives the lattice state (coverage) for the 
!   cluster i1. For i3 = 1 it gives the entity number occupying site i2, 
!   and for i3 it gives the dentate number
allocate(clusterlatticestate(nclusters,clusternmxsites,2))

! clusteradsorbpos(i1,i2,i3) gives the site occupied by dentate i3 of 
! molecule i2 for the cluster i1. It is the inverse mapping of 
! clusterlatticestate
allocate(clusteradsorbpos(nclusters,clusternmxsites,maxdent))

! clusterenrg(i1) gives the energetic contribution for cluster i1
allocate(clusterenrg(nclusters))

! nspecclusterparticip(i1) gives the number of clusters in which species i1 
!   participates. If a species participates in one cluster two times, this is 
!   doubly counted (or multiply-counted in general)
allocate(nspecclusterparticip(0:nsurfspecs))
! specclusterparticip(i1,i2,i3) for i3 = 1 it gives the i2^th cluster in 
!   which species i1 participates, and for i3 = 2 it gives the molecule number 
!   with which this species appears in that cluster
allocate(specclusterparticip(0:nsurfspecs,clusternmxsites*nclusters,2))

! clustergraphmultipl(i1) gives the stoichiometric multiplicity of cluster i1
allocate(clustergraphmultipl(nclusters))

! clusterSC(i1) gives the occurrences of each cluster
allocate(clusterOcc(nclusters))
! clusterlevels(i1,i2) gives the graph level of the pattern corresponding 
!   to the cluster i1 if we start from the sites of molecule i2
allocate(clusterlevels(nclusters,clusternmxsites))

! clusternanglessequenc(i1,3*(i2-2):3*i2) gives the site sequence i2 in angle definitions for cluster i1
! It gives the number of such defined sequences for the second argument = 0
allocate(clusternanglessequenc(nclusters,0:3*clusternmxangles))
! clusterangles(i1) gives the angle values for the angle definitions for cluster i1
allocate(clusterangles(nclusters,clusternmxangles))

! clusternomirrorimgs(i1) is .false. by default allowing mirror images of the cluster i1 to be detected. A value 
! of .true. prevents mirror images from being detected
allocate(clusternomirrorimgs(nclusters))

! clusterabslorientat(i1) is .true. if absolute orientation is required for cluster i1. Combined with
! the prevention of mirror image detection, this feature allows for only translations of the pattern to
! be detected. By default clusterabslorientat is .false.
allocate(clusterabslorientat(nclusters))
! clusterabslorientat(i1,1:2) gives the vertexes of the graph that serve as the "orientation key". The
! code calculates the angle between vector (x2-x1,y2-y1) and (1,0) and checks whether it is equal to that
! specified in clusterorientationangles(i1) that refers to the cluster i1.
allocate(clusterorientationedge(nclusters,2))
allocate(clusterorientationangles(nclusters))

! Initialization of variables
defaultnaming = .false.
clusteropen = .false.
energopen = .false.
clustermaxlevel = 0
icluster = 0
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

! The first pass has not found any problems in the syntax. Begin parsing the mechanism

rewind(ienergread)
irec = 0
io = 0

do while (io >= 0)

    call getrecord(ienergread,recinput,irec,io)
    
    call break_words(recinput,' ',remchar,nwords,iw)
    
    if (nwords == 0) cycle
            
    if (striccompare(recinput(iw(1):iw(2)),'end_energetics')) then

        exit
        
    elseif (striccompare(recinput(iw(1):iw(2)),'cluster')) then
    
        call parse_cluster(.false.,recinput,nwords,irec,io,icluster,iw,clusteropen)
                
    endif
        
enddo

close(ienergread)

! From the energetics specification data create arrays:
!   nclusterparticip
!   clusterparticip
!   clusteradsorbpos

allocate(startsites(clusternmxsites))

! Construct the species'-participation and pattern-levels arrays
do kspec = 0,nsurfspecs
    nspecclusterparticip(kspec) = 0
enddo

do i = 1,nclusters

    ! Create the inverse mappings that give the positions of all dentates of the 
    ! species participating in the cluster
    do j = 1,clusternsites(i) ! For every site of the pattern
        
        if (clusterlatticestate(i,j,1) == -1) then
            cycle
        endif
        
        ! Lattice state
        iadsorb = clusterlatticestate(i,j,1) ! molecule number
        if (iadsorb <= 0 .or. iadsorb > clusterspecs(i,0)) then
            moreinfo = 'The entity number '// trim(int2str(iadsorb)) // ' is out of range 1:' // &
                       trim(int2str(clusterspecs(i,0))) // & 
                       ' in the lattice state of cluster ' // trim(int2str(i)) // &
                       ' (' // trim(clusternames(i)) // ').'
            call error(8025)
        endif
        identat = clusterlatticestate(i,j,2) ! molecule's dentate number
        if (clusteradsorbpos(i,iadsorb,identat) .ne. 0) then
            moreinfo = 'Dentate ' // trim(int2str(identat)) // ' of entity number '// trim(int2str(iadsorb)) // &
                       ' has already been assigned site ' // trim(int2str(clusteradsorbpos(i,iadsorb,identat))) // &
                       ' in the lattice state of cluster ' // trim(int2str(i)) // &
                       ' (' // trim(clusternames(i)) // ').'
            call error(8029)        
        endif
        clusteradsorbpos(i,iadsorb,identat) = j
                
    enddo
    
    ! Check that all dentates of each species have been assigned a site
    do j = 1,clusterspecs(i,0)
        kspec = clusterspecs(i,j)
        do kdent = 1,surfspecsdent(kspec)
            if (clusteradsorbpos(i,j,kdent) == 0) then
                moreinfo = 'The reactant entity number '// trim(int2str(j)) // ' (species ' // &
                           trim(surfspecsnames(kspec)) // ') has at least one dentate which has not ' // &
                           'been assigned a site in the lattice state of cluster ' // trim(int2str(i)) // &
                           ' (' // trim(clusternames(i)) // ').'
                call error(8026)
            endif
        enddo
    enddo

    do j = 1,clusterspecs(i,0)
        kspec = clusterspecs(i,j)
        
        ! Populate species' participation arrays 
        nspecclusterparticip(kspec) = nspecclusterparticip(kspec) + 1
        specclusterparticip(kspec,nspecclusterparticip(kspec),1) = i
        specclusterparticip(kspec,nspecclusterparticip(kspec),2) = j
        
        ! Define as start sites the ones covered by this molecule 
        nstart = 0
        do m = 1,clusternsites(i)
            mmolec = clusterlatticestate(i,m,1)
            if (mmolec == j) then
                nstart = nstart + 1
                startsites(clusterlatticestate(i,m,2)) = m
            endif
        enddo
        
        if (nstart /= surfspecsdent(kspec)) then
            write(iwrite,*) 'DEBUGGING: Error in finding the d-level of pattern ' // trim(int2str(i)) // '.'
            stop
        endif
        
        call find_max_level(clusternsites(i),clusterneigh(i,:,:),dmax,nstart,startsites)
        
        clusterlevels(i,j) = dmax
        if (clustermaxlevel < dmax) clustermaxlevel = dmax
    enddo
    
enddo

! Output information for all the clusters

write(iwrite,'(/,a)') '    Number of clusters: ' // trim(int2str(nclusters))
write(iwrite,'(/,a)') '    Maximum number of sites involved in a cluster: ' // trim(int2str(clusternmxsites))

maxcharsname = 0
do i = 1,nclusters
    if (len_trim(clusternames(i)) > maxcharsname) then
        maxcharsname = len_trim(clusternames(i))
    endif
enddo

! Briefly report the clusters
write(iwrite,'(/,a)') '    Clusters:'
do i = 1,nclusters
    write(iwrite,'(/,6x,a4,1x,a,":"' // int2str(maxcharsname+1-len_trim(clusternames(i))) // 'x,)', &
                    advance='no') trim(int2str(i)) // '.', trim(clusternames(i))
    
    write(iwrite,'(1x,"Mult = ",a6)',advance='no') int2str(clustergraphmultipl(i))
    write(iwrite,'(1x,"ECI = ",a25)',advance='no') dbl2str(clusterenrg(i))
    write(iwrite,'(1x,"Entities: ",' // int2str(clusternmxsites) // '(1x,a),)', &
                    advance='no') (trim(surfspecsnames(clusterspecs(i,j))),j=1,clusterspecs(i,0))
enddo

write(iwrite,'(a)') ''
write(iwrite,'(/,a)') 'Finished reading energetics input.'

return

end subroutine read_energetics_setup

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine parse_cluster(firstpass,recinput,nwords,irec,io,icluster,iw,clusteropen)

use constants_module, only: maxwords2, lengthrecinp, remchar, imechread,       &
                            cenergfname
use error_module, only: moreinfo, error
use parser_module, only: int2str, striccompare, str_check_expression, str2int, &
                         findinstrv, break_words, getrecord
use simulation_setup_module, only: surfspecsdent, surfspecsnames, nsurfspecs

implicit none

integer iadsorb, identat, io, irec, icluster, ispec, k, m
integer nv(25), tmpv(2), nsites
integer nwords, nwords2, ivariant, j, j1, j2
integer iw(maxwords2), iw2(maxwords2)

character(lengthrecinp) recinput, tmpstr
character(nnam0) basename, varname

logical clusteropen
logical readnsites, readstypes, readgeom, readgraphmult, readlatstate, readneighb
logical readeng, canspecifyvariant, firstpass, validexp

! Parses a single cluster in energetics_input.dat

! Level-1 keywords allowed in this subroutine:
!  -Normal keywords recognized in this scope
!    sites
!    lattice_state
!    neighboring
!    end_cluster
!  -The following keywords direct the program flow to subroutine parse_stypes_geom_eng
!    site_types
!    angles
!    cluster_eng
!  -The following keywords direct the program flow to subroutine parse_variant
!    variant
!  -The following keywords are recognized as invalid in this context and produce and error
!    cluster

  ! Depending on the value of the logical variable firstpass the subroutine can be run in one of 
  ! two "modes": a first pass mode and a full parsing mode. The former just checks for consistency
  ! of input and finds the max bounds needed for the arrays that encode the mechanism. The latter
  ! populates these arrays. The user MUST run a first pass of the mechanism before the full parsing.

  readnsites = .false.
  readstypes = .false.
  readlatstate = .false.
  readneighb = .false.
  readeng = .false.
  readgeom = .false.
  readgraphmult = .false.
  canspecifyvariant = .true.

  icluster = icluster + 1
  ivariant = 0
  if (nwords == 2) then
      basename = recinput(iw(3):iw(4))
  elseif (nwords == 1) then
      basename = 'Cluster_' // trim(int2str(icluster))
  else
      call error(8002)
  endif

  if (.not.firstpass) then
          clusternames(icluster) = basename
  endif

  clusteropen = .true.
  do while (io >= 0)

      call getrecord(imechread,recinput,irec,io)
      
      call break_words(recinput,' ',remchar,nwords,iw)

      if (nwords == 0) then
          cycle
      endif
      
      if (striccompare(recinput(iw(1):iw(2)),'sites')) then

          ! If a first pass is requested ...
          if (firstpass) then

              if (readnsites) then
                  moreinfo = 'Keyword ' // recinput(iw(1):iw(2))                 &
                             // ' has already been parsed; yet, it '             &
                             // 'appears again in line ' // trim(int2str(irec))  &
                             // ' in ' // trim(cenergfname) // '.'
                  call error(100)
              endif
          
              if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
                  nsites = str2int(recinput(iw(3):iw(4)))
                  clusternmxsites = max(nsites,clusternmxsites)
              else
                  moreinfo = 'Invalid expression "'                              &
                             // recinput(iw(1):iw(2*nwords)) // '" in line '     &
                             // trim(int2str(irec)) // '.'
                  call error(8005)
              endif

              readnsites = .true.
              cycle
              
          endif

          ! Full parsing is requested:

          clusternsites(icluster) = str2int(recinput(iw(3):iw(4)))
          
          readnsites = .true.
          
      elseif (striccompare(recinput(iw(1):iw(2)),'lattice_state')) then
          
          ! For a first pass ...
          if (firstpass) then

              if (readlatstate) then
                  moreinfo = 'Keyword ' // recinput(iw(1):iw(2))                 &
                             // ' has already been parsed; yet, it '             &
                             // 'appears again in line ' // trim(int2str(irec))  &
                             // ' in ' // trim(cenergfname) // '.'
                  call error(100)
              endif
              
              ! ... we check if this keyword is preceded by the sites keyword and then
              ! skip this part of the input ...
              if (.not.readnsites) then
                  moreinfo = 'Occurred while parsing line '                      &
                             // trim(int2str(irec)) // '.'
                  call error(8006)
              endif
              
              j = 0
              do while (j < nsites .and. io >= 0)
              
                  call getrecord(imechread,recinput,irec,io)
                  call break_words(recinput,' ',remchar,nwords,iw)
                  if (nwords == 0) cycle
                  
                  j = j + 1
                  
                  validexp = .false.
                  if (nwords == 3) then
                      if (str_check_expression(recinput,nwords,iw,'I4/AA/I4',.true.)) then
                          if (.not.(striccompare(recinput(iw(3):iw(4)),'&'))) then
                              validexp = .true.
                          endif
                      else
                          if (striccompare(recinput(iw(1):iw(2)),'&') & 
                              .and. striccompare(recinput(iw(3):iw(4)),'&') & 
                              .and. striccompare(recinput(iw(5):iw(6)),'&') ) then
                              validexp = .true.
                          endif
                      endif
                  endif
                  
                  if (.not.validexp) then
                      moreinfo = 'Invalid expression "'                          &
                                 // recinput(iw(1):iw(2*nwords))                 &
                                 // '" in line ' // trim(int2str(irec)) // '.'
                      call error(8007)
                  endif
                  
              enddo
              
              readlatstate = .true.
              cycle
              
          endif
                  
          ! Full parsing requested...
          j = 0
          do while (j < clusternsites(icluster) .and. io >= 0)
          
              call getrecord(imechread,recinput,irec,io)
              call break_words(recinput,' ',remchar,nwords,iw)
              if (nwords == 0) cycle
                          
              j = j + 1
              
              ! If the "&" symbol appears, this means unpecified site state; take care of it and continue
              if (striccompare(recinput(iw(1):iw(2)),'&')) then
                  clusterlatticestate(icluster,j,1) = -1
                  clusterlatticestate(icluster,j,2) = -1
                  cycle
              endif

              iadsorb = str2int(recinput(iw(1):iw(2)))
              ispec = findinstrv(surfspecsnames,0,0,nsurfspecs,recinput(iw(3):iw(4)))
            if (ispec < 0) then
                moreinfo = 'Encountered unknown surface species "' // recinput(iw(3):iw(4)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(8004)    
            endif
            identat = str2int(recinput(iw(5):iw(6)))

            if (identat < 0 .or. identat > surfspecsdent(ispec)) then
                moreinfo = 'Dentate number ' // trim(int2str(identat)) // ' for species ' // recinput(iw(3):iw(4)) // &
                           ' is out of bounds 1:' // trim(int2str(surfspecsdent(ispec))) // ' in line ' // trim(int2str(irec)) // '.'
                call error(8028)
            endif

            clusterlatticestate(icluster,j,1) = iadsorb
            clusterlatticestate(icluster,j,2) = identat
            
            if (clusterspecs(icluster,iadsorb) >= 0) then
            
                ! Check if this adsorbate/entity has previously been encountered
                if (clusterspecs(icluster,iadsorb) /= ispec) then
                    moreinfo = 'Adsorbate/entity ' // trim(int2str(iadsorb)) // ' has previously been declared as ' // &
                               trim(surfspecsnames(clusterspecs(icluster,iadsorb))) // & 
                               ' species but now appears as ' // trim(surfspecsnames(ispec)) // &
                               ' in line ' // trim(int2str(irec)) // '.'
                    call error(8008)
                endif
            else
                
                clusterspecs(icluster,0) = clusterspecs(icluster,0) + 1
                clusterspecs(icluster,iadsorb) = ispec                            
                
            endif
            
        enddo
        
        readlatstate = .true.
                
    elseif (striccompare(recinput(iw(1):iw(2)),'neighboring')) then
        
        ! For a first pass...
        if (firstpass) then
        
            nv = 0 ! remark: nv is a vector
            
            do j = 1,nwords-1
                
                tmpstr = recinput(iw(2*j+1):iw(2*j+2))
                call break_words(trim(tmpstr),'-',remchar,nwords2,iw2)
                
                if (.not.str_check_expression(tmpstr,nwords2,iw2,'2I4',.true.)) then
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(8009)
                endif                    
                
                j1 = str2int(tmpstr(iw2(1):iw2(2)))
                j2 = str2int(tmpstr(iw2(3):iw2(4)))
                
                if (j1 == j2) then
                    moreinfo = 'Encountered "' // recinput(iw2(1):iw2(2)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(8010)
                endif
                
                tmpv = (/j1,j2/)
                do k = 1,2                                    
                    if (tmpv(k) > nsites .or. tmpv(k) < 1) then
                        moreinfo = 'Site number ' // trim(int2str(tmpv(k))) // &
                                   ' referenced in neighboring input in line ' // trim(int2str(irec)) // &
                                   ' is out of range 1:' // trim(int2str(nsites)) // &
                                   ' for this cluster.'
                        call error(8009)
                    endif
                    
                enddo
                
                nv(j1) = nv(j1) + 1
                nv(j2) = nv(j2) + 1
                clusternmxcoord = max(clusternmxcoord,nv(j1),nv(j2))
                
            enddo
            
            readneighb = .true.
            cycle
            
        endif

        ! Full parsing requested
        do j = 1,nwords-1
            
            tmpstr = recinput(iw(2*j+1):iw(2*j+2))
            call break_words(trim(tmpstr),'-',remchar,nwords2,iw2)
            
            j1 = str2int(tmpstr(iw2(1):iw2(2)))
            j2 = str2int(tmpstr(iw2(3):iw2(4)))
            
            clusterneigh(icluster,j1,0) = clusterneigh(icluster,j1,0) + 1
            clusterneigh(icluster,j1,clusterneigh(icluster,j1,0)) = j2

            clusterneigh(icluster,j2,0) = clusterneigh(icluster,j2,0) + 1
            clusterneigh(icluster,j2,clusterneigh(icluster,j2,0)) = j1
            
            ! Check if the j1-j2 pair has already been encountered in the neighboring list
            do k = 1,clusterneigh(icluster,j1,0)
                do m = k+1,clusterneigh(icluster,j1,0)
                    if (clusterneigh(icluster,j1,k) == clusterneigh(icluster,j1,m)) then
                        moreinfo = 'Pair ' // trim(int2str(j1)) // '-' // trim(int2str(j2)) // &
                                   ' is referenced twice in neighboring input in line ' // trim(int2str(irec)) // '.'
                        call error(8011)
                    endif
                enddo
            enddo
                        
        enddo
                
        ! Check whether a site has been left unconnected
        do j = 1,clusternsites(icluster)
            if (clusterneigh(icluster,j,0) == 0) then
                moreinfo = 'Site number ' // trim(int2str(j)) // &
                           ' is not referenced in neighboring input in line ' // trim(int2str(irec)) // '.'
                call error(8012)
            endif
        enddo
        
        readneighb = .true.
        
    elseif (striccompare(recinput(iw(1):iw(2)),'site_types') .or. &
            striccompare(recinput(iw(1):iw(2)),'graph_multiplicity') .or. &
            striccompare(recinput(iw(1):iw(2)),'cluster_eng') .or. &
            striccompare(recinput(iw(1):iw(2)),'angles') .or. &
            striccompare(recinput(iw(1):iw(2)),'no_mirror_images') .or. &
            striccompare(recinput(iw(1):iw(2)),'absl_orientation')) then
        
        if (ivariant > 0) then
            moreinfo = 'Invalid specification in line ' // trim(int2str(irec)) // '.'
            call error(8054)
        endif

        canspecifyvariant = .false.
                
        call parse_stypes_geom_eng(firstpass,recinput,irec,io,nwords,iw,icluster,nsites, &
            readnsites,readneighb,readstypes,readgeom,readgraphmult,readeng)
 
    elseif (striccompare(recinput(iw(1):iw(2)),'variant')) then
        
        ivariant = ivariant + 1

        readstypes = .false.
        readeng = .false.
        readgeom = .false.
        readgraphmult = .false.

        if (firstpass) then

            if (.not.canspecifyvariant) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call error(8027)
            endif
                    
            if (ivariant == 1) then
                ! If this is the first variant found, check whether the information
                ! about sites, initial and final state and neighboring structure has been read
                call check_level1_cluster_completeness(irec,icluster,basename,nsites, &
                        readnsites,readlatstate,readneighb)
            endif
        
            readeng = .false.
        
        endif        
        
        call parse_variant(firstpass,ivariant,varname,basename,recinput,irec,io,nwords,iw,icluster,nsites, &
               readnsites,readneighb,readstypes,readeng,readgeom,readgraphmult)
        
        if (firstpass) then
            call check_level2_cluster_completeness(irec,icluster,basename,ivariant,varname,readeng)
        endif
        
    elseif ( striccompare(recinput(iw(1):iw(2)),'end_cluster') .and. clusteropen ) then
        
        clusteropen = .false.
        
        if (firstpass) then

            if (ivariant > 0) then
                nclusters = nclusters + ivariant
            else
                nclusters = nclusters + 1
            endif
            
        endif
        
        exit
    
    elseif (striccompare(recinput(iw(1):iw(2)),'cluster')) then

        moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
        call error(8023)

    else
        moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
        call error(8100)
    endif
        
enddo ! loop for reading info for one cluster

if (ivariant == 0 .and. firstpass) then   
    call check_level1_cluster_completeness(irec,icluster,basename,nsites, &
                        readnsites,readlatstate,readneighb)
    call check_level2_cluster_completeness(irec,icluster,basename,ivariant,varname,readeng)
endif

return

end subroutine parse_cluster

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine parse_stypes_geom_eng(firstpass, recinput, irec, io, nwords, iw,    &
                                 icluster, nsites, readnsites, readneighb,     &
                                 readstypes, readgeom, readgraphmult, readeng)

use constants_module, only: maxwords2, lengthrecinp, remchar, cenergfname
use parser_module, only: str_is_integer, int2str, striccompare, str2int,       &
                         str_check_expression, dbl2str, str2dbl, findinstrv,   &
                         break_words
use error_module, only: error, moreinfo
use lattice_setup_module, only: nsitetypes, sitetypenames
implicit none

integer io, irec, icluster, isite, istype, nsites
integer j, k, m, p
integer nwords
integer iw(maxwords2)
integer nwords2
integer iw2(maxwords2)
integer nwords3
integer iw3(maxwords2)

integer sites123(3)

character(lengthrecinp) recinput, tmpstr1, tmpstr2

logical firstpass, readnsites, readneighb, readeng, readstypes, readgeom,      &
        readgraphmult, neighfound !, readdeltah, readperatio
!, isreversible, readpreexp, 
logical tmplog1, tmplog2

! Level-1 keywords allowed in this subroutine:
!  -Normal keywords recognized in this scope
!    site_types
!    cluster_eng

if (striccompare(recinput(iw(1):iw(2)),'site_types')) then
    
    if (firstpass) then
        
        if (.not. readnsites) then
            moreinfo = 'Occurred while parsing line '                          &
                       // trim(int2str(irec)) // '.'
            call error(8031)
        endif
        
        if (nwords-1 == nsites) then
            readstypes = .true.
            return
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords))  &
                       // '" in line ' // trim(int2str(irec)) // '.'
            call error(8013)
        endif        
        
    endif
    
    do j = 1,nwords-1
        
        if (str_is_integer(recinput(iw(2*j+1):iw(2*j+2)),4)) then
            istype = str2int(recinput(iw(2*j+1):iw(2*j+2)))
            if (istype > nsitetypes) then
                moreinfo = 'Encountered site type ' // trim(int2str(istype)) // ' for site ' // trim(int2str(j)) // &
                           ' in line ' // trim(int2str(irec)) // '.'
                call error(8014)
            endif
        else
            istype = findinstrv(sitetypenames,1,1,nsitetypes,recinput(iw(2*j+1):iw(2*j+2)))                                
            if (istype < 0) then
                moreinfo = 'Encountered site type name "' // recinput(iw(2*j+1):iw(2*j+2)) // '" for site ' // trim(int2str(j)) // &
                           ' in line ' // trim(int2str(irec)) // '.'
                call error(8014)
            endif
        endif
        
        clusterstype(icluster,j) = istype
                                            
    enddo
    
    readstypes = .true.
    
    return
    
elseif (striccompare(recinput(iw(1):iw(2)),'graph_multiplicity')) then

    ! If a first pass is requested ...
    if (firstpass) then

        if (readgraphmult) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(cenergfname) // '.'
            call error(100)
        endif
    
        if (.not. str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(80)
        endif

        readgraphmult = .true.
        return
        
    endif

    ! Full parsing is requested:

    clustergraphmultipl(icluster) = str2int(recinput(iw(3):iw(4)))
    
    readgraphmult = .true.

elseif (striccompare(recinput(iw(1):iw(2)),'angles')) then
    
    if (firstpass) then

        if (.not. readnsites) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(8032)
        endif

        if (.not. readneighb) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(8033)
        endif
        
        do j = 2,nwords ! Check validity of each expression; it should be site1-site2-site3:angle_value
            
            tmpstr1 = recinput(iw(2*j-1):iw(2*j))
            call break_words(trim(tmpstr1),':',remchar,nwords2,iw2)
            
            tmpstr2 = tmpstr1(iw2(1):iw2(2))
            call break_words(trim(tmpstr2),'-',remchar,nwords3,iw3)
            
            tmplog1 = str_check_expression(tmpstr1,nwords2,iw2,'AA/R8',.true.) .or. &
                      str_check_expression(tmpstr1,nwords2,iw2,'AA/I4',.true.)
            tmplog2 = str_check_expression(tmpstr2,nwords3,iw3,'I4/I4/I4',.true.)
            
            if (tmplog1 .and. tmplog2) then
                clusternmxangles = clusternmxangles + 1
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(8034)
            endif

        enddo
        
        readgeom = .true.
        return
        
    endif
            
    do j = 2,nwords ! Check validity of each expression
        
        tmpstr1 = recinput(iw(2*j-1):iw(2*j))
        call break_words(trim(tmpstr1),':',remchar,nwords2,iw2)
        
        tmpstr2 = tmpstr1(iw2(1):iw2(2))
        call break_words(trim(tmpstr2),'-',remchar,nwords3,iw3)
        
        sites123(1) = str2int(tmpstr2(iw3(1):iw3(2)))
        sites123(2) = str2int(tmpstr2(iw3(3):iw3(4)))
        sites123(3) = str2int(tmpstr2(iw3(5):iw3(6)))
            
        do k = 1,3
            if (sites123(k) < 1 .or. sites123(k) > clusternsites(icluster)) then
                moreinfo = 'Site number ' // trim(int2str(sites123(k))) // &
                           ' referenced in angle input in line ' // trim(int2str(irec)) // &
                           ' is out of range 1:' // trim(int2str(nsites)) // &
                           ' for this cluster.'
                call error(8035)
            endif
        enddo
        
        if (sites123(1) == sites123(2) .or. sites123(1) == sites123(3) .or. sites123(2) == sites123(3)) then
            moreinfo = 'Invalid repetition of sites in angle specification for a cluster in line ' // trim(int2str(irec))
            call error(8035)
        endif
        
        ! See if the neighboring implied by the angle specification is compatible with that of the cluster
        do k = 1,2
            m = k+1
            neighfound = .false.
            do p = 1,clusterneigh(icluster,sites123(k),0)
                if (sites123(m) == clusterneigh(icluster,sites123(k),p)) then
                    neighfound = .true.
                    exit
                endif
            enddo
            
            if (.not.neighfound) then
                moreinfo = 'Site number ' // trim(int2str(sites123(k))) // &
                           ' neighbors with site number ' // trim(int2str(sites123(m))) // & 
                           ' in angle input in line ' // trim(int2str(irec)) // '.' // &
                           ' This is invalid for this cluster.'
                call error(8036)
            endif
                
        enddo
        
        clusternanglessequenc(icluster,0) = clusternanglessequenc(icluster,0) + 1
        m = clusternanglessequenc(icluster,0)
        clusternanglessequenc(icluster,3*m-2) = sites123(1)
        clusternanglessequenc(icluster,3*m-1) = sites123(2)
        clusternanglessequenc(icluster,3*m) = sites123(3)
        
        clusterangles(icluster,m) = str2dbl(tmpstr1(iw2(3):iw2(4)))
    
    enddo
    
elseif (striccompare(recinput(iw(1):iw(2)),'cluster_eng')) then

    if (firstpass) then

        if (str_check_expression(recinput,nwords,iw,'AA/R8',.true.)) then
            readeng = .true.
            return
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(8018)
        endif
        
    endif
        
    clusterenrg(icluster) = str2dbl(recinput(iw(3):iw(4)))

    readeng = .true.
    
    return

elseif (striccompare(recinput(iw(1):iw(2)),'no_mirror_images')) then

    if (firstpass) then
        return
    endif
    
    clusternomirrorimgs(icluster) = .true.

elseif (striccompare(recinput(iw(1):iw(2)),'absl_orientation')) then
    
    if (firstpass) then

        if (.not. readnsites) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(8050)
        endif

        if (.not. readneighb) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(8051)
        endif
        
        ! Check validity of the expression; it should be site1-site2:angle_value
            
        tmpstr1 = recinput(iw(3):iw(4))
        call break_words(trim(tmpstr1),':',remchar,nwords2,iw2)
            
        tmpstr2 = tmpstr1(iw2(1):iw2(2))
        call break_words(trim(tmpstr2),'-',remchar,nwords3,iw3)
            
        tmplog1 = str_check_expression(tmpstr1,nwords2,iw2,'AA/R8',.true.) .or. &
                    str_check_expression(tmpstr1,nwords2,iw2,'AA/I4',.true.)
        tmplog2 = str_check_expression(tmpstr2,nwords3,iw3,'I4/I4',.true.)
            
        if (.not.(tmplog1 .and. tmplog2)) then
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(8052)
        endif

        return
        
    endif
            
    ! Check validity of the expression
        
    tmpstr1 = recinput(iw(3):iw(4))
    call break_words(trim(tmpstr1),':',remchar,nwords2,iw2)
        
    tmpstr2 = tmpstr1(iw2(1):iw2(2))
    call break_words(trim(tmpstr2),'-',remchar,nwords3,iw3)
        
    sites123(1) = str2int(tmpstr2(iw3(1):iw3(2)))
    sites123(2) = str2int(tmpstr2(iw3(3):iw3(4)))
            
    do k = 1,2
        if (sites123(k) < 1 .or. sites123(k) > clusternsites(icluster)) then
            moreinfo = 'Site number ' // trim(int2str(sites123(k))) // &
                        ' referenced in absolute orientation input in line ' // trim(int2str(irec)) // &
                        ' is out of range 1:' // trim(int2str(nsites)) // &
                        ' for this cluster.'
            call error(8053)
        endif
    enddo
        
    if (sites123(1) == sites123(2)) then
        moreinfo = 'Invalid repetition of sites in absolute orientation specification for a cluster in line ' // trim(int2str(irec))
        call error(8053)
    endif
        
    ! See if the neighboring implied by the angle specification is compatible with that of the cluster
    neighfound = .false.
    do p = 1,clusterneigh(icluster,sites123(1),0)
        if (sites123(2) == clusterneigh(icluster,sites123(1),p)) then
            neighfound = .true.
            exit
        endif
    enddo
            
    if (.not.neighfound) then
        moreinfo = 'Site number ' // trim(int2str(sites123(k))) // &
                    ' neighbors with site number ' // trim(int2str(sites123(m))) // & 
                    ' in absolute orientation input in line ' // trim(int2str(irec)) // '.' // &
                    ' This is invalid for this cluster.'
        call error(8053)
    endif
                        
    clusterabslorientat(icluster) = .true.
    clusterorientationedge(icluster,1) = sites123(1)
    clusterorientationedge(icluster,2) = sites123(2)
        
    clusterorientationangles(icluster) = str2dbl(tmpstr1(iw2(3):iw2(4)))
        
else
    
    moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
    call error(8100)
    
endif

end subroutine parse_stypes_geom_eng

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine parse_variant(firstpass, ivariant, varname, basename, recinput,     &
                         irec, io, nwords, iw, icluster, nsites, readnsites,   &
                         readneighb, readstypes, readeng, readgeom,            &
                         readgraphmult)

use constants_module, only: maxwords2, lengthrecinp, remchar, imechread
use parser_module, only: int2str, striccompare, break_words, getrecord
use error_module, only: moreinfo, error


implicit none

integer i, io, irec, icluster, nsites
integer ivariant, j, k, k1, k2, dk, nwords
integer iw(maxwords2)

character(lengthrecinp) recinput
character(nnam0) basename, varname

logical firstpass, readeng, readstypes, readgeom, readgraphmult, readnsites, readneighb

if (.not. firstpass) then

    if (ivariant > 1) then

        k1 = icluster
        k2 = icluster
        dk = 1
        icluster = icluster + 1   

        do k = k1,k2
               
            clusternsites(k+dk) = clusternsites(k)
            
            do i = 1,clusternsites(k)
                do j = 0,clusterneigh(k,i,0)
                    clusterneigh(k+dk,i,j) = clusterneigh(k,i,j)
                enddo
            enddo
            
            do i = 1,clusternsites(k)
                do j = 1,2
                    clusterlatticestate(k+dk,i,j) = clusterlatticestate(k,i,j)
                enddo
            enddo
            
            do i = 0,clusternsites(k)
                clusterspecs(k+dk,i) = clusterspecs(k,i)
            enddo
            
        enddo

    endif

endif

if (nwords > 1) then
    varname = recinput(iw(3):iw(4))
else
    varname = 'var' // int2str(ivariant)
endif

if (.not. firstpass) then
    clusternames(icluster) = trim(basename) // '_' // varname
endif

do while (io >= 0)

    call getrecord(imechread,recinput,irec,io)
    
    call break_words(recinput,' ',remchar,nwords,iw)

    if (nwords > 0) then

        if (striccompare(recinput(iw(1):iw(2)),'site_types') .or. &
            striccompare(recinput(iw(1):iw(2)),'graph_multiplicity') .or. &
            striccompare(recinput(iw(1):iw(2)),'cluster_eng') .or. &
            striccompare(recinput(iw(1):iw(2)),'angles') .or. &
            striccompare(recinput(iw(1):iw(2)),'no_mirror_images') .or. &
            striccompare(recinput(iw(1):iw(2)),'absl_orientation')) then

            call parse_stypes_geom_eng(firstpass,recinput,irec,io,nwords,iw,icluster,nsites, &
                    readnsites,readneighb,readstypes,readgeom,readgraphmult,readeng)
        
        elseif (striccompare(recinput(iw(1):iw(2)),'end_variant')) then
            
            return
            
        else
            
            moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
            call error(8100)
            
        endif
                
    endif
    
enddo

return

end subroutine parse_variant

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_level1_cluster_completeness(irec,icluster,clustername,nsites, &
            readnsites,readlatstate,readneighb)

use parser_module, only: int2str
use error_module, only: moreinfo, error

implicit none

integer irec, icluster, nsites

character(nnam0) clustername

logical clusternotdone, readnsites, readlatstate, readneighb

! Check for completeness of information:
clusternotdone = .false.
moreinfo = 'Missing information for step ' // trim(int2str(icluster)) // ' (' // &
            trim(clustername) // '):'
if (.not.readnsites) then
    clusternotdone = .true.
    moreinfo = trim(moreinfo) // ' number of sites,'
endif
if (.not.readlatstate) then
    clusternotdone = .true.
    moreinfo = trim(moreinfo) // ' initial state,'
endif
if ((.not.readneighb) .and. nsites > 1) then
    clusternotdone = .true.
    moreinfo = trim(moreinfo) // ' neighboring structure,'
endif
moreinfo = trim(moreinfo) // ' has(have) not been read. ' // &
           'Occurred while parsing line ' // trim(int2str(irec)) // '.'

if (clusternotdone) then
    call error(8024)
endif

return

end subroutine check_level1_cluster_completeness

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_level2_cluster_completeness(irec,icluster,clustername,ivariant,varname,readeng)

use parser_module, only: int2str
use error_module, only: moreinfo, error
implicit none

integer irec, icluster, ivariant

logical clusternotdone, readeng
character(nnam0) varname, clustername

! Check for completeness of information:
clusternotdone = .false.
moreinfo = 'Missing information for cluster ' // trim(int2str(icluster)) // ' (' // &
            trim(clustername) // ')'

if (ivariant > 0) then
    moreinfo = trim(moreinfo) // ' variant ' // trim(int2str(ivariant)) // ' (' // &
            trim(varname) // ')'
endif
moreinfo = trim(moreinfo) // ':'

if (.not.readeng) then
    clusternotdone = .true.
    moreinfo = trim(moreinfo) // ' cluster energy,'
endif

moreinfo = trim(moreinfo) // ' has(have) not been read. ' // &
           'Occurred while parsing line ' // trim(int2str(irec)) // '.'

if (clusternotdone) then
    call error(8024)
endif

end subroutine check_level2_cluster_completeness

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function factorial(n)

integer i, n, p

p = 1
do i = 1, n
	p = p * i
enddo

factorial = p
return

end function factorial

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module energetics_setup_module
