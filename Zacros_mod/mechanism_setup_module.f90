! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module mechanism_setup_module

use constants_module
use error_module
use parser_module
use graph_functions_module

use simulation_setup_module

implicit none

integer nelemsteps, nSAparams, nreversible, elemstepnmxsites, elemstepnmxcoord, elemstepnmxgas, elemstepnmxangles

integer, allocatable :: reverselemstep(:) 
integer, allocatable :: elemstepnsites(:) ! Number of species involved in step i
integer, allocatable :: elemstepreactnts(:,:)
integer, allocatable :: elemstepproducts(:,:)
integer, allocatable :: elemstepneigh(:,:,:)
integer, allocatable :: elemstepstype(:,:)
integer, allocatable :: elemsteplatticestatein(:,:,:)
integer, allocatable :: elemsteplatticestatefn(:,:,:)
integer, allocatable :: elemstepadsorbposin(:,:,:)
integer, allocatable :: elemstepadsorbposfn(:,:,:)
integer, allocatable :: elemstepgases(:,:)
integer, allocatable :: elemstepnanglessequenc(:,:)
integer, allocatable :: elemsteporientationedge(:,:)

integer, allocatable :: nspecparticip(:)
integer, allocatable :: specparticip(:,:,:)

integer, allocatable :: elemsteplevels(:,:)

integer(8), allocatable :: elemstep_noccur(:)

real(8), allocatable :: preexp(:,:)
real(8), allocatable :: acteng(:)
real(8), allocatable :: omega(:)
real(8), allocatable :: elemstep_avgtime(:)
real(8), allocatable :: elemstepangles(:,:)
real(8), allocatable :: elemsteporientationangles(:)
real(8), allocatable :: propCountvec(:)

logical, allocatable :: preexpisconst(:)
logical, allocatable :: elemstepnomirrorimgs(:)
logical, allocatable :: elemstepabslorientat(:)

character(nnam0), allocatable :: elemstepnames(:)

contains

  subroutine cleanup_msm
    ! cleans up module globals

    if(allocated(reverselemstep))         deallocate(reverselemstep)
    if(allocated(elemstepnsites))         deallocate(elemstepnsites)
    if(allocated(elemstepreactnts))       deallocate(elemstepreactnts)
    if(allocated(elemstepproducts))       deallocate(elemstepproducts)
    if(allocated(elemstepneigh))          deallocate(elemstepneigh)
    if(allocated(elemstepstype))          deallocate(elemstepstype)
    if(allocated(elemsteplatticestatein)) deallocate(elemsteplatticestatein)
    if(allocated(elemsteplatticestatefn)) deallocate(elemsteplatticestatefn)
    if(allocated(elemstepadsorbposin))    deallocate(elemstepadsorbposin)
    if(allocated(elemstepadsorbposfn))    deallocate(elemstepadsorbposfn)
    if(allocated(elemstepgases))          deallocate(elemstepgases)
    if(allocated(elemstepnanglessequenc)) deallocate(elemstepnanglessequenc)
    if(allocated(elemstepnomirrorimgs))   deallocate(elemstepnomirrorimgs)
    if(allocated(nspecparticip))          deallocate(nspecparticip)
    if(allocated(specparticip))           deallocate(specparticip)
    if(allocated(elemsteplevels))         deallocate(elemsteplevels)
    if(allocated(elemstep_noccur))        deallocate(elemstep_noccur)
    if(allocated(preexp))                 deallocate(preexp)
    if(allocated(acteng))                 deallocate(acteng)
    if(allocated(omega))                  deallocate(omega)
    if(allocated(elemstep_avgtime))       deallocate(elemstep_avgtime)
    if(allocated(elemstepangles))         deallocate(elemstepangles)
    if(allocated(preexpisconst))          deallocate(preexpisconst)
    if(allocated(elemstepnames))          deallocate(elemstepnames)

  end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine read_mechanism_setup()

implicit none

integer i, iadsorb, identat, io, irec, istep, istype, j, m, mmolec
integer k, nwords, nstart, kspec, kdent, dmax, q, setvbuf3f
integer maxcharsname, iw(maxwords2)
integer, allocatable :: startsites(:)

real(8) preexpfac, ktemp

character(lengthrecinp) recinput

logical genoutop, stepopen, mechopen, defaultnaming, mechfex

! Parses the simulation setup file simulation_input.dat

! Level-1 keywords allowed in this subroutine:
!    mechanism
!    end_mechanism
!    step
!    reversible_step

! The program first performs a parsing of the entire input to check whether the 
! blocks are opened and closed properly and collect information that will be used 
! for determining the bounds of the arrays used for mechanism specification.

! Check if the lattice specification file exists
mechfex = .false.
inquire(file=cmechfname,exist=mechfex)
if (.not.mechfex) then
    moreinfo = 'Mechanism specification file ' // trim(cmechfname) // ' does not exist.'
    call error(1)
endif

inquire(iwrite,opened=genoutop)
if (.not.genoutop) then
    open(unit=iwrite,file=trim(cgenoutfname),status='unknown',position='append')
    q = setvbuf3f(iwrite,1,0) ! set line-buffered behavior
endif

write(iwrite,'(/,a)') 'Mechanism setup:'
write(iwrite,'(a)')   '~~~~~~~~~~~~~~~~'

open(unit=imechread,file=trim(cmechfname),status='old')

! Perform an initial pass to check for proper syntax and opening and closing
! of blocks and find the following:
!    number of elementary steps: nelemsteps
!    the max number of sites involved in a step: elemstepnmxsites
!    the max coordination number of sites involved in a step: elemstepnmxcoord
!    the max number of gas species participating in a step: elemstepnmxgas

istep = 0
mechopen = .false.
stepopen = .false.
nelemsteps = 0
nreversible = 0
elemstepnmxsites = 0
elemstepnmxcoord = 0
elemstepnmxgas = 0
elemstepnmxangles = 0

irec = 0
io = 0

call getrecord(imechread,recinput,irec,io)

do while (io >= 0)
    
    call break_words(recinput,' ',remchar,nwords,iw)
    
    if (nwords == 0) then
        call getrecord(imechread,recinput,irec,io)
        cycle
    endif        
        
    if (striccompare(recinput(iw(1):iw(2)),'mechanism')) then
    
        if (mechopen) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(csimfname) // '.'
            call error(100)
        endif
    
        mechopen = .true.

    elseif (striccompare(recinput(iw(1):iw(2)),'end_mechanism')) then
    
        mechopen = .false.
        exit

    elseif ( striccompare(recinput(iw(1):iw(2)),'step') .or. &
        striccompare(recinput(iw(1):iw(2)),'reversible_step') ) then
        
        if (.not. mechopen) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(3001)
        endif
        
        call parse_elemstep(.true.,recinput,nwords,irec,io,istep,iw,stepopen)            
                        
    else
        
        moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
        call error(3100)

    endif

    call getrecord(imechread,recinput,irec,io)

enddo

if (stepopen) then
    moreinfo = 'Elementary step information not completed in line ' // trim(int2str(irec)) // '.'
    call error(3024)
endif
if (mechopen) then
    call warning(3100)
endif
if (nelemsteps == 0) then
    moreinfo = 'File does not appear to contain any elementary step information.'
    call error(3001)
endif

nSAparams = nelemsteps				! This is how many parameters we are doing sensitivity analysis on

! elemstepnames(i1) gives the name of elementary step i1
allocate(elemstepnames(nelemsteps))

! elemstep_noccur and elemstep_avgtime contain information on the elementary step stetistics
allocate(elemstep_noccur(0:nelemsteps))
allocate(elemstep_avgtime(0:nelemsteps))

allocate(propCountvec(1:nelemsteps))

! reverselemstep(i1) gives the reverse step of elementary step i1. If step i1 is
! irreversible reverselemstep(i1) == 0
allocate(reverselemstep(nelemsteps))

! elemstepnsites(i1) gives the number of sites involved in elementary step i1
allocate(elemstepnsites(nelemsteps))

! elemstepreactnts(i1,i2) gives the i2^th molecule species number involved
!   in elementary step i1. If i2 = 0 it gives the total number of molecules.
!   Similarly for elemstepproducts.
allocate(elemstepreactnts(nelemsteps,0:elemstepnmxsites))
allocate(elemstepproducts(nelemsteps,0:elemstepnmxsites))

! elemstepneigh(i1,i2,i3) gives the i3^th neighbor of the i2^th site
!   involved in reaction i1. If i3 = 0 it gives the total number of
!   neighbors for site i2. Note that this structure is the subgraph (reaction pattern)
!   that will be searched for in the lattice, thereby identifying a lattice process.
allocate(elemstepneigh(nelemsteps,elemstepnmxsites,0:elemstepnmxcoord))
allocate(elemstepstype(nelemsteps,elemstepnmxsites))

! elemsteplatticestatein(i1,i2,i3) gives the lattice state (coverage) for the 
!   elementary step i1. For i3 = 1 it gives the entity number occupying site i2, 
!   and for i3 = 2 it gives the dentate number
allocate(elemsteplatticestatein(nelemsteps,elemstepnmxsites,2))
allocate(elemsteplatticestatefn(nelemsteps,elemstepnmxsites,2))

! elemstepadsorbposin(i1,i2,i3) gives the site occupied by dentate i3 of 
! molecule i2 for the elementary step i1. It is the inverse mapping of 
! elemsteplatticestatein
allocate(elemstepadsorbposin(nelemsteps,elemstepnmxsites,maxdent))
allocate(elemstepadsorbposfn(nelemsteps,elemstepnmxsites,maxdent))

! elemstepgases(i1,i2) gives the affected gas species number for i2 odd
!   and the corresponding stoichiometric coefficient for i2 even. For i2 = 0
!   it gives the total number of gas species affected by elementary step i1
allocate(elemstepgases(nelemsteps,0:2*elemstepnmxgas))

! preexp(i1) and acteng(i1) give preexponentials and activation energies for 
!   elementary step i1
allocate(preexp(nelemsteps,0:7))
allocate(preexpisconst(nelemsteps))
allocate(acteng(nelemsteps))

! omega(i1) gives the proximity factor for elementary step i1
allocate(omega(nelemsteps))

! nspecparticip(i1) gives the number of elementary steps in which species i1 
!   participates. If a species participates in one step two times, this is 
!   doubly counted (or multiply-counted in general)
allocate(nspecparticip(0:nsurfspecs))
! specparticip(i1,i2,i3) for i3 = 1 it gives the i2^th elementary step in 
!   which species i1 participates, and for i3 = 2 it gives the molecule number 
!   with which this species appears in that elementary step
allocate(specparticip(0:nsurfspecs,elemstepnmxsites*nelemsteps,2))

! elemsteplevels(i1,i2) gives the graph level of the pattern corresponding 
!   to the elementary step i1 if we start from the sites of molecule i2
allocate(elemsteplevels(nelemsteps,elemstepnmxsites))

! elemstepnanglessequenc(i1,3*(i2-2):3*i2) gives the site sequence i2 in angle definitions for elementary step i1
! It gives the number of such defined sequences for the second argument = 0
allocate(elemstepnanglessequenc(nelemsteps,0:3*elemstepnmxangles))
! elemstepangles(i1) gives the angle values for the angle definitions for elementary step i1
allocate(elemstepangles(nelemsteps,elemstepnmxangles))

! elemstepnomirrorimgs(i1) is .false. by default allowing mirror images of the elementary step i1 to be detected. A value 
! of .true. prevents mirror images from being detected
allocate(elemstepnomirrorimgs(nelemsteps))

! elemstepabslorientat(i1) is .true. if absolute orientation is required for elementary step i1. Combined with
! the prevention of mirror image detection, this feature allows for only translations of the pattern to
! be detected. By default elemstepabslorientat is .false.
allocate(elemstepabslorientat(nelemsteps))
! elemstepabslorientat(i1,1:2) gives the vertexes of the graph that serve as the "orientation key". The
! code calculates the angle between vector (x2-x1,y2-y1) and (1,0) and checks whether it is equal to that
! specified in elemsteporientationangles(i1) that refers to the elementary step i1.
allocate(elemsteporientationedge(nelemsteps,2))
allocate(elemsteporientationangles(nelemsteps))

! Initialization of variables
defaultnaming = .false.
stepopen = .false.
mechopen = .false.
istep = 0
do i = 0,nelemsteps
    elemstep_noccur(i) = 0_8
    elemstep_avgtime(i) = 0.d0
enddo

do i = 1,nelemsteps
	propCountvec(i) = 0.d0
enddo

do i = 1,nelemsteps
    elemstepnames(i) = 'No_name'
    elemstepnsites(i) = 0
    do j = 1,elemstepnmxsites
        do k = 0,elemstepnmxcoord
            elemstepneigh(i,j,k) = 0
        enddo
    enddo
    elemstepnanglessequenc(i,0) = 0
    elemstepnomirrorimgs(i) = .false.
    elemstepabslorientat(i) = .false.
    elemsteporientationedge(i,1:2) = 0
    elemsteporientationangles(i) = 0.d0
    do j = 1,elemstepnmxangles
        elemstepangles(i,j) = 0.d0
        elemstepnanglessequenc(i,3*j-2:3*j) = 0
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

! The first pass has not found any problems in the syntax. Begin parsing the mechanism

rewind(imechread)
irec = 0
io = 0

do while (io >= 0)

    call getrecord(imechread,recinput,irec,io)
    
    call break_words(recinput,' ',remchar,nwords,iw)
    
    if (nwords == 0) cycle
            
    if (striccompare(recinput(iw(1):iw(2)),'end_mechanism')) then

        exit
        
    elseif (striccompare(recinput(iw(1):iw(2)),'step') .or. &
            striccompare(recinput(iw(1):iw(2)),'reversible_step')) then
    
        call parse_elemstep(.false.,recinput,nwords,irec,io,istep,iw,stepopen)
                
    endif
        
enddo

close(imechread)

! From the mechanism specification data create arrays:
!   nspecparticip
!   specparticip
!   specenerg
!   elemstepadsorbposin
!   elemstepadsorbposfn

allocate(startsites(elemstepnmxsites))

! Construct the species'-participation and pattern-levels arrays
! (referring only to reactants)
do kspec = 0,nsurfspecs
    nspecparticip(kspec) = 0
enddo

do i = 1,nelemsteps

    ! Create the inverse mappings that give the positions of all dentates of the 
    ! species participating in the elementary step
    do j = 1,elemstepnsites(i) ! For every site of the pattern
    
        ! Initial state
        iadsorb = elemsteplatticestatein(i,j,1) ! molecule number
        if (iadsorb <= 0 .or. iadsorb > elemstepreactnts(i,0)) then
            moreinfo = 'The reactant entity number '// trim(int2str(iadsorb)) // ' is out of range 1:' // &
                       trim(int2str(elemstepreactnts(i,0))) // & 
                       ' in the initial state of elementary step ' // trim(int2str(i)) // '.'
            call error(3025)
        endif
        identat = elemsteplatticestatein(i,j,2) ! molecule's dentate number
        if (elemstepadsorbposin(i,iadsorb,identat) .ne. 0) then
            moreinfo = 'Dentate ' // trim(int2str(identat)) // ' of entity number '// trim(int2str(iadsorb)) // &
                       ' has already been assigned site ' // trim(int2str(elemstepadsorbposin(i,iadsorb,identat))) // &
                       ' in the initial state of elementary step ' // trim(int2str(i)) // &
                       ' (' // trim(elemstepnames(i)) // ').'
            call error(3044)        
        endif
        elemstepadsorbposin(i,iadsorb,identat) = j
        
        ! Final state
        iadsorb = elemsteplatticestatefn(i,j,1) ! molecule number
        if (iadsorb <= 0 .or. iadsorb > elemstepproducts(i,0)) then
            moreinfo = 'The product entity number '// trim(int2str(iadsorb)) // ' is out of range 1:' // &
                       trim(int2str(elemstepproducts(i,0))) // & 
                       ' in the final state of elementary step ' // trim(int2str(i)) // '.'
            call error(3025)
        endif
        identat = elemsteplatticestatefn(i,j,2) ! molecule's dentate number
        if (elemstepadsorbposfn(i,iadsorb,identat) .ne. 0) then
            moreinfo = 'Dentate ' // trim(int2str(identat)) // ' of entity number '// trim(int2str(iadsorb)) // &
                       ' has already been assigned site ' // trim(int2str(elemstepadsorbposfn(i,iadsorb,identat))) // &
                       ' in the final state of elementary step ' // trim(int2str(i)) // &
                       ' (' // trim(elemstepnames(i)) // ').'
            call error(3044)        
        endif
        elemstepadsorbposfn(i,iadsorb,identat) = j
        
    enddo
    
    ! Check that all dentates of each species have been assigned a site
    do j = 1,elemstepreactnts(i,0)
        kspec = elemstepreactnts(i,j)
        do kdent = 1,surfspecsdent(kspec)
            if (elemstepadsorbposin(i,j,kdent) == 0) then
                moreinfo = 'The reactant entity number '// trim(int2str(j)) // ' (species ' // &
                           trim(surfspecsnames(kspec)) // ') has at least one dentate which has not ' // &
                           'been assigned a site in the initial state of elementary step ' // trim(int2str(i)) // &
                           ' (' // trim(elemstepnames(i)) // ').'
                call error(3026)
            endif
        enddo
    enddo
    do j = 1,elemstepproducts(i,0)
        kspec = elemstepproducts(i,j)
        do kdent = 1,surfspecsdent(kspec)
            if (elemstepadsorbposfn(i,j,kdent) == 0) then
                moreinfo = 'The product entity number '// trim(int2str(j)) // ' (species ' // &
                           trim(surfspecsnames(kspec)) // ') has at least one dentate which has not ' // &
                           'been assigned a site in the final state of elementary step ' // trim(int2str(i)) // &
                           ' (' // trim(elemstepnames(i)) // ').'
                call error(3026)
            endif
        enddo
    enddo
    
    
    do j = 1,elemstepreactnts(i,0)
        kspec = elemstepreactnts(i,j)
        
        ! Populate species' participation arrays 
        nspecparticip(kspec) = nspecparticip(kspec) + 1
        specparticip(kspec,nspecparticip(kspec),1) = i
        specparticip(kspec,nspecparticip(kspec),2) = j
        
        ! Define as start sites the ones covered by this molecule 
        nstart = 0
        do m = 1,elemstepnsites(i)
            mmolec = elemsteplatticestatein(i,m,1)
            if (mmolec == j) then
                nstart = nstart + 1
                startsites(elemsteplatticestatein(i,m,2)) = m
            endif
        enddo
        
        if (nstart /= surfspecsdent(kspec)) then
            write(iwrite,*) 'DEBUGGING: Error in finding the d-level of pattern ' // trim(int2str(i)) // '.'
            stop
        endif
        
        call find_max_level(elemstepnsites(i),elemstepneigh(i,:,:),dmax,nstart,startsites)
        
        elemsteplevels(i,j) = dmax
    enddo

enddo

! Output information for all the elementary steps

write(iwrite,'(/,a)') '    Number of elementary steps: ' // trim(int2str(nelemsteps))
write(iwrite,'(/,a)') '    Maximum number of sites involved in a step: ' // trim(int2str(elemstepnmxsites))

maxcharsname = 0
do i = 1,nelemsteps
    if (len_trim(elemstepnames(i)) > maxcharsname) then
        maxcharsname = len_trim(elemstepnames(i))
    endif
enddo

! Briefly report the reaction network
write(iwrite,'(/,a)') '    Reaction network:'
do i = 1,nelemsteps
    write(iwrite,'(/,6x,a4,1x,a,":"' // int2str(maxcharsname+3-len(trim(elemstepnames(i)))) // 'x,)', &
                    advance='no') trim(int2str(i)) // '.', trim(elemstepnames(i))
    
    call calculate_preex_factor(i,temp,preexpfac)
    
    if (reverselemstep(i) == 0 .or. reverselemstep(i) > i) then ! forward step of a reversible reaction or irreversible step
        write(iwrite,'(1x,"A(Tini) = ",ES11.4E2,3x,"Ea = ",F5.2,3x)',advance='no')  preexpfac, acteng(i)
        ktemp = preexpfac*dexp(-acteng(i)/enrgconv/kboltz/temp)
		
		do j = 1,elemstepgases(i,0)
            if (elemstepgases(i,2*j) < 0) then
                kspec = elemstepgases(i,2*j-1)
                ktemp = ktemp*gasmolfracs(j)
            endif
        enddo
		
        write(iwrite,'(1x,"k(Tini) = ",ES11.4E2,2x,"Reaction:",1x)',advance='no')  ktemp
    
    else
        write(iwrite,'(1x,"A(Tini) = ",ES11.4E2,40x,"Reaction:",1x)',advance='no')  preexpfac
    endif
    

    ! Reactants
    ! Gas
    do j = 1,elemstepgases(i,0)
        if (elemstepgases(i,2*j) < 0) then
            kspec = elemstepgases(i,2*j-1)
            if (elemstepgases(i,2*j) > 1) then
                write(iwrite,'(a,1x)',advance='no') trim(int2str(elemstepgases(i,2*j)))
            endif
            write(iwrite,'(a,"  +  ")',advance='no') trim(gasspecsnames(kspec))
        endif
    enddo
    ! Surface
    do j = 1,elemstepreactnts(i,0)
        kspec = elemstepreactnts(i,j)
        write(iwrite,'(a)',advance='no') trim(surfspecsnames(kspec)) // '('
        do kdent = 1,surfspecsdent(kspec)
            istype = elemstepstype(i,elemstepadsorbposin(i,j,kdent))
            if (istype > 0) then
                write(iwrite,'(a)',advance='no') trim(sitetypenames(istype))
            else
                write(iwrite,'(a)',advance='no') '.'
            endif
            if (kdent == surfspecsdent(kspec)) then
                write(iwrite,'(a)',advance='no') ')'
            else
                write(iwrite,'(a)',advance='no') ','
            endif
        enddo
        if (j == elemstepreactnts(i,0)) then
            write(iwrite,'(a)',advance='no') '  ->  '
        else
            write(iwrite,'(a)',advance='no') '  +  '
        endif
    enddo
    
    ! Products
    ! Gas
    do j = 1,elemstepgases(i,0)
        if (elemstepgases(i,2*j) > 0) then
            kspec = elemstepgases(i,2*j-1)
            if (elemstepgases(i,2*j) > 1) then
                write(iwrite,'(a,1x)',advance='no') trim(int2str(elemstepgases(i,2*j)))
            endif
            write(iwrite,'(a,"  +  ")',advance='no') trim(gasspecsnames(kspec))
        endif
    enddo
    ! Surface
    do j = 1,elemstepproducts(i,0)
        kspec = elemstepproducts(i,j)
        write(iwrite,'(a)',advance='no') trim(surfspecsnames(kspec)) // '('
        do kdent = 1,surfspecsdent(kspec)
            istype = elemstepstype(i,elemstepadsorbposfn(i,j,kdent))
            if (istype > 0) then
                write(iwrite,'(a)',advance='no') trim(sitetypenames(istype))
            else
                write(iwrite,'(a)',advance='no') '.'
            endif
            if (kdent == surfspecsdent(kspec)) then
                write(iwrite,'(a)',advance='no') ')'
            else
                write(iwrite,'(a)',advance='no') ','
            endif
        enddo
        if (j /= elemstepproducts(i,0)) then
            write(iwrite,'(a)',advance='no') '  +  '
        endif
    enddo

enddo

write(iwrite,'(a)') ''
write(iwrite,'(/,a)') 'Finished reading mechanism input.'

return

end subroutine read_mechanism_setup

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine parse_elemstep(firstpass,recinput,nwords,irec,io,istep,iw,stepopen)

implicit none

integer i, iadsorb, identat, io, irec, istep, ispec, k, m
integer elemstepngas, nv(25), tmpv(2), nsites
integer nwords, nwords2, ivariant, j, j1, j2
integer iw(maxwords2), iw2(maxwords2)

character(lengthrecinp) recinput, tmpstr
character(nnam0) basename, varname

logical stepopen, isreversible
logical readgasrctprod, readnsites, readstypes, readinitlstate, readfinalstate, readneighb
logical readpreexp, readactiveng, readomega, readperatio, canspecifyvariant, firstpass !, readdeltah

! Parses a single elementary step in simulation_input.dat

! Level-1 keywords allowed in this subroutine:
!  -Normal keywords recognized in this scope
!    gas_reacs_prods
!    sites
!    initial
!    final
!    neighboring
!    end_step
!    end_reversible_step
!  -The following keywords direct the program flow to subroutine parse_stypes_ratepars
!    site_types
!    pre_expon
!    pe_ratio
!    activ_eng
!    delta_eng
!    energetics
!    energetics_rev
!  -The following keywords direct the program flow to subroutine parse_variant
!    variant
!  -The following keywords are recognized as invalid in this context and produce and error
!    step
!    reversible_step

! Depending on the value of the logical variable firstpass the subroutine can be run in one of 
! two "modes": a first pass mode and a full parsing mode. The former just checks for consistency
! of input and finds the max bounds needed for the arrays that encode the mechanism. The latter
! populates these arrays. The user MUST run a first pass of the mechanism before the full parsing.

readgasrctprod = .false.
readnsites = .false.
readstypes = .false.
readinitlstate = .false.
readfinalstate = .false.
readneighb = .false.
readpreexp = .false.
readactiveng = .false.
readperatio = .false.
canspecifyvariant = .true.

istep = istep + 1
ivariant = 0
if (nwords == 2) then
    basename = recinput(iw(3):iw(4))
elseif (nwords == 1) then
    basename = 'Elem_Step_' // trim(int2str(istep))
else
    call error(3002)
endif

if (striccompare(recinput(iw(1):iw(2)),'reversible_step')) then
    isreversible = .true.
    nreversible = nreversible + 1
    if (.not.firstpass) then
        reverselemstep(istep) = istep+1
        reverselemstep(istep+1) = istep
        elemstepnames(istep) = trim(basename) // '_fwd'
        elemstepnames(istep+1) = trim(basename) // '_rev'
    endif
else
    isreversible = .false.
    if (.not.firstpass) then
            reverselemstep(istep) = 0
            elemstepnames(istep) = basename
    endif
endif

stepopen = .true.
do while (io >= 0)

    call getrecord(imechread,recinput,irec,io)
    
    call break_words(recinput,' ',remchar,nwords,iw)

    if (nwords == 0) then
        cycle
    endif
    
    if (striccompare(recinput(iw(1):iw(2)),'gas_reacs_prods')) then
                            
        ! If a first pass is requested ...
        if (firstpass) then

            if (readgasrctprod) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(cmechfname) // '.'
                call error(100)
            endif

            elemstepngas = (nwords-1)/2
            
            tmpstr = 'AA'
            do i = 1,elemstepngas
                tmpstr = trim(tmpstr) // '/AA/I4'
            enddo

            if (str_check_expression(recinput,nwords,iw,trim(tmpstr),.true.)) then 
                ! ... we only check input consistency and update elemstepnmxgas.
                elemstepnmxgas = max(elemstepngas,elemstepnmxgas)
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(3003)
            endif
            
            readgasrctprod = .true.
            cycle
            
        endif
        
        ! Full parsing is requested:
        elemstepgases(istep,0) = (nwords-1)/2
        do i = 1,elemstepgases(istep,0)
            ispec = findinstrv(gasspecsnames,1,1,ngasspecs,recinput(iw(4*i-1):iw(4*i)))
            if (ispec < 0) then
                moreinfo = 'Encountered unknown gas species "' // recinput(iw(4*i-1):iw(4*i)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(3004)    
            endif
            elemstepgases(istep,2*i-1) = ispec
            elemstepgases(istep,2*i) = str2int(recinput(iw(4*i+1):iw(4*i+2)))

            ! If the step is reversible, add info for the (istep+1)th step
            ! for which the gas species coefficient is the opposite of that of the current step
            if (isreversible) then
                elemstepgases(istep+1,0) = elemstepgases(istep,0)
                elemstepgases(istep+1,2*i-1) = ispec
                elemstepgases(istep+1,2*i) = -elemstepgases(istep,2*i)
            endif

        enddo            
        
    elseif (striccompare(recinput(iw(1):iw(2)),'sites')) then

        ! If a first pass is requested ...
        if (firstpass) then

            if (readnsites) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(cmechfname) // '.'
                call error(100)
            endif
        
            if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
                nsites = str2int(recinput(iw(3):iw(4)))
                elemstepnmxsites = max(nsites,elemstepnmxsites)
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(3005)
            endif

            readnsites = .true.
            cycle
            
        endif

        ! Full parsing is requested:

        elemstepnsites(istep) = str2int(recinput(iw(3):iw(4)))
        
        ! For the reverse step, same number of sites
        if (isreversible) then
            elemstepnsites(istep+1) = elemstepnsites(istep)
        endif

        readnsites = .true.
        
    elseif (striccompare(recinput(iw(1):iw(2)),'initial')) then
        
        ! For a first pass ...
        if (firstpass) then

            if (readinitlstate) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(cmechfname) // '.'
                call error(100)
            endif
            
            ! ... we check if this keyword is preceded by the sites keyword and then
            ! skip this part of the input ...
            if (.not.readnsites) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call error(3006)
            endif
            
            j = 0
            do while (j < nsites .and. io >= 0)
            
                call getrecord(imechread,recinput,irec,io)
                call break_words(recinput,' ',remchar,nwords,iw)
                if (nwords == 0) cycle
                
                j = j + 1

                if (.not.str_check_expression(recinput,nwords,iw,'I4/AA/I4',.true.)) then
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(3007)
                endif
                
            enddo
            
            readinitlstate = .true.
            cycle
            
        endif
                
        ! Full parsing requested...
        j = 0
        do while (j < elemstepnsites(istep) .and. io >= 0)
        
            call getrecord(imechread,recinput,irec,io)
            call break_words(recinput,' ',remchar,nwords,iw)
            if (nwords == 0) cycle
            
            j = j + 1
            
            iadsorb = str2int(recinput(iw(1):iw(2)))
            ispec = findinstrv(surfspecsnames,0,0,nsurfspecs,recinput(iw(3):iw(4)))
            if (ispec < 0) then
                moreinfo = 'Encountered unknown surface species "' // recinput(iw(3):iw(4)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(3004)    
            endif
            identat = str2int(recinput(iw(5):iw(6)))

            if (identat < 0 .or. identat > surfspecsdent(ispec)) then
                moreinfo = 'Dentate number ' // trim(int2str(identat)) // ' for species ' // recinput(iw(3):iw(4)) // &
                           ' is out of bounds 1:' // trim(int2str(surfspecsdent(ispec))) // ' in line ' // trim(int2str(irec)) // '.'
                call error(3028)
            endif

            elemsteplatticestatein(istep,j,1) = iadsorb
            elemsteplatticestatein(istep,j,2) = identat
            
            if (elemstepreactnts(istep,iadsorb) >= 0) then
            
                ! Check if this adsorbate/entity has previously been encountered
                if (elemstepreactnts(istep,iadsorb) /= ispec) then
                    moreinfo = 'Adsorbate/entity ' // trim(int2str(iadsorb)) // ' has previously been declared as ' // &
                               trim(surfspecsnames(elemstepreactnts(istep,iadsorb))) // & 
                               ' species but now appears as ' // trim(surfspecsnames(ispec)) // &
                               ' in line ' // trim(int2str(irec)) // '.'
                    call error(3008)
                endif
            else
                
                elemstepreactnts(istep,0) = elemstepreactnts(istep,0) + 1
                elemstepreactnts(istep,iadsorb) = ispec                            
                
            endif
            
            ! The inital state of the forward step is the final state of the reverse step
            if (isreversible) then
                elemsteplatticestatefn(istep+1,j,1) = iadsorb
                elemsteplatticestatefn(istep+1,j,2) = identat
                elemstepproducts(istep+1,0) = elemstepreactnts(istep,0)
                elemstepproducts(istep+1,iadsorb) = ispec                            
            endif

        enddo
        
        readinitlstate = .true.
        
    elseif (striccompare(recinput(iw(1):iw(2)),'final')) then

        ! For a first pass ...
        if (firstpass) then

            if (readfinalstate) then
                moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
                'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(cmechfname) // '.'
                call error(100)
            endif
            
            ! ... we check if this keyword is preceded by the sites keyword and then
            ! skip this part of the input ...
            if (.not.readnsites) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call error(3006)
            endif
            
            j = 0
            do while (j < nsites .and. io >= 0)
            
                call getrecord(imechread,recinput,irec,io)
                call break_words(recinput,' ',remchar,nwords,iw)
                if (nwords == 0) cycle
                
                j = j + 1

                if (.not.str_check_expression(recinput,nwords,iw,'I4/AA/I4',.true.)) then
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(3007)
                endif
                
            enddo
            
            readfinalstate = .true.
            cycle
            
        endif
                
        ! Full parsing requested...
        j = 0
        do while (j < elemstepnsites(istep) .and. io >= 0)
        
            call getrecord(imechread,recinput,irec,io)
            call break_words(recinput,' ',remchar,nwords,iw)
            if (nwords == 0) cycle
            
            j = j + 1
            
            iadsorb = str2int(recinput(iw(1):iw(2)))
            ispec = findinstrv(surfspecsnames,0,0,nsurfspecs,recinput(iw(3):iw(4)))
            if (ispec < 0) then
                moreinfo = 'Encountered unknown surface species "' // recinput(iw(3):iw(4)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(3004)    
            endif
            identat = str2int(recinput(iw(5):iw(6)))

            if (identat < 0 .or. identat > surfspecsdent(ispec)) then
                moreinfo = 'Dentate number ' // trim(int2str(identat)) // ' for species ' // recinput(iw(3):iw(4)) // &
                           ' is out of bounds 1:' // trim(int2str(surfspecsdent(ispec))) // ' in line ' // trim(int2str(irec)) // '.'
                call error(3028)
            endif

            elemsteplatticestatefn(istep,j,1) = iadsorb
            elemsteplatticestatefn(istep,j,2) = identat
            
            if (elemstepproducts(istep,iadsorb) >= 0) then
            
                ! Check if this adsorbate/entity has previously been encountered
                if (elemstepproducts(istep,iadsorb) /= ispec) then
                    moreinfo = 'Adsorbate/entity ' // trim(int2str(iadsorb)) // ' has previously been declared as ' // &
                               trim(surfspecsnames(elemstepproducts(istep,iadsorb))) // & 
                               ' species but now appears as ' // trim(surfspecsnames(ispec)) // &
                               ' in line ' // trim(int2str(irec)) // '.'
                    call error(3008)
                endif
            else
                
                elemstepproducts(istep,0) = elemstepproducts(istep,0) + 1
                elemstepproducts(istep,iadsorb) = ispec                            
                
            endif
            
            ! The final state of the forward step is the inital state of the reverse step
            if (isreversible) then
                elemsteplatticestatein(istep+1,j,1) = iadsorb
                elemsteplatticestatein(istep+1,j,2) = identat
                elemstepreactnts(istep+1,0) = elemstepproducts(istep,0)
                elemstepreactnts(istep+1,iadsorb) = ispec                            
            endif

        enddo
        
        readfinalstate = .true.
        
    elseif (striccompare(recinput(iw(1):iw(2)),'neighboring')) then
        
        ! For a first pass...
        if (firstpass) then
        
            nv = 0 ! remark: nv is a vector
            
            do j = 1,nwords-1
                
                tmpstr = recinput(iw(2*j+1):iw(2*j+2))
                call break_words(trim(tmpstr),'-',remchar,nwords2,iw2)
                
                if (.not.str_check_expression(tmpstr,nwords2,iw2,'2I4',.true.)) then
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(3009)
                endif                    
                
                j1 = str2int(tmpstr(iw2(1):iw2(2)))
                j2 = str2int(tmpstr(iw2(3):iw2(4)))
                
                if (j1 == j2) then
                    moreinfo = 'Encountered "' // recinput(iw2(1):iw2(2)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(3010)
                endif
                
                tmpv = (/j1,j2/)
                do k = 1,2                                    
                    if (tmpv(k) > nsites .or. tmpv(k) < 1) then
                        moreinfo = 'Site number ' // trim(int2str(tmpv(k))) // &
                                   ' referenced in neighboring input in line ' // trim(int2str(irec)) // &
                                   ' is out of range 1:' // trim(int2str(nsites)) // &
                                   ' for this elementary step.'
                        call error(3009)
                    endif
                    
                enddo
                
                nv(j1) = nv(j1) + 1
                nv(j2) = nv(j2) + 1
                elemstepnmxcoord = max(elemstepnmxcoord,nv(j1),nv(j2))
                
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
            
            elemstepneigh(istep,j1,0) = elemstepneigh(istep,j1,0) + 1
            elemstepneigh(istep,j1,elemstepneigh(istep,j1,0)) = j2

            elemstepneigh(istep,j2,0) = elemstepneigh(istep,j2,0) + 1
            elemstepneigh(istep,j2,elemstepneigh(istep,j2,0)) = j1
            
            ! Check if the j1-j2 pair has already been encountered in the neighboring list
            do k = 1,elemstepneigh(istep,j1,0)
                do m = k+1,elemstepneigh(istep,j1,0)
                    if (elemstepneigh(istep,j1,k) == elemstepneigh(istep,j1,m)) then
                        moreinfo = 'Pair ' // trim(int2str(j1)) // '-' // trim(int2str(j2)) // &
                                   ' is referenced twice in neighboring input in line ' // trim(int2str(irec)) // '.'
                        call error(3011)
                    endif
                enddo
            enddo
            
            ! Neighboring is the same for the reverse step
            if (isreversible) then
                elemstepneigh(istep+1,j1,0) = elemstepneigh(istep,j1,0)
                elemstepneigh(istep+1,j1,elemstepneigh(istep+1,j1,0)) = j2

                elemstepneigh(istep+1,j2,0) = elemstepneigh(istep,j2,0)
                elemstepneigh(istep+1,j2,elemstepneigh(istep+1,j2,0)) = j1
            endif
            
        enddo
                
        ! Check whether a site has been left unconnected
        do j = 1,elemstepnsites(istep)
            if (elemstepneigh(istep,j,0) == 0) then
                moreinfo = 'Site number ' // trim(int2str(j)) // &
                           ' is not referenced in neighboring input in line ' // trim(int2str(irec)) // '.'
                call error(3012)
            endif
        enddo
        
        readneighb = .true.
        
    elseif (striccompare(recinput(iw(1):iw(2)),'site_types') .or. &
            striccompare(recinput(iw(1):iw(2)),'pre_expon') .or. &
            striccompare(recinput(iw(1):iw(2)),'pe_ratio') .or. &
            striccompare(recinput(iw(1):iw(2)),'prox_factor') .or. &
            striccompare(recinput(iw(1):iw(2)),'activ_eng') .or. &
            striccompare(recinput(iw(1):iw(2)),'angles') .or. &
            striccompare(recinput(iw(1):iw(2)),'no_mirror_images') .or. &
            striccompare(recinput(iw(1):iw(2)),'absl_orientation')) then
        
        if (ivariant > 0) then
            moreinfo = 'Invalid specification in line ' // trim(int2str(irec)) // '.'
            call error(3033)
        endif
        
        canspecifyvariant = .false.
                
        call parse_stypes_ratepars(firstpass,recinput,irec,io,nwords,iw,istep,nsites,isreversible, &
                             readnsites,readneighb,readstypes,readpreexp,&
                             readactiveng,readomega,readperatio)
 
    elseif (striccompare(recinput(iw(1):iw(2)),'variant')) then
        
        ivariant = ivariant + 1

        if (firstpass) then

            if (.not.canspecifyvariant) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call error(3027)
            endif
                    
            if (ivariant == 1) then
                ! If this is the first variant found, check whether the information
                ! about sites, initial and final state and neighboring structure has been read
                call check_level1_elemstep_completeness(irec,istep,basename,nsites, &
                        readnsites,readinitlstate,readfinalstate,readneighb)
            endif
        
            readpreexp = .false.
            readactiveng = .false.
            readperatio = .false.
        
        endif        
        
        call parse_variant(firstpass,ivariant,varname,basename,recinput,irec,io,nwords,iw,istep,nsites, &
               isreversible,readnsites,readneighb,readstypes,readpreexp,readactiveng,readomega,readperatio)
        
        if (firstpass) then
            call check_level2_elemstep_completeness(irec,istep,basename,ivariant,varname, &
                   isreversible,readpreexp,readactiveng,readperatio)
        endif
        
    elseif ( ( striccompare(recinput(iw(1):iw(2)),'end_step') .or. &
               striccompare(recinput(iw(1):iw(2)),'end_reversible_step') ) .and. stepopen ) then
        
        stepopen = .false.
        
        if (firstpass) then
            if ( (      isreversible .and. striccompare(recinput(iw(1):iw(2)),'end_step')           ) .or. &
                 ( .not.isreversible .and. striccompare(recinput(iw(1):iw(2)),'end_reversible_step')) ) then
                moreinfo = 'Encountered ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
                call error(3022)
            endif

            if (ivariant > 0) then
                nelemsteps = nelemsteps + ivariant
                if (isreversible) nelemsteps = nelemsteps + ivariant
            else
                nelemsteps = nelemsteps + 1
                if (isreversible) nelemsteps = nelemsteps + 1
            endif
            
        endif
        
        exit
    
    elseif (striccompare(recinput(iw(1):iw(2)),'step') .or. &
            striccompare(recinput(iw(1):iw(2)),'reversible_step')) then

        moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
        call error(3023)

    else
        moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
        call error(3100)
    endif
        
enddo ! loop for reading info for one elementary step

if (ivariant == 0 .and. firstpass) then   
    call check_level1_elemstep_completeness(irec,istep,basename,nsites, &
            readnsites,readinitlstate,readfinalstate,readneighb)
    call check_level2_elemstep_completeness(irec,istep,basename,ivariant,varname, &
           isreversible,readpreexp,readactiveng,readperatio)
endif

if (isreversible) then
    istep = istep + 1
endif

return

end subroutine parse_elemstep

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine parse_stypes_ratepars(firstpass,recinput,irec,io,nwords,iw,istep,nsites,isreversible, &
                                 readnsites,readneighb, &                                 
                                 readstypes,readpreexp,readactiveng,readomega,readperatio)

implicit none

integer io, irec, istep, istype, nsites
integer j, k, m, p, nwords
integer iw(maxwords2)
integer nwords2
integer iw2(maxwords2)
integer nwords3
integer iw3(maxwords2)
integer sites123(3)

character(lengthrecinp) recinput, tmpstr1, tmpstr2

logical firstpass, isreversible, readnsites, readneighb
logical readpreexp, readactiveng, readomega, readperatio, readstypes, neighfound
logical tmplog1, tmplog2

! Level-1 keywords allowed in this subroutine:
!  -Normal keywords recognized in this scope
!    site_types
!    pre_expon
!    pe_ratio
!    activ_eng

if (striccompare(recinput(iw(1):iw(2)),'site_types')) then
    
    if (firstpass) then
        
        if (nwords-1 == nsites) then
            readstypes = .true.
            return
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(3013)
        endif        
        
    endif
    
    do j = 1,nwords-1
        
        if (str_is_integer(recinput(iw(2*j+1):iw(2*j+2)),4)) then
            istype = str2int(recinput(iw(2*j+1):iw(2*j+2)))
            if (istype > nsitetypes) then
                moreinfo = 'Encountered site type ' // trim(int2str(istype)) // ' for site ' // trim(int2str(j)) // &
                           ' in line ' // trim(int2str(irec)) // '.'
                call error(3014)
            endif
        else
            istype = findinstrv(sitetypenames,1,1,nsitetypes,recinput(iw(2*j+1):iw(2*j+2)))                                
            if (istype < 0) then
                moreinfo = 'Encountered site type name "' // recinput(iw(2*j+1):iw(2*j+2)) // '" for site ' // trim(int2str(j)) // &
                           ' in line ' // trim(int2str(irec)) // '.'
                call error(3014)
            endif
        endif
        
        elemstepstype(istep,j) = istype
        
        ! Same site types for reverse step
        if (isreversible) then
            elemstepstype(istep+1,j) = istype
        endif
                                    
    enddo
    
    readstypes = .true.
    
    return
    
elseif (striccompare(recinput(iw(1):iw(2)),'pre_expon')) then
    
    if (firstpass) then

        if ( str_check_expression(recinput,nwords,iw,'AA/R8',.true.) .or. &
             str_check_expression(recinput,nwords,iw,'AA/7R8',.true.) ) then
            readpreexp = .true.
            return
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(3015)
        endif
        
    endif
    
    if (str_check_expression(recinput,nwords,iw,'AA/R8',.true.)) then
        
        if (readperatio .and. (preexpisconst(istep) .eqv. .false.) ) then
            moreinfo = 'Pre-exponential and pre-exponential ratios must both be specified as either constants or ' // & 
                       'temperature dependent functions. Occurred while parsing in line ' // trim(int2str(irec)) // '.'
            call error(3100)
        endif
        
        preexpisconst(istep) = .true.
        preexp(istep,0) = str2dbl(recinput(iw(3):iw(4)))
        
        ! The reverse pre-exponential will be computed as Afwd/PEratio
        if (isreversible) then
            preexpisconst(istep+1) = .true.
            preexp(istep+1,0) = preexp(istep+1,0)*preexp(istep,0)
        endif

        readpreexp = .true.
        
    else
        
        if (readperatio .and. (preexpisconst(istep) .eqv. .true.) ) then
            moreinfo = 'Pre-exponential and pre-exponential ratios must both be specified as either constants or ' // & 
                       'temperature dependent functions. Occurred while parsing in line ' // trim(int2str(irec)) // '.'
            call error(3100)
        endif
        
        preexpisconst(istep) = .false.
        preexp(istep,1) = str2dbl(recinput(iw(3):iw(4)))
        preexp(istep,2) = str2dbl(recinput(iw(5):iw(6)))
        preexp(istep,3) = str2dbl(recinput(iw(7):iw(8)))
        preexp(istep,4) = str2dbl(recinput(iw(9):iw(10)))
        preexp(istep,5) = str2dbl(recinput(iw(11):iw(12)))
        preexp(istep,6) = str2dbl(recinput(iw(13):iw(14)))
        preexp(istep,7) = str2dbl(recinput(iw(15):iw(16)))
                    
        ! The reverse pre-exponential will be computed as Afwd/PEratio
        if (isreversible) then
            preexpisconst(istep+1) = .false.
            do j = 1,7
                preexp(istep+1,j) = preexp(istep+1,j) + preexp(istep,j)
            enddo
        endif

        readpreexp = .true.

    endif
    
    return

elseif (striccompare(recinput(iw(1):iw(2)),'pe_ratio')) then
   
    if (firstpass) then

        if (.not.isreversible) then
            moreinfo = 'Invalid keyword "' // recinput(iw(1):iw(2)) // '" in line ' // trim(int2str(irec)) // '.'       
            call error(3016)
        endif

        if (str_check_expression(recinput,nwords,iw,'AA/R8',.true.) .or. &
             str_check_expression(recinput,nwords,iw,'AA/7R8',.true.) ) then
            readperatio = .true.
            return
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(3017)
        endif
        
    endif

    if (str_check_expression(recinput,nwords,iw,'AA/R8',.true.)) then
        
        if (readpreexp .and. (preexpisconst(istep+1) .eqv. .false.)) then
            moreinfo = 'Pre-exponential and pre-exponential ratios must both be specified as either constants or ' // & 
                       'temperature dependent functions. Occurred while parsing in line ' // trim(int2str(irec)) // '.'
            call error(3100)
        endif

        preexpisconst(istep) = .true.
        preexpisconst(istep+1) = .true.
        preexp(istep+1,0) = preexp(istep+1,0)/str2dbl(recinput(iw(3):iw(4)))
                
        readperatio = .true.
    
    else
    
        if (readpreexp .and. (preexpisconst(istep) .eqv. .true.)) then
            moreinfo = 'Pre-exponential and pre-exponential ratios must both be specified as either constants or ' // & 
                       'temperature dependent functions. Occurred while parsing in line ' // trim(int2str(irec)) // '.'
            call error(3100)
        endif
        
        preexpisconst(istep) = .false.
        preexpisconst(istep+1) = .false.
        preexp(istep+1,1) = preexp(istep+1,1) - str2dbl(recinput(iw(3):iw(4)))
        preexp(istep+1,2) = preexp(istep+1,2) - str2dbl(recinput(iw(5):iw(6)))
        preexp(istep+1,3) = preexp(istep+1,3) - str2dbl(recinput(iw(7):iw(8)))
        preexp(istep+1,4) = preexp(istep+1,4) - str2dbl(recinput(iw(9):iw(10)))
        preexp(istep+1,5) = preexp(istep+1,5) - str2dbl(recinput(iw(11):iw(12)))
        preexp(istep+1,6) = preexp(istep+1,6) - str2dbl(recinput(iw(13):iw(14)))
        preexp(istep+1,7) = preexp(istep+1,7) - str2dbl(recinput(iw(15):iw(16)))
            
    endif
    
    return
    
elseif (striccompare(recinput(iw(1):iw(2)),'activ_eng')) then

    if (firstpass) then

        if (str_check_expression(recinput,nwords,iw,'AA/R8',.true.)) then
            readactiveng = .true.
            return
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(3018)
        endif
        
    endif
        
    acteng(istep) = str2dbl(recinput(iw(3):iw(4)))

    if (acteng(istep) < 0.d0) then
        moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
        call warning(3029)
    endif

    readactiveng = .true.
    
    ! MA -- The following is not portable. It doesbn't look like it's used as a
    ! marker anyhow, so it should be safe to comment  out.
    !  ! The reverse activation energy will be computed as Efwd - DeltaH
    !  if (isreversible) then
    !      acteng(istep+1) = 0.d0/0.d0 ! NaN for reverse activation energy
    !      ! since everything will be computed by the energetics model and the forward activation energy
    !  endif

    return

elseif (striccompare(recinput(iw(1):iw(2)),'prox_factor')) then

    if (firstpass) then

        if (str_check_expression(recinput,nwords,iw,'AA/R8',.true.)) then
            readactiveng = .true.
            return
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(3031)
        endif
        
    endif
        
    omega(istep) = str2dbl(recinput(iw(3):iw(4)))

    if (omega(istep) < 0.d0 .or. omega(istep) > 1.d0) then
        moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
        call error(3032)
    endif

    readomega = .true.
    
    return

elseif (striccompare(recinput(iw(1):iw(2)),'angles')) then
    
    if (firstpass) then

        if (.not. readnsites) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(3034)
        endif

        if (.not. readneighb) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(3035)
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
                elemstepnmxangles = elemstepnmxangles + 1
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(3036)
            endif

        enddo
        
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
            if (sites123(k) < 1 .or. sites123(k) > elemstepnsites(istep)) then
                moreinfo = 'Site number ' // trim(int2str(sites123(k))) // &
                           ' referenced in angle input in line ' // trim(int2str(irec)) // &
                           ' is out of range 1:' // trim(int2str(nsites)) // &
                           ' for this elementary step.'
                call error(3037)
            endif
        enddo
        
        if (sites123(1) == sites123(2) .or. sites123(1) == sites123(3) .or. sites123(2) == sites123(3)) then
            moreinfo = 'Invalid repetition of sites in angle specification for an elementary step in line ' // trim(int2str(irec))
            call error(3037)
        endif
        
        ! See if the neighboring implied by the angle specification is compatible with that of the elementary step
        do k = 1,2
            m = k+1
            neighfound = .false.
            do p = 1,elemstepneigh(istep,sites123(k),0)
                if (sites123(m) == elemstepneigh(istep,sites123(k),p)) then
                    neighfound = .true.
                    exit
                endif
            enddo
            
            if (.not.neighfound) then
                moreinfo = 'Site number ' // trim(int2str(sites123(k))) // &
                           ' neighbors with site number ' // trim(int2str(sites123(m))) // & 
                           ' in angle input in line ' // trim(int2str(irec)) // '.' // &
                           ' This is invalid for this elementary step.'
                call error(3038)
            endif
                
        enddo
        
        elemstepnanglessequenc(istep,0) = elemstepnanglessequenc(istep,0) + 1
        m = elemstepnanglessequenc(istep,0)
        elemstepnanglessequenc(istep,3*m-2) = sites123(1)
        elemstepnanglessequenc(istep,3*m-1) = sites123(2)
        elemstepnanglessequenc(istep,3*m) = sites123(3)
        
        elemstepangles(istep,m) = str2dbl(tmpstr1(iw2(3):iw2(4)))

        ! Angle definition same for the reverse step
        if (isreversible) then
            elemstepnanglessequenc(istep+1,0) = m
            elemstepnanglessequenc(istep+1,3*m-2) = elemstepnanglessequenc(istep,3*m-2)
            elemstepnanglessequenc(istep+1,3*m-1) = elemstepnanglessequenc(istep,3*m-1)
            elemstepnanglessequenc(istep+1,3*m) = elemstepnanglessequenc(istep,3*m)
            elemstepangles(istep+1,m) = elemstepangles(istep,m)
        endif
        
    enddo

elseif (striccompare(recinput(iw(1):iw(2)),'no_mirror_images')) then

    if (firstpass) then
        return
    endif
    
    elemstepnomirrorimgs(istep) = .true.
    
    ! Same setting for the reverse step
    if (isreversible) then
        elemstepnomirrorimgs(istep+1) = .true.
    endif

elseif (striccompare(recinput(iw(1):iw(2)),'absl_orientation')) then
    
    if (firstpass) then

        if (.not. readnsites) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(3040)
        endif

        if (.not. readneighb) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(3041)
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
            call error(3042)
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
        if (sites123(k) < 1 .or. sites123(k) > elemstepnsites(istep)) then
            moreinfo = 'Site number ' // trim(int2str(sites123(k))) // &
                        ' referenced in absolute orientation input in line ' // trim(int2str(irec)) // &
                        ' is out of range 1:' // trim(int2str(nsites)) // &
                        ' for this elementary step.'
            call error(3043)
        endif
    enddo
        
    if (sites123(1) == sites123(2)) then
        moreinfo = 'Invalid repetition of sites in absolute orientation specification for an elementary step in line ' // trim(int2str(irec))
        call error(3043)
    endif
        
    ! See if the neighboring implied by the angle specification is compatible with that of the elementary step
    neighfound = .false.
    do p = 1,elemstepneigh(istep,sites123(1),0)
        if (sites123(2) == elemstepneigh(istep,sites123(1),p)) then
            neighfound = .true.
            exit
        endif
    enddo
            
    if (.not.neighfound) then
        moreinfo = 'Site number ' // trim(int2str(sites123(k))) // &
                    ' neighbors with site number ' // trim(int2str(sites123(m))) // & 
                    ' in absolute orientation input in line ' // trim(int2str(irec)) // '.' // &
                    ' This is invalid for this elementary step.'
        call error(3043)
    endif
                        
    elemstepabslorientat(istep) = .true.
    elemsteporientationedge(istep,1) = sites123(1)
    elemsteporientationedge(istep,2) = sites123(2)
        
    elemsteporientationangles(istep) = str2dbl(tmpstr1(iw2(3):iw2(4)))

    ! Same absolute orientation setting for the reverse step
    if (isreversible) then
        
        elemstepabslorientat(istep+1) = elemstepabslorientat(istep)
        
        elemsteporientationedge(istep+1,1) = elemsteporientationedge(istep,1)
        elemsteporientationedge(istep+1,2) = elemsteporientationedge(istep,2)

        elemsteporientationangles(istep+1) = elemsteporientationangles(istep)
        
    endif

else
    
    moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
    call error(3100)
    
endif

end subroutine parse_stypes_ratepars

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine parse_variant(firstpass,ivariant,varname,basename,recinput,irec,io,nwords,iw,istep,nsites, &
                         isreversible,readnsites,readneighb,readstypes,readpreexp, &
                         readactiveng,readomega,readperatio) !,readdeltah

implicit none

integer i, io, irec, istep, nsites
integer ivariant, j, k, k1, k2, dk, nwords
integer iw(maxwords2)

character(lengthrecinp) recinput
character(nnam0) basename, varname

logical firstpass, isreversible, readpreexp, readactiveng, readomega
logical readperatio, readstypes, readnsites, readneighb

if (.not. firstpass) then

    if (ivariant > 1) then
        if (isreversible) then
            k1 = istep
            k2 = istep+1
            dk = 2
            istep = istep + 2
            reverselemstep(istep) = istep+1
            reverselemstep(istep+1) = istep
        else
            k1 = istep
            k2 = istep
            dk = 1
            istep = istep + 1   
            reverselemstep(istep) = 0
        endif

        do k = k1,k2
        
            preexpisconst(k+dk) = preexpisconst(k)
        
            elemstepnsites(k+dk) = elemstepnsites(k)
            
            do i = 1,elemstepnsites(k)
                do j = 0,elemstepneigh(k,i,0)
                    elemstepneigh(k+dk,i,j) = elemstepneigh(k,i,j)
                enddo
            enddo
            
            do i = 1,elemstepnsites(k)
                do j = 1,2
                    elemsteplatticestatein(k+dk,i,j) = elemsteplatticestatein(k,i,j)
                    elemsteplatticestatefn(k+dk,i,j) = elemsteplatticestatefn(k,i,j)
                enddo
            enddo
            
            do i = 0,elemstepnsites(k)
                elemstepreactnts(k+dk,i) = elemstepreactnts(k,i)
                elemstepproducts(k+dk,i) = elemstepproducts(k,i)
            enddo
            
            do i = 0,2*elemstepgases(k,0)
                elemstepgases(k+dk,i) = elemstepgases(k,i)
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
    if (isreversible) then
        elemstepnames(istep) = trim(basename) // '_fwd_' // varname
        elemstepnames(istep+1) = trim(basename) // '_rev_' // varname
    else
        elemstepnames(istep) = trim(basename) // '_' // varname
    endif
endif

do while (io >= 0)

    call getrecord(imechread,recinput,irec,io)
    
    call break_words(recinput,' ',remchar,nwords,iw)

    if (nwords > 0) then

        if (striccompare(recinput(iw(1):iw(2)),'site_types') .or. &
            striccompare(recinput(iw(1):iw(2)),'pre_expon') .or. &
            striccompare(recinput(iw(1):iw(2)),'pe_ratio') .or. &
            striccompare(recinput(iw(1):iw(2)),'prox_factor') .or. &
            striccompare(recinput(iw(1):iw(2)),'activ_eng') .or. &
            striccompare(recinput(iw(1):iw(2)),'angles') .or. &
            striccompare(recinput(iw(1):iw(2)),'no_mirror_images') .or. &
            striccompare(recinput(iw(1):iw(2)),'absl_orientation')) then

            call parse_stypes_ratepars(firstpass,recinput,irec,io,nwords,iw,istep,nsites,isreversible, &
                                       readnsites,readneighb,readstypes,readpreexp,readactiveng, &
                                       readomega,readperatio) !,readdeltah
        
        elseif (striccompare(recinput(iw(1):iw(2)),'end_variant')) then
            
            return
            
        else
            
            moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
            call error(3100)
            
        endif
                
    endif
    
enddo

return

end subroutine parse_variant

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_level1_elemstep_completeness(irec,istep,elstepname,nsites, &
            readnsites,readinitlstate,readfinalstate,readneighb)

implicit none

integer irec, istep, nsites

character(nnam0) elstepname

logical stepnotdone, readnsites, readinitlstate, readfinalstate, readneighb

! Check for completeness of information:
stepnotdone = .false.
moreinfo = 'Missing information for step ' // trim(int2str(istep)) // ' (' // &
            trim(elstepname) // '):'
if (.not.readnsites) then
    stepnotdone = .true.
    moreinfo = trim(moreinfo) // ' number of sites,'
endif
if (.not.readinitlstate) then
    stepnotdone = .true.
    moreinfo = trim(moreinfo) // ' initial state,'
endif
if (.not.readfinalstate) then
    stepnotdone = .true.
    moreinfo = trim(moreinfo) // ' final state,'
endif
if ((.not.readneighb) .and. nsites > 1) then
    stepnotdone = .true.
    moreinfo = trim(moreinfo) // ' neighboring structure,'
endif
moreinfo = trim(moreinfo) // ' has(have) not been read. ' // &
           'Occurred while parsing line ' // trim(int2str(irec)) // '.'

if (stepnotdone) then
    call error(3024)
endif

return

end subroutine check_level1_elemstep_completeness

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_level2_elemstep_completeness(irec,istep,elstepname,ivariant,varname, &
                   isreversible,readpreexp,readactiveng,readperatio) !,readdeltah

implicit none

integer irec, istep, ivariant

logical stepnotdone, readpreexp, readactiveng, readperatio
logical isreversible
character(nnam0) varname, elstepname

! Check for completeness of information:
stepnotdone = .false.
moreinfo = 'Missing information for step ' // trim(int2str(istep)) // ' (' // &
            trim(elstepname) // ')'
if (ivariant > 0) then
    moreinfo = trim(moreinfo) // ' variant ' // trim(int2str(ivariant)) // ' (' // &
            trim(varname) // ')'
endif
moreinfo = trim(moreinfo) // ':'

if (.not.readpreexp) then
    stepnotdone = .true.
    moreinfo = trim(moreinfo) // ' pre-exponential,'
endif
if (isreversible .and. .not.readperatio) then
    stepnotdone = .true.
    moreinfo = trim(moreinfo) // ' pre-exponential fwd/rev ratio for reversible step,'
endif
if (.not.readactiveng) then
    stepnotdone = .true.
    moreinfo = trim(moreinfo) // ' activation energy,'
endif
moreinfo = trim(moreinfo) // ' has(have) not been read. ' // &
           'Occurred while parsing line ' // trim(int2str(irec)) // '.'

if (stepnotdone) then
    call error(3024)
endif

end subroutine check_level2_elemstep_completeness
                   
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

real(8) function tempfun(tempini,tempramp,timereq,fcntype)

implicit none

integer, intent(in) :: fcntype
real(8), intent(in) :: tempini
real(8), intent(in) :: tempramp
real(8), intent(in) :: timereq


select case (fcntype)
    
    case (1) ! linear ramp (which is the only supported dependence currently)

        tempfun = max(tiny(1.d0),tempini+tempramp*timereq)
    
    case default
    
        write(iwrite,*) 'ERROR: Temperature dependence can only be linear'
        stop

end select

end function tempfun

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

real(8) function timeoftemp0(tempini,tempramp,fcntype)

implicit none

integer, intent(in) :: fcntype
real(8), intent(in) :: tempini
real(8), intent(in) :: tempramp

! Function that gives the time when Temperature = 0 K. Useful in simulated 
! annealing  calculations.

select case (fcntype)
    
    case (1) ! linear ramp (which is the only supported dependence currently)
        
    if (tempramp > 0.d0) then
        timeoftemp0 = huge(1.d0)
    else
        timeoftemp0 = -tempini/tempramp
    endif
    
    case default
    
        write(iwrite,*) 'ERROR: Temperature dependence can only be linear'
        stop

end select

end function timeoftemp0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine calculate_preex_factor(jelemstep,tempcur,preexpfac)

use simulation_setup_module, only: pres, maxtime, temp, tpdsim, tramp,         &
                                   gasmolfracs
implicit none

integer, intent(in) :: jelemstep
real(8), intent(in) :: tempcur
real(8), intent(out) :: preexpfac
real(8) tempfin, tempmin, tempmax

integer i


if (preexpisconst(jelemstep)) then

    preexpfac = preexp(jelemstep,0)

else
    
    tempfin = tempfun(temp,tramp,maxtime,1)
    tempmin = min(temp,tempfin)
    tempmax = max(temp,tempfin)
    
    if (tempcur < tempmin) then
        preexpfac = dexp( -( &
                        preexp(jelemstep,1)*dlog(tempmin) + &
                        preexp(jelemstep,2)*1.d0/tempmin + &
                        preexp(jelemstep,3) + &
                        preexp(jelemstep,4)*tempmin + &
                        preexp(jelemstep,5)*tempmin**2 + &
                        preexp(jelemstep,6)*tempmin**3 + &
                        preexp(jelemstep,7)*tempmin**4 &
                        ) )
    elseif (tempcur > tempmax) then
        preexpfac = dexp( -( &
                        preexp(jelemstep,1)*dlog(tempmax) + &
                        preexp(jelemstep,2)*1.d0/tempmax + &
                        preexp(jelemstep,3) + &
                        preexp(jelemstep,4)*tempmax + &
                        preexp(jelemstep,5)*tempmax**2 + &
                        preexp(jelemstep,6)*tempmax**3 + &
                        preexp(jelemstep,7)*tempmax**4 &
                        ) )
    
    else    
        preexpfac = dexp( -( &
                        preexp(jelemstep,1)*dlog(tempcur) + &
                        preexp(jelemstep,2)*1.d0/tempcur + &
                        preexp(jelemstep,3) + &
                        preexp(jelemstep,4)*tempcur + &
                        preexp(jelemstep,5)*tempcur**2 + &
                        preexp(jelemstep,6)*tempcur**3 + &
                        preexp(jelemstep,7)*tempcur**4 &
                        ) )
    endif
    
    if (.not.tpdsim) then
        preexp(jelemstep,0) = preexpfac
        preexpisconst(jelemstep) = .true.
    endif
    
endif

do i = 1,elemstepgases(jelemstep,0)
    if (elemstepgases(jelemstep,2*i) < 0) then
        preexpfac = preexpfac*pres*gasmolfracs(elemstepgases(jelemstep,2*i-1))
    endif
enddo

return

end subroutine calculate_preex_factor

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module mechanism_setup_module
