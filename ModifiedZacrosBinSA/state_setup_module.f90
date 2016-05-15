! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module state_setup_module

use constants_module
use error_module
use parser_module
use graph_functions_module

use simulation_setup_module
use lattice_setup_module

implicit none

integer nadsorb

integer nentrmult, dentmxcoord
integer, allocatable :: specseedmult(:)
integer, allocatable :: nseedmult(:)
integer, allocatable :: specseedmultneigh(:,:,:)
integer, allocatable :: specseedmultstypes(:,:)

integer nentronsites
integer, allocatable :: speciesonsites(:)
integer, allocatable :: specseedonsites(:,:)

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine read_state_setup()

implicit none

integer io, irec, i, j, nwords, q, setvbuf3f
integer iw(maxwords2)
integer ientrmult, ientronsites

character(lengthrecinp) recinput

logical statefex, genoutop, inistateopen, entryopen

! Parses the simulation setup file state_input.dat

! Level-1 keywords allowed in this subroutine:
!    initial_state
!    end_initial_state
!    seed_multiple
!    seed_on_sites

! The program first performs a parsing of the entire input to check whether the 
! blocks are opened and closed properly and collect information that will be used 
! for determining the bounds of the arrays used for initial state specification.

! Check if the initial state specification file exists
! This is not a mandatory file, so if it doesn't exist the program exits and 
! the simulation initializes with an empty lattice
statefex = .false.
nentrmult = 0
nentronsites = 0
inquire(file=cstatefname,exist=statefex)
if (.not.statefex) then
    return
endif

inquire(iwrite,opened=genoutop)
if (.not.genoutop) then
    open(unit=iwrite,file=trim(cgenoutfname),status='unknown',position='append')
    q = setvbuf3f(iwrite,1,0) ! set line-buffered behavior
endif

write(iwrite,'(/,a)') 'Initial state setup:'
write(iwrite,'(a)')   '~~~~~~~~~~~~~~~~~~~~'

open(unit=istateread,file=trim(cstatefname),status='old')

! Perform an initial pass to check for proper syntax and opening and closing
! of blocks and find the following:
!    number of seed_multiple entries: nentrmult
!    the max coordination number of dentates: dentmxcoord
!    number of seed_on_sites entries: nentronsites

ientrmult = 0
ientronsites = 0
inistateopen = .false.
entryopen = .false.
! nentrmult and nentronsites have already been set to zero
dentmxcoord = 0

irec = 0
io = 0

call getrecord(istateread,recinput,irec,io)

do while (io >= 0)
    
    call break_words(recinput,' ',remchar,nwords,iw)
    
    if (nwords == 0) then
        call getrecord(istateread,recinput,irec,io)
        cycle
    endif        
        
    if (striccompare(recinput(iw(1):iw(2)),'initial_state')) then
    
        if (inistateopen) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(cstatefname) // '.'
            call error(100)
        endif
    
        inistateopen = .true.

    elseif (striccompare(recinput(iw(1):iw(2)),'end_initial_state')) then
    
        inistateopen = .false.
        exit

    elseif ( striccompare(recinput(iw(1):iw(2)),'seed_multiple') .or. &
             striccompare(recinput(iw(1):iw(2)),'seed_on_sites') ) then
        
        if (.not. inistateopen) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(4001)
        endif
        
        call parse_specseeding(.true.,recinput,nwords,irec,io,ientrmult,ientronsites,iw,entryopen)
                        
    else
        
        moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
        call error(4100)

    endif

    call getrecord(istateread,recinput,irec,io)

enddo

if (entryopen) then
    moreinfo = 'Elementary step information not completed in line ' // trim(int2str(irec)) // '.'
    call error(4013)
endif
if (inistateopen) then
    call warning(4100)
endif
if (nentrmult + nentronsites == 0) then
    moreinfo = 'File does not appear to contain any initial state information.'
    call error(4001)
endif

! specseedmult(i1) gives the species that will be seeded according to the 
! multiple-seeding instruction i1
allocate(specseedmult(nentrmult))

! specseedmult(i1) gives the number of adsorbate molecules that will be seeded 
! according to the multiple-seeding instruction i1
allocate(nseedmult(nentrmult))

! specseedmultneigh(i1,i2,i3) gives the i3^th neighbor of the i2^th dentate
!   in the i1^th seed_multiple entry. If i3 = 0 it gives the total number of
!   neighbors for dentate i2. Note that this structure is the subgraph
!   that will be searched for in the total lattice.
allocate(specseedmultneigh(nentrmult,maxdent,0:dentmxcoord))

! specseedmultstypes(i1,i2) gives the site type of the i2^th dentate
! for the species appearing in the multiple-seeding instruction i1
allocate(specseedmultstypes(nentrmult,maxdent))

! speciesonsites(i1) gives the species that will be seeded according to the 
! site-explicit-seeding instruction i1
allocate(speciesonsites(nentronsites))

! specseedonsites(i1,i2) gives the site type of the i2^th dentate
! for the species appearing in the site-explicit-seeding instruction i1
allocate(specseedonsites(nentronsites,maxdent))

! Initialization of variables
do i = 1,nentrmult
    specseedmult(i) = 0
    nseedmult(i) = 0
    do j = 1,maxdent
        specseedmultneigh(i,j,0) = 0
        specseedmultstypes(i,j) = -1
    enddo
enddo

do i = 1,nentronsites
    speciesonsites(i) = 0
    do j = 1,maxdent
        specseedonsites(i,j) = -1
    enddo
enddo

! The first pass has not found any problems in the syntax. Begin parsing the state setup instructions

ientrmult = 0
ientronsites = 0
entryopen = .false.

rewind(istateread)
irec = 0
io = 0

do while (io >= 0)

    call getrecord(istateread,recinput,irec,io)
    
    call break_words(recinput,' ',remchar,nwords,iw)
    
    if (nwords == 0) cycle
            
    if (striccompare(recinput(iw(1):iw(2)),'end_initial_state')) then

        exit
        
    elseif (striccompare(recinput(iw(1):iw(2)),'seed_multiple') .or. &
            striccompare(recinput(iw(1):iw(2)),'seed_on_sites')) then
    
        call parse_specseeding(.false.,recinput,nwords,irec,io,ientrmult,ientronsites,iw,entryopen)
                
    endif
        
enddo

close(istateread)

if (nentrmult > 0) then
    write(iwrite,'(/,a)') trim(int2str(nentrmult)) // ' "multiple" seeding instructions for:'
    do i = 1,nentrmult
        write(iwrite,'(/,a)') '  ' // trim(int2str(nseedmult(i))) //  &
                              ' adparticles of species ' // trim(surfspecsnames(specseedmult(i)))
    enddo
endif

if (nentronsites > 0) then
    write(iwrite,'(/,a)',advance='no') trim(int2str(nentronsites)) // ' "on-sites" seeding instructions:'
    do i = 1,nentronsites
        write(iwrite,'(/,a)',advance='no') '  ' // trim(surfspecsnames(speciesonsites(i))) // ' adparticle to be seeded on site(s) '

        do j=1,surfspecsdent(speciesonsites(i))
            write(iwrite,'(a,1x)',advance='no') trim(int2str(specseedonsites(i,j)))
        enddo
        
    enddo
endif

write(iwrite,'(a)') '  '
write(iwrite,'(/,a)') 'Finished reading initial state input.'

return

end subroutine read_state_setup

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine parse_specseeding(firstpass,recinput,nwords,irec,io,ientrmult,ientronsites,iw,entryopen)  

implicit none

integer i, irec, io, ientrmult, ientronsites, ispec, istype, j, j1, j2, k, m, nwords, nwords2
integer nv(25), tmpv(2), iw(maxwords2), iw2(maxwords2)

character(lengthrecinp) recinput, tmpstr

logical firstpass, entryopen, readneighb, readstypes

! Parses a single seeding instruction in state_input.dat

! Level-1 keywords allowed in this subroutine:
!  -Normal keywords valid in a seed_multiple block
!    neighboring
!    site_types
!    end_seed_multiple
!  -The seed_on_sites instruction is a single line
!  -The following keywords are recognized as invalid in this context and produce and error
!    seed_multiple
!    seed_on_sites

! Depending on the value of the logical variable firstpass the subroutine can be run in one of 
! two "modes": a first pass mode and a full parsing mode. The former just checks for consistency
! of input and finds the max bounds needed for the arrays that encode the mechanism. The latter
! populates these arrays. The user MUST run a first pass of the mechanism before the full parsing.


ispec = findinstrv(surfspecsnames,0,1,nsurfspecs,recinput(iw(3):iw(4)))
if (firstpass .and. ispec < 0) then
    moreinfo = 'Encountered unknown surface species "' // recinput(iw(3):iw(4)) // '" in line ' // trim(int2str(irec)) // '.'
    call error(4002)    
endif

if (striccompare(recinput(iw(1):iw(2)),'seed_on_sites')) then

    if (str_check_expression(recinput,nwords,iw,'2AA/' // trim(int2str(surfspecsdent(ispec))) // 'I4',.true.)) then

        ientronsites = ientronsites + 1
        if (firstpass) then
            nentronsites = ientronsites
            return
        endif

        speciesonsites(ientronsites) = ispec
        do i = 1,surfspecsdent(ispec)
            k = str2int(recinput(iw(2*i+3):iw(2*i+4)))
            if (1 <= k .and. k <= nsites) then
                specseedonsites(ientronsites,i) = k
            else
                moreinfo = 'Site number  ' // recinput(iw(2*i+3):iw(2*i+4)) // ' in line ' // trim(int2str(irec)) // ' is ' // &
                           'out of bounds 1:' // int2str(nsites) // '.'
                call error(4003)
            endif
        enddo
        
        return

    else

        moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
        call error(4004)

    endif

elseif (striccompare(recinput(iw(1):iw(2)),'seed_multiple')) then

    if (str_check_expression(recinput,nwords,iw,'2AA/I4',.true.)) then
        if (str2int(recinput(iw(5):iw(6))) <= 0) then
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(4005)            
        endif
        ientrmult = ientrmult + 1
        readneighb = .false.
        readstypes = .false.
        if (firstpass) then 
            nentrmult = nentrmult + 1
        else
            specseedmult(ientrmult) = ispec
            nseedmult(ientrmult) = str2int(recinput(iw(5):iw(6)))
        endif
    else
        moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
        call error(4005)
    endif

    entryopen = .true.
    do while (io >= 0)

        call getrecord(istateread,recinput,irec,io)
        
        call break_words(recinput,' ',remchar,nwords,iw)

        if (nwords == 0) then
            cycle
        endif
        
        if (striccompare(recinput(iw(1):iw(2)),'neighboring')) then
            
            ! For a first pass...
            if (firstpass) then
            
                if (surfspecsdent(ispec) == 1) then
                    moreinfo = 'Encountered "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(4006)
                endif
                
                nv = 0 ! remark: nv is a vector
                
                do j = 1,nwords-1
                    
                    tmpstr = recinput(iw(2*j+1):iw(2*j+2))
                    call break_words(trim(tmpstr),'-',remchar,nwords2,iw2)
                    
                    if (.not.str_check_expression(tmpstr,nwords2,iw2,'2I4',.true.)) then
                        moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                        call error(4007)
                    endif                    
                    
                    j1 = str2int(tmpstr(iw2(1):iw2(2)))
                    j2 = str2int(tmpstr(iw2(3):iw2(4)))
                    
                    if (j1 == j2) then
                        moreinfo = 'Encountered "' // recinput(iw2(1):iw2(2)) // '" in line ' // trim(int2str(irec)) // '.'
                        call error(4008)
                    endif
                    
                    tmpv = (/j1,j2/)
                    do k = 1,2                                    
                        if (tmpv(k) > surfspecsdent(ispec) .or. tmpv(k) < 1) then
                            moreinfo = 'Site number ' // trim(int2str(tmpv(k))) // &
                                       ' referenced in neighboring input in line ' // trim(int2str(irec)) // &
                                       ' is out of range 1:' // trim(int2str(surfspecsdent(ispec))) // &
                                       ' for this elementary step.'
                            call error(4007)
                        endif
                        
                    enddo
                    
                    nv(j1) = nv(j1) + 1
                    nv(j2) = nv(j2) + 1
                    dentmxcoord = max(dentmxcoord,nv(j1),nv(j2))
                    
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
                
                specseedmultneigh(ientrmult,j1,0) = specseedmultneigh(ientrmult,j1,0) + 1
                specseedmultneigh(ientrmult,j1,specseedmultneigh(ientrmult,j1,0)) = j2

                specseedmultneigh(ientrmult,j2,0) = specseedmultneigh(ientrmult,j2,0) + 1
                specseedmultneigh(ientrmult,j2,specseedmultneigh(ientrmult,j2,0)) = j1
                
                ! Check if the j1-j2 pair has already been encountered in the neighboring list
                do k = 1,specseedmultneigh(ientrmult,j1,0)
                    do m = k+1,specseedmultneigh(ientrmult,j1,0)
                        if (specseedmultneigh(ientrmult,j1,k) == specseedmultneigh(ientrmult,j1,m)) then
                            moreinfo = 'Pair ' // trim(int2str(j1)) // '-' // trim(int2str(j2)) // &
                                       ' is referenced twice in neighboring input in line ' // trim(int2str(irec)) // '.'
                            call error(4009)
                        endif
                    enddo
                enddo
            
            enddo
            
            readneighb = .true.
       
        elseif (striccompare(recinput(iw(1):iw(2)),'site_types')) then
            
            if (firstpass) then
                
                if (nwords-1 == surfspecsdent(ispec)) then
                    readstypes = .true.
                    cycle
                else
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(4010)
                endif        
                
            endif
            
            do j = 1,nwords-1
                
                if (str_is_integer(recinput(iw(2*j+1):iw(2*j+2)),4)) then
                    istype = str2int(recinput(iw(2*j+1):iw(2*j+2)))
                    if (istype > nsitetypes) then
                        moreinfo = 'Encountered site type ' // trim(int2str(istype)) // ' for dentate ' // trim(int2str(j)) // &
                                   ' in line ' // trim(int2str(irec)) // '.'
                        call error(4011)
                    endif
                else
                    istype = findinstrv(sitetypenames,1,1,nsitetypes,recinput(iw(2*j+1):iw(2*j+2)))                                
                    if (istype < 0) then
                        moreinfo = 'Encountered site type name "' // recinput(iw(2*j+1):iw(2*j+2)) // '" for site ' // trim(int2str(j)) // &
                                   ' in line ' // trim(int2str(irec)) // '.'
                        call error(4011)
                    endif
                endif
                
                specseedmultstypes(ientrmult,j) = istype
                                                            
            enddo
            
            readstypes = .true.
            
        elseif (striccompare(recinput(iw(1):iw(2)),'end_seed_multiple')) then
            
            entryopen = .false.
            
            exit
        
        elseif (striccompare(recinput(iw(1):iw(2)),'seed_multiple') .or. &
                striccompare(recinput(iw(1):iw(2)),'seed_on_sites')) then

            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(4012)

        else
            moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
            call error(4100)
        endif
    
    enddo
    
endif

if ((.not.readneighb) .and. surfspecsdent(ispec) > 1) then ! dentate neighboring not specified
    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
    call error(4101)
endif

return

end subroutine parse_specseeding

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module state_setup_module
