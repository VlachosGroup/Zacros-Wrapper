! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module lattice_setup_module

use constants_module, only: nnam0
implicit none

integer nsites, nsitetypes, maxcoord ! in lattice
integer ncells, ncellsites ! in lattice cell
integer ncellsrep(2)
integer, allocatable :: cellsitetype(:) ! site type in cell
integer, allocatable :: cellsiteneighsf(:,:) ! neighboring sites in self cell 
integer, allocatable :: cellsiteneighnn(:,:) ! neighboring sites in north cell 
integer, allocatable :: cellsiteneighne(:,:) ! neighboring sites in north-east cell 
integer, allocatable :: cellsiteneighee(:,:) ! neighboring sites in east cell 
integer, allocatable :: cellsiteneighse(:,:) ! neighboring sites in south-east cell 
integer, allocatable :: cellsiteneighss(:,:) ! neighboring sites in south cell 
integer, allocatable :: cellsiteneighsw(:,:) ! neighboring sites in south-west cell 
integer, allocatable :: cellsiteneighww(:,:) ! neighboring sites in west cell 
integer, allocatable :: cellsiteneighnw(:,:) ! neighboring sites in north-west cell 
integer, allocatable :: sitetype(:) ! site type in lattice
integer, allocatable :: siteneighb1(:,:) ! 1st nearest neighbors in lattice

real(8) vcell(2,2) ! This gives the unit cell
! vcell = [ v1x v2x
!           v1y x2y]

real(8) v1box(2) ! These give the cell-vectors defining the simulation box 
real(8) v2box(2)

real(8), allocatable :: xfcellsite(:), yfcellsite(:) ! fractional coordinates of sites in a single cell
real(8), allocatable :: xsite(:), ysite(:) ! coordinates of each site in lattice

character(nnam0), allocatable :: sitetypenames(:)

contains

  subroutine cleanup_lsm
    ! Cleanup module globals.

    if(allocated(cellsitetype))    deallocate(cellsitetype)
    if(allocated(cellsiteneighsf)) deallocate(cellsiteneighsf)
    if(allocated(cellsiteneighnn)) deallocate(cellsiteneighnn)
    if(allocated(cellsiteneighne)) deallocate(cellsiteneighne)
    if(allocated(cellsiteneighee)) deallocate(cellsiteneighee)
    if(allocated(cellsiteneighse)) deallocate(cellsiteneighse)
    if(allocated(cellsiteneighss)) deallocate(cellsiteneighss)
    if(allocated(cellsiteneighsw)) deallocate(cellsiteneighsw)
    if(allocated(cellsiteneighww)) deallocate(cellsiteneighww)
    if(allocated(cellsiteneighnw)) deallocate(cellsiteneighnw)
    if(allocated(sitetype))        deallocate(sitetype)
    if(allocated(siteneighb1))     deallocate(siteneighb1)
    if(allocated(xfcellsite))      deallocate(xfcellsite)
    if(allocated(yfcellsite))      deallocate(yfcellsite)
    if(allocated(xsite))           deallocate(xsite)
    if(allocated(ysite))           deallocate(ysite)
    if(allocated(sitetypenames))   deallocate(sitetypenames)

  end subroutine cleanup_lsm

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine read_lattice_setup()

use constants_module, only: clatfname, cgenoutfname, lengthrecinp, ilatwrite,  &
                            maxwords2, ilatread, remchar, iwrite, clatoutfname
use error_module, only: moreinfo, error, warning
use parser_module, only: int2str, striccompare, dbl2str, getrecord, break_words
implicit none

integer i, j, io, irec, nwords, setvbuf3f, q
integer iw(maxwords2)

real(8) aunit

character(lengthrecinp) recinput

logical genoutop, readlat, periodiclat, latfex

! Parses the lattice setup file lattice_input.dat

! Level-1 keywords allowed in this subroutine:
!    lattice

! The flow of the program is directed to one of two subroutines
! construct_periodic_lattice or read_explicit_lattice:
!
! - The former constructs a periodic lattice from unit cell data. These can be
!   either default data (for triangular, rectangular or hexagonal lattices with
!   a single site type) or data read from syntaxed input in the lattice_input.dat
!   file. The two cases are differentiated by the value of the first argument 
!   of the subroutine ( = 0 for default unit cell, = 1 for supplied data)
!
! - The latter constructs a lattice from data read explicitly from the 
!   lattice_input.dat file.

! Check if the lattice specification file exists
latfex = .false.
inquire(file=clatfname, exist=latfex)
if (.not.latfex) then
    moreinfo = 'Lattice specification file ' // trim(clatfname) // ' does not exist.'
    call error(1)
endif

inquire(iwrite,opened=genoutop)
if (.not.genoutop) then
    open(unit=iwrite,file=trim(cgenoutfname),status='unknown',position='append')
    q = setvbuf3f(iwrite,1,0) ! set line-buffered behavior
endif

write(iwrite,'(/,a)') 'Lattice setup:'
write(iwrite,'(a)')   '~~~~~~~~~~~~~~'

open(unit=ilatread,file=trim(clatfname),status='old')

moreinfo = ''
readlat = .false.
periodiclat = .true.

irec = 0
io = 0

call getrecord(ilatread,recinput,irec,io)

do while (io >= 0)

    call break_words(recinput,' ',remchar,nwords,iw)

    if (nwords == 0) then
        call getrecord(ilatread,recinput,irec,io)
        cycle
    
    elseif (striccompare(recinput(iw(1):iw(2)),'lattice')) then
    
        if (nwords /= 2) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(2001)                
        else
            
            if (striccompare(recinput(iw(3):iw(4)),'default_choice')) then                
                write(iwrite,'(/,a)') '    Using one of the default lattice structures.'
                call construct_periodic_lattice(0,irec,io)
            
            elseif (striccompare(recinput(iw(3):iw(4)),'periodic_cell')) then
                write(iwrite,'(/,a)') '    Constructing a periodic lattice from unit cell data.'
                call construct_periodic_lattice(1,irec,io)
            
            elseif (striccompare(recinput(iw(3):iw(4)),'explicit')) then
                write(iwrite,'(/,a)') '    Explicit lattice structure information provided.'            
                call read_explicit_lattice(irec,io,periodiclat)
                
            else
                moreinfo = 'Unknown lattice directive ' // recinput(iw(3):iw(4)) // &
                           ' in line ' // trim(int2str(irec)) // '.'
                call error(2001)                

            endif
                
            readlat = .true.
            exit

        endif
    
    else ! lattice input must start with the keyword lattice, any other expression is invalid at this point
        moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
        call error(2001)
    
    endif
    
    call getrecord(ilatread,recinput,irec,io)

enddo

close(ilatread)

if (io < 0) then ! finish keyword not encountered
    call warning(2100)
endif

if (.not.readlat) then
    moreinfo = 'File does not appear to contain any lattice structure information.'
    call error(2001)
endif

v1box = (/0.d0,0.d0/)
v2box = (/0.d0,0.d0/)
if (periodiclat) then
    aunit = dabs(vcell(1,1)*vcell(2,2) - vcell(1,2)*vcell(2,1))
    v1box = ncellsrep(1)*vcell(:,1)
    v2box = ncellsrep(2)*vcell(:,2)
    write(iwrite,'(/,a)') '    Unit cell area: ' // trim(dbl2str(aunit))
    write(iwrite,'(/,a)') '    Lattice cell: ' // trim(int2str(ncellsrep(1))) // ' x ' // trim(int2str(ncellsrep(2)))
    write(iwrite,'(/,a)') '    Surface area: ' // trim(dbl2str(aunit*ncellsrep(1)*ncellsrep(2)))
endif

write(iwrite,'(/,a)') '    Number of lattice sites: ' // trim(int2str(nsites))
write(iwrite,'(/,a)') '    Number of site types: ' // trim(int2str(nsitetypes))

if (.not. allocated(sitetypenames)) then
    write(iwrite,'(/,a)',advance='no') '    Site type names not defined, automatic strings assigned.'
    allocate(sitetypenames(nsitetypes))
    do i = 1,nsitetypes
        sitetypenames(i) = 'StTp' // trim(int2str(i))
    enddo
endif

write(iwrite,'(/,a)') '    Site type names and number of sites of that type:'
do i = 1,nsitetypes
    write(iwrite,'(6x,a)') trim(sitetypenames(i)) // ' (' // trim(int2str(sum(sitetype,mask = sitetype == i)/i)) // ')'
enddo

  write(iwrite,'(/,a)') '    Maximum coordination number: ' // trim(int2str(maxcoord))
  write(iwrite,'(/,a)') '    Lattice structure written to ' // trim(clatoutfname)
write(iwrite,'(/,a)') 'Finished reading lattice input.'

open(unit=ilatwrite,file=clatoutfname,status='unknown')

write(ilatwrite,'(i10,1x,2ES32.16,'//int2str(maxcoord+2)//'i10)') &
       0,v1box(1),v1box(2),0,0,(0,j=1,maxcoord)
write(ilatwrite,'(i10,1x,2ES32.16,'//int2str(maxcoord+2)//'i10)') &
       0,v2box(1),v2box(2),0,0,(0,j=1,maxcoord)

do i = 1,nsites
    write(ilatwrite,'(i10,1x,2ES32.16,'//int2str(maxcoord+2)//'i10)') &
           i,xsite(i),ysite(i),sitetype(i),siteneighb1(i,0), &
           (siteneighb1(i,j),j=1,siteneighb1(i,0)), &
           (0,j=1,maxcoord-siteneighb1(i,0))
enddo

close(ilatwrite)

return

end subroutine read_lattice_setup

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine construct_periodic_lattice(indic,irec,io)

implicit none

integer i, irec, io, j, k, m, cntr, sntr, nntr, indic
integer icellneigh, jcellneigh, neighc

real(8) cellpos(2)

! Constucts a periodic lattice from unit cell data. The latter are supplied 
! by either the default_lattices subroutine (which encompasses the three planar
! regular lattices, triangular, rectangular and hexagonal, with single site types)
! or the read_periodic_lattice_cell subroutine in which the unit cell data is 
! parsed from the lattice_input.dat file.

! Note that this subroutine does not read anything from lattice_input.dat itself;
! it calls one of the aforementioned subroutines for this purpose.

select case (indic)
    case (0)
        call default_lattices(irec,io)
    case (1)
        call read_periodic_lattice_cell(irec,io)
end select

ncells = ncellsrep(1)*ncellsrep(2)
nsites = ncells*ncellsites

allocate(xsite(nsites))
allocate(ysite(nsites))
allocate(sitetype(nsites))
allocate(siteneighb1(nsites,0:maxcoord))

!write(*,*) ((cellsiteneighsf(i,j),j=0,ncellsites),i=1,ncellsites)
!write(*,*) ' '
!write(*,*) ((cellsiteneighnn(i,j),j=0,ncellsites),i=1,ncellsites)
!write(*,*) ' '
!write(*,*) ((cellsiteneighne(i,j),j=0,ncellsites),i=1,ncellsites)
!write(*,*) ' '
!write(*,*) ((cellsiteneighee(i,j),j=0,ncellsites),i=1,ncellsites)
!write(*,*) ' '
!write(*,*) ((cellsiteneighse(i,j),j=0,ncellsites),i=1,ncellsites)
!write(*,*) ' '
!write(*,*) ((cellsiteneighss(i,j),j=0,ncellsites),i=1,ncellsites)
!write(*,*) ' '
!write(*,*) ((cellsiteneighsw(i,j),j=0,ncellsites),i=1,ncellsites)
!write(*,*) ' '
!write(*,*) ((cellsiteneighww(i,j),j=0,ncellsites),i=1,ncellsites)
!write(*,*) ' '
!write(*,*) ((cellsiteneighnw(i,j),j=0,ncellsites),i=1,ncellsites)
!write(*,*) ' '
!write(*,*) nsites,ncellsites,ncellsrep(1),ncellsrep(2)

do i = 1,ncellsrep(1)
    do j = 1,ncellsrep(2)
        
        cntr = (i-1)*ncellsrep(2) + j ! cell counter
        cellpos = (i-1)*vcell(:,1) + (j-1)*vcell(:,2) ! cell position
        
        do k = 1,ncellsites
            
            sntr = ncellsites*(cntr-1) + k
            
            ! x-y coordinates of the site
            xsite(sntr) = xfcellsite(k)*vcell(1,1) + yfcellsite(k)*vcell(1,2) + CellPos(1)
            ysite(sntr) = xfcellsite(k)*vcell(2,1) + yfcellsite(k)*vcell(2,2) + CellPos(2)

            ! site type
            sitetype(sntr) = cellsitetype(k)
            
            ! Site neighbors
            nntr = 0

            ! in-cell neighbors
            do m = 1,cellsiteneighsf(k,0)
                nntr = nntr + 1
                siteneighb1(sntr,nntr) = ncellsites*(cntr-1) + cellsiteneighsf(k,m)
            enddo

            ! neighbors in north neighboring cell
            do m = 1,cellsiteneighnn(k,0)
                nntr = nntr + 1
                ! position of neighboring cell
                icellneigh = i
                jcellneigh = mod(j,ncellsrep(2)) + 1
                neighc = (icellneigh-1)*ncellsrep(2) + jcellneigh
                siteneighb1(sntr,nntr) = ncellsites*(neighc-1) + cellsiteneighnn(k,m)
            enddo
            
            ! neighbors in north-east neighboring cell
            do m = 1,cellsiteneighne(k,0)
                nntr = nntr + 1
                ! position of neighboring cell
                icellneigh = mod(i,ncellsrep(1)) + 1
                jcellneigh = mod(j,ncellsrep(2)) + 1
                neighc = (icellneigh-1)*ncellsrep(2) + jcellneigh
                siteneighb1(sntr,nntr) = ncellsites*(neighc-1) + cellsiteneighne(k,m)
            enddo
            
            ! neighbors in east neighboring cell
            do m = 1,cellsiteneighee(k,0)
                nntr = nntr + 1
                ! position of neighboring cell
                icellneigh = mod(i,ncellsrep(1)) + 1
                jcellneigh = j
                neighc = (icellneigh-1)*ncellsrep(2) + jcellneigh
                siteneighb1(sntr,nntr) = ncellsites*(neighc-1) + cellsiteneighee(k,m)
            enddo

            ! neighbors in south-east neighboring cell
            do m = 1,cellsiteneighse(k,0)
                nntr = nntr + 1
                ! position of neighboring cell
                icellneigh = mod(i,ncellsrep(1)) + 1
                jcellneigh = mod(j-2+ncellsrep(2),ncellsrep(2)) + 1
                neighc = (icellneigh-1)*ncellsrep(2) + jcellneigh
                siteneighb1(sntr,nntr) = ncellsites*(neighc-1) + cellsiteneighse(k,m)
            enddo
           
            ! neighbors in south neighboring cell
            do m = 1,cellsiteneighss(k,0)
                nntr = nntr + 1
                ! position of neighboring cell
                icellneigh = i
                jcellneigh = mod(j-2+ncellsrep(2),ncellsrep(2)) + 1
                neighc = (icellneigh-1)*ncellsrep(2) + jcellneigh
                siteneighb1(sntr,nntr) = ncellsites*(neighc-1) + cellsiteneighss(k,m)
            enddo

            ! neighbors in south-west neighboring cell
            do m = 1,cellsiteneighsw(k,0)
                nntr = nntr + 1
                ! position of neighboring cell
                icellneigh = mod(i-2+ncellsrep(1),ncellsrep(1)) + 1
                jcellneigh = mod(j-2+ncellsrep(2),ncellsrep(2)) + 1
                neighc = (icellneigh-1)*ncellsrep(2) + jcellneigh
                siteneighb1(sntr,nntr) = ncellsites*(neighc-1) + cellsiteneighsw(k,m)
            enddo
            
            ! neighbors in west neighboring cell
            do m = 1,cellsiteneighww(k,0)
                nntr = nntr + 1
                ! position of neighboring cell
                icellneigh = mod(i-2+ncellsrep(1),ncellsrep(1)) + 1
                jcellneigh = j
                neighc = (icellneigh-1)*ncellsrep(2) + jcellneigh
                siteneighb1(sntr,nntr) = ncellsites*(neighc-1) + cellsiteneighww(k,m)
            enddo
            
            ! neighbors in north-west neighboring cell
            do m = 1,cellsiteneighnw(k,0)
                nntr = nntr + 1
                ! position of neighboring cell
                icellneigh = mod(i-2+ncellsrep(1),ncellsrep(1)) + 1
                jcellneigh = mod(j,ncellsrep(2)) + 1
                neighc = (icellneigh-1)*ncellsrep(2) + jcellneigh
                siteneighb1(sntr,nntr) = ncellsites*(neighc-1) + cellsiteneighnw(k,m)
            enddo
            
            siteneighb1(sntr,0) = nntr
            
        enddo
        
    enddo
enddo

! Check self-consistency of lattice structure
call check_lattice()

return

end subroutine construct_periodic_lattice

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine default_lattices(irec,io)

use constants_module, only: maxwords2, lengthrecinp, remchar, ilatread,        &
                            clatfname
use parser_module, only: striccompare, str2int, int2str, str2dbl,              &
                         str_check_expression, break_words, getrecord
use error_module, only: error, moreinfo
implicit none

integer nwords, i, io, irec, j, jneigh, lattind, icoord
integer iw(maxwords2)

real(8) alfa

logical readlat

character(lengthrecinp) recinput

! Parses a default lattice specification in lattice_input.dat

! Level-1 keywords allowed:
!    triangular_periodic
!    rectangular_periodic
!    hexagonal_periodic
!    end_lattice

readlat = .false.

io = 0

call getrecord(ilatread,recinput,irec,io)

do while (io >= 0)

    call break_words(recinput,' ',remchar,nwords,iw)

    if (nwords == 0) then

        call getrecord(ilatread,recinput,irec,io)
        cycle
        
    elseif (striccompare(recinput(iw(1):iw(2)),'end_lattice')) then

        exit

    elseif (striccompare(recinput(iw(1):iw(2)),'lattice')) then

        moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
        'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
        call error(100)

    elseif (str_check_expression(recinput,nwords,iw,'AA/R8/2I4',.true.)) then
        
        if (readlat) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(2002)
        endif
        
        if (striccompare(recinput(iw(1):iw(2)),'triangular_periodic')) then
            lattind = 3
        elseif (striccompare(recinput(iw(1):iw(2)),'rectangular_periodic')) then
            lattind = 4
        elseif (striccompare(recinput(iw(1):iw(2)),'hexagonal_periodic')) then
            lattind = 6
        else    
            moreinfo = 'Unknown default lattice keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
            call error(2003)                
        endif
        
        alfa = str2dbl(recinput(iw(3):iw(4)))
        ncellsrep = (/ str2int(recinput(iw(5):iw(6))), &
                       str2int(recinput(iw(7):iw(8))) /)                    

        readlat = .true.

    else    
        moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
        call error(2003)                
    endif
        
    call getrecord(ilatread,recinput,irec,io)

enddo

if (.not.readlat) then
    moreinfo = 'Directive for default lattice not encountered.'
    call error(2003)
endif

! Lattice selection: regular lattices with a single site type
nsitetypes = 1
select case (lattind)

    case(3)
        ncellsites = 4
        maxcoord = 3

    case(4)
        ncellsites = 1
        maxcoord = 4

    case(6)
        ncellsites = 2
        maxcoord = 6

endselect

allocate(xfcellsite(ncellsites))
allocate(yfcellsite(ncellsites))
allocate(cellsitetype(ncellsites))
allocate(cellsiteneighsf(ncellsites,0:ncellsites))
allocate(cellsiteneighnn(ncellsites,0:ncellsites))
allocate(cellsiteneighne(ncellsites,0:ncellsites))
allocate(cellsiteneighee(ncellsites,0:ncellsites))
allocate(cellsiteneighse(ncellsites,0:ncellsites))
allocate(cellsiteneighss(ncellsites,0:ncellsites))
allocate(cellsiteneighsw(ncellsites,0:ncellsites))
allocate(cellsiteneighww(ncellsites,0:ncellsites))
allocate(cellsiteneighnw(ncellsites,0:ncellsites))

select case (lattind)

    case(3)

        cellsitetype = (/ 1, 1, 1, 1 /)

        vcell = reshape((/ alfa*dsqrt(3.d0), 0.d0, 0.d0, 3.d0*alfa /), (/ 2, 2 /))

        xfcellsite = (/ 0.0d0, 1.0d0/2.0d0, 1.0d0/2.0d0, 0.0d0/)
        yfcellsite = (/ 0.0d0, 1.0d0/6.0d0, 1.0d0/2.0d0, 2.0d0/3.0d0 /)

        cellsiteneighsf(1,0) = 1
        cellsiteneighsf(1,1) = 2
        cellsiteneighsf(2,0) = 2
        cellsiteneighsf(2,1) = 1
        cellsiteneighsf(2,2) = 3
        cellsiteneighsf(3,0) = 2
        cellsiteneighsf(3,1) = 2
        cellsiteneighsf(3,2) = 4
        cellsiteneighsf(4,0) = 1
        cellsiteneighsf(4,1) = 3

        cellsiteneighnn(1,0) = 0
        cellsiteneighnn(2,0) = 0
        cellsiteneighnn(3,0) = 0
        cellsiteneighnn(4,0) = 1
        cellsiteneighnn(4,1) = 1

        cellsiteneighne(1,0) = 0
        cellsiteneighne(2,0) = 0
        cellsiteneighne(3,0) = 0
        cellsiteneighne(4,0) = 0

        cellsiteneighee(1,0) = 0
        cellsiteneighee(2,0) = 1
        cellsiteneighee(2,1) = 1
        cellsiteneighee(3,0) = 1
        cellsiteneighee(3,1) = 4
        cellsiteneighee(4,0) = 0

        cellsiteneighse(1,0) = 0
        cellsiteneighse(2,0) = 0
        cellsiteneighse(3,0) = 0
        cellsiteneighse(4,0) = 0

    case(4)

        cellsitetype = (/ 1 /)

        vcell = reshape((/ alfa, 0.d0, 0.d0, alfa /), (/ 2, 2 /))

        xfcellsite = (/ 0.0d0 /)
        yfcellsite = (/ 0.0d0 /)

        cellsiteneighsf(1,0) = 0

        cellsiteneighnn(1,0) = 1
        cellsiteneighnn(1,1) = 1

        cellsiteneighne(1,0) = 0

        cellsiteneighee(1,0) = 1
        cellsiteneighee(1,1) = 1

        cellsiteneighse(1,0) = 0

    case(6)

        cellsitetype = (/ 1, 1 /)

        vcell = reshape((/ alfa*dsqrt(3.d0), 0.d0, 0.d0, alfa /), (/ 2, 2 /))

        xfcellsite = (/ 0.0d0, 0.5d0 /)
        yfcellsite = (/ 0.0d0, 0.5d0 /)

        cellsiteneighsf(1,0) = 1
        cellsiteneighsf(1,1) = 2
        cellsiteneighsf(2,0) = 1
        cellsiteneighsf(2,1) = 1

        cellsiteneighnn(1,0) = 1
        cellsiteneighnn(1,1) = 1
        cellsiteneighnn(2,0) = 2
        cellsiteneighnn(2,1) = 1
        cellsiteneighnn(2,2) = 2

        cellsiteneighne(1,0) = 0
        cellsiteneighne(2,0) = 1
        cellsiteneighne(2,1) = 1

        cellsiteneighee(1,0) = 0
        cellsiteneighee(2,0) = 1
        cellsiteneighee(2,1) = 1

        cellsiteneighse(1,0) = 0
        cellsiteneighse(2,0) = 0

endselect

! Construct the neighboring lists for the south, south-west, west and north-west
! cells
do i = 1,ncellsites

    cellsiteneighss(i,0) = 0
    cellsiteneighsw(i,0) = 0
    cellsiteneighww(i,0) = 0
    cellsiteneighnw(i,0) = 0

enddo

do i = 1,ncellsites

    do j = 1,cellsiteneighnn(i,0)
        jneigh = cellsiteneighnn(i,j)
        cellsiteneighss(jneigh,0) = cellsiteneighss(jneigh,0) + 1
        cellsiteneighss(jneigh,cellsiteneighss(jneigh,0)) = i
    enddo

    do j = 1,cellsiteneighne(i,0)
        jneigh = cellsiteneighne(i,j)
        cellsiteneighsw(jneigh,0) = cellsiteneighsw(jneigh,0) + 1
        cellsiteneighsw(jneigh,cellsiteneighsw(jneigh,0)) = i
    enddo

    do j = 1,cellsiteneighee(i,0)
        jneigh = cellsiteneighee(i,j)
        cellsiteneighww(jneigh,0) = cellsiteneighww(jneigh,0) + 1
        cellsiteneighww(jneigh,cellsiteneighww(jneigh,0)) = i
    enddo

    do j = 1,cellsiteneighse(i,0)
        jneigh = cellsiteneighse(i,j)
        cellsiteneighnw(jneigh,0) = cellsiteneighnw(jneigh,0) + 1
        cellsiteneighnw(jneigh,cellsiteneighnw(jneigh,0)) = i
    enddo

enddo

! Find maximum coordination number

maxcoord = 0

do i = 1,ncellsites

    icoord = cellsiteneighsf(i,0) + &
        cellsiteneighnn(i,0) + cellsiteneighne(i,0) + cellsiteneighee(i,0) + cellsiteneighse(i,0) + &
        cellsiteneighss(i,0) + cellsiteneighsw(i,0) + cellsiteneighww(i,0) + cellsiteneighnw(i,0)

    if (icoord > maxcoord) then
        maxcoord = icoord
    endif

enddo

! Check self-consistency of unit cell
call check_unit_cell()

return

end subroutine default_lattices

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine read_periodic_lattice_cell(irec,io)

use constants_module, only: maxwords2, lengthrecinp, remchar, ilatread,        &
                            clatfname
use parser_module, only: striccompare, str_check_expression, str2dbl, str2int, &
                         int2str, findinstrv, getrecord, break_words
use error_module, only: error, moreinfo, warning
implicit none

integer i, io, irec, isite, jneigh, icoord, istype, nwords, nwords2
integer iw(maxwords2),iw2(maxwords2)

logical readlat, readcellvecs, readrepcell, readnstypes, readncellsites
logical readcellsitetypes, readcellsitecoords, readneighborstr, readstypenames

character(lengthrecinp) recinput, tmpstr

! Parses an explicit lattice specification in lattice_input.dat

! Level-1 keywords allowed:
!    cell_vectors
!    repeat_cell
!    n_site_types
!    site_type_names
!    n_cell_sites
!    site_types
!    site_coordinates
!    neighboring_structure
!    end_neighboring_structure
!    end_lattice

readcellvecs = .false.
readrepcell = .false.
readnstypes = .false.
readstypenames = .false.
readncellsites = .false.
readcellsitetypes = .false.
readcellsitecoords = .false.
readneighborstr = .false.

io = 0

call getrecord(ilatread,recinput,irec,io)

do while (io >= 0)

    call break_words(recinput,' ',remchar,nwords,iw)

    if (nwords == 0) then
    
        call getrecord(ilatread,recinput,irec,io)
        cycle
    
    elseif (striccompare(recinput(iw(1):iw(2)),'lattice')) then

        moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
        'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
        call error(100)

    elseif (striccompare(recinput(iw(1):iw(2)),'end_lattice')) then
        
        exit

    elseif (striccompare(recinput(iw(1):iw(2)),'cell_vectors')) then            

        if (readcellvecs) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif
            
        if (nwords > 1)    then
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(2030)
        endif
            
        i = 0
        do while (i < 2 .and. io >= 0)
            call getrecord(ilatread,recinput,irec,io)
            call break_words(recinput,' ',remchar,nwords,iw)

            if (nwords == 0) then
                cycle
            elseif (str_check_expression(recinput,nwords,iw,'2R8',.true.)) then
                i = i + 1
                vcell(1,i) = str2dbl(recinput(iw(1):iw(2)))
                vcell(2,i) = str2dbl(recinput(iw(3):iw(4)))
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(2031)
            endif    
        enddo
        
        if (i < 2) then 
            moreinfo = 'Not all unit cell vectors were specified in line ' // trim(int2str(irec)) // '.'
            call error(2106)
        endif
        
        readcellvecs = .true.

    elseif (striccompare(recinput(iw(1):iw(2)),'repeat_cell')) then

        if (readrepcell) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif

        if (str_check_expression(recinput,nwords,iw,'AA/2I4',.true.)) then
            ncellsrep = (/ str2int(recinput(iw(3):iw(4))), &
                           str2int(recinput(iw(5):iw(6))) /)                    
            readrepcell = .true.            
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(2032)
        endif

    elseif (striccompare(recinput(iw(1):iw(2)),'n_site_types')) then

        if (readnstypes) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif

        if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
            nsitetypes = str2int(recinput(iw(3):iw(4)))
            readnstypes = .true.            
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(2033)
        endif

    elseif (striccompare(recinput(iw(1):iw(2)),'site_type_names')) then

        if (readstypenames) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif

        if (.not.readnstypes) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(2034)
        endif
            
        if (nwords >= nsitetypes + 1) then
        
            if (nwords > nsitetypes + 1) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call warning(2001)
            endif
        
            allocate(sitetypenames(nsitetypes))
            do i = 1,nsitetypes
                sitetypenames(i) = recinput(iw(2*i+1):iw(2*i+2))
            enddo
            readstypenames = .true.

        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(2035)

        endif

    elseif (striccompare(recinput(iw(1):iw(2)),'n_cell_sites')) then

        if (readncellsites) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif

        if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
            ncellsites = str2int(recinput(iw(3):iw(4)))
            readncellsites = .true.
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(2036)
        endif

    elseif (striccompare(recinput(iw(1):iw(2)),'site_types')) then

        if (readcellsitetypes) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif
        
        if (.not.readncellsites) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(2037)
        elseif (.not.readnstypes) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(2038)
        endif

        allocate(cellsitetype(ncellsites))
        
        if (str_check_expression(recinput,nwords,iw,'AA/' // trim(int2str(ncellsites)) // 'I4',.true.)) then
        
            do i = 1,ncellsites
                istype = str2int(recinput(iw(2*i+1):iw(2*i+2)))
                if (istype < 1 .or. istype > nsitetypes) then
                    moreinfo = 'Site type ' // trim(int2str(istype)) // ' is out of bounds 1:' // &
                    trim(int2str(nsitetypes)) // ' in line ' // trim(int2str(irec)) // '.'
                    call error(2040)
                endif
                cellsitetype(i) = istype
            enddo
        
        elseif (str_check_expression(recinput,nwords,iw,trim(int2str(ncellsites+1)) // 'AA',.true.) &
            .and. readstypenames) then
        
            do i = 1,ncellsites
                istype = findinstrv(sitetypenames,1,1,nsitetypes,recinput(iw(2*i+1):iw(2*i+2)))
                if (istype < 0) then
                    moreinfo = 'Invalid site type name "' // recinput(iw(2*i+1):iw(2*i+2)) // &
                               '" in line ' // trim(int2str(irec)) // '.'
                    call error(2040)
                endif
                cellsitetype(i) = istype
            enddo            
            
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(2039)
            
        endif

        readcellsitetypes = .true.

    elseif (striccompare(recinput(iw(1):iw(2)),'site_coordinates')) then

        if (readcellsitecoords) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif
        
        if (.not.readncellsites) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(2041)
        else
            allocate(xfcellsite(ncellsites))
            allocate(yfcellsite(ncellsites))
        endif

        i = 0
        do
            call getrecord(ilatread,recinput,irec,io)
            call break_words(recinput,' ',remchar,nwords,iw)

            if (nwords == 0) then
                cycle
            elseif (str_check_expression(recinput,nwords,iw,'2R8',.true.)) then
                i = i + 1
                xfcellsite(i) = str2dbl(recinput(iw(1):iw(2)))
                yfcellsite(i) = str2dbl(recinput(iw(3):iw(4)))
                if (i == ncellsites) exit
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(2042)
            endif    
        enddo

        readcellsitecoords = .true.

    elseif (striccompare(recinput(iw(1):iw(2)),'neighboring_structure')) then

        if (readneighborstr) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif
        
        if (.not.readncellsites) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(2043)
        else
            allocate(cellsiteneighsf(ncellsites,0:ncellsites))
            allocate(cellsiteneighnn(ncellsites,0:ncellsites))
            allocate(cellsiteneighne(ncellsites,0:ncellsites))
            allocate(cellsiteneighee(ncellsites,0:ncellsites))
            allocate(cellsiteneighse(ncellsites,0:ncellsites))
            allocate(cellsiteneighss(ncellsites,0:ncellsites))
            allocate(cellsiteneighsw(ncellsites,0:ncellsites))
            allocate(cellsiteneighww(ncellsites,0:ncellsites))
            allocate(cellsiteneighnw(ncellsites,0:ncellsites))
            do isite = 1,ncellsites
                cellsiteneighsf(isite,0) = 0
                cellsiteneighnn(isite,0) = 0
                cellsiteneighne(isite,0) = 0
                cellsiteneighee(isite,0) = 0
                cellsiteneighse(isite,0) = 0
                cellsiteneighss(isite,0) = 0
                cellsiteneighsw(isite,0) = 0
                cellsiteneighww(isite,0) = 0
                cellsiteneighnw(isite,0) = 0
            enddo
        endif
        
        do
            call getrecord(ilatread,recinput,irec,io)

            call break_words(recinput,' ',remchar,nwords,iw)

            if (nwords == 0) then            
                cycle
                
            elseif (nwords == 2) then

                tmpstr = recinput(iw(1):iw(2))
                call break_words(trim(tmpstr),'-',remchar,nwords2,iw2)

                if (str_check_expression(tmpstr,nwords2,iw2,'2I4',.true.)) then
                    isite = str2int(tmpstr(iw2(1):iw2(2)))
                    jneigh = str2int(tmpstr(iw2(3):iw2(4)))
                else
                    moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(2044)
                endif                

                if (striccompare(recinput(iw(3):iw(4)),'self')) then

                    cellsiteneighsf(isite,0) = cellsiteneighsf(isite,0) + 1
                    cellsiteneighsf(isite,cellsiteneighsf(isite,0)) = jneigh
                    
                    cellsiteneighsf(jneigh,0) = cellsiteneighsf(jneigh,0) + 1
                    cellsiteneighsf(jneigh,cellsiteneighsf(jneigh,0)) = isite

                elseif (striccompare(recinput(iw(3):iw(4)),'north')) then

                    cellsiteneighnn(isite,0) = cellsiteneighnn(isite,0) + 1
                    cellsiteneighnn(isite,cellsiteneighnn(isite,0)) = jneigh
                    
                    cellsiteneighss(jneigh,0) = cellsiteneighss(jneigh,0) + 1
                    cellsiteneighss(jneigh,cellsiteneighss(jneigh,0)) = isite
                    
                elseif (striccompare(recinput(iw(3):iw(4)),'northeast')) then

                    cellsiteneighne(isite,0) = cellsiteneighne(isite,0) + 1
                    cellsiteneighne(isite,cellsiteneighne(isite,0)) = jneigh
                    
                    cellsiteneighsw(jneigh,0) = cellsiteneighsw(jneigh,0) + 1
                    cellsiteneighsw(jneigh,cellsiteneighsw(jneigh,0)) = isite

                elseif (striccompare(recinput(iw(3):iw(4)),'east')) then

                    cellsiteneighee(isite,0) = cellsiteneighee(isite,0) + 1
                    cellsiteneighee(isite,cellsiteneighee(isite,0)) = jneigh
                    
                    cellsiteneighww(jneigh,0) = cellsiteneighww(jneigh,0) + 1
                    cellsiteneighww(jneigh,cellsiteneighww(jneigh,0)) = isite

                elseif (striccompare(recinput(iw(3):iw(4)),'southeast')) then

                    cellsiteneighse(isite,0) = cellsiteneighse(isite,0) + 1
                    cellsiteneighse(isite,cellsiteneighse(isite,0)) = jneigh
                    
                    cellsiteneighnw(jneigh,0) = cellsiteneighnw(jneigh,0) + 1
                    cellsiteneighnw(jneigh,cellsiteneighnw(jneigh,0)) = isite

                else

                    moreinfo = 'Invalid directive "' // recinput(iw(3):iw(4)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(2044)                

                endif

            elseif (nwords == 1 .and. striccompare(recinput(iw(1):iw(2)),'end_neighboring_structure')) then            
                
                exit
    
            else
                
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(2044)                
            
            endif
            
        enddo
        
        readneighborstr = .true.

    else    
        
        moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
        call error(2100)
    
    endif

    call getrecord(ilatread,recinput,irec,io)

enddo

readlat = readcellvecs .and. readrepcell .and. readnstypes .and. readncellsites .and. &
          readcellsitetypes .and.readcellsitecoords .and. readneighborstr

! Check for completeness of input
moreinfo = ''
if (.not.readlat) then

    if (.not. readcellvecs) moreinfo = trim(moreinfo) // 'Real array cell_vectors for periodic cell lattice not encountered.' // char(13) // char(10)
    if (.not. readrepcell) moreinfo = trim(moreinfo) // 'Integer vector repeat_cell for periodic cell lattice not encountered.' // char(13) // char(10)
    if (.not. readnstypes) moreinfo = trim(moreinfo) // 'Integer n_site_types for periodic cell lattice not encountered.' // char(13) // char(10)
    if (.not. readncellsites) moreinfo = trim(moreinfo) // 'Integer n_cell_sites for periodic cell lattice not encountered.' // char(13) // char(10)
    if (.not. readcellsitetypes) moreinfo = trim(moreinfo) // 'Integer vector site_types for periodic cell lattice not encountered.' // char(13) // char(10)
    if (.not. readcellsitecoords) moreinfo = trim(moreinfo) // 'Real array site_coordinates for periodic cell lattice not encountered.' // char(13) // char(10)
    if (.not. readneighborstr) moreinfo = trim(moreinfo) // 'Section neighboring_structure for periodic cell lattice not encountered.' // char(13) // char(10)

    call error(2045)
    
endif

! Find maximum coordination number

maxcoord = 0

do i = 1,ncellsites

    icoord = cellsiteneighsf(i,0) + &
        cellsiteneighnn(i,0) + cellsiteneighne(i,0) + cellsiteneighee(i,0) + cellsiteneighse(i,0) + &
        cellsiteneighss(i,0) + cellsiteneighsw(i,0) + cellsiteneighww(i,0) + cellsiteneighnw(i,0)

    if (icoord > maxcoord) then
        maxcoord = icoord
    endif

enddo

! Check self-consistency of unit cell
call check_unit_cell()

return

end subroutine read_periodic_lattice_cell

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine read_explicit_lattice(irec,io,periodiclat)

use constants_module, only: maxwords2, lengthrecinp, remchar, ilatread,        &
                            clatfname
use error_module, only: error, moreinfo, warning
use parser_module, only: int2str, striccompare, str2int, str2dbl, findinstrv,  &
                         str_is_integer, str_check_expression, getrecord,      &
                         break_words
implicit none

integer i, io, irec, istype, nsitesread, j, isite, jneigh
integer nmxcoord, nwords
integer iw(maxwords2)

logical readnsites, readnstypes, readstypenames, readmaxcoord, readlatstruct,  &
        readcellvecs, periodiclat

character(lengthrecinp) recinput

! Parses an explicit lattice specification in lattice_input.dat

! Level-1 keywords allowed:
!    cell_vectors
!    n_sites
!    max_coord
!    n_site_types
!    site_type_names
!    lattice_structure
!    end_lattice_structure
!    end_lattice

periodiclat = .false. ! if cell vector information is parsed,
                      ! this will be set to true
                      
readnsites = .false.
readmaxcoord = .false.
readnstypes = .false.
readstypenames = .false.
readlatstruct = .false.
readcellvecs = .false.

io = 0
maxcoord = 0

call getrecord(ilatread,recinput,irec,io)

do while (io >= 0)

    call break_words(recinput,' ',remchar,nwords,iw)

    if (nwords == 0) then
    
        call getrecord(ilatread,recinput,irec,io)
        cycle

    elseif (striccompare(recinput(iw(1):iw(2)),'lattice')) then

        moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
        'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
        call error(100)

    elseif (striccompare(recinput(iw(1):iw(2)),'end_lattice')) then
        
        exit
    
    elseif (striccompare(recinput(iw(1):iw(2)),'cell_vectors')) then            

        if (readcellvecs) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif
            
        if (nwords > 1)    then
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(2030)
        endif
            
        i = 0
        do while (i < 2 .and. io >= 0)
            call getrecord(ilatread,recinput,irec,io)
            call break_words(recinput,' ',remchar,nwords,iw)

            if (nwords == 0) then
                cycle
            elseif (str_check_expression(recinput,nwords,iw,'2R8',.true.)) then
                i = i + 1
                vcell(1,i) = str2dbl(recinput(iw(1):iw(2)))
                vcell(2,i) = str2dbl(recinput(iw(3):iw(4)))
            else
                moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(2031)
            endif    
        enddo
        
        if (i < 2) then 
            moreinfo = 'Not all unit cell vectors were specified in line ' // trim(int2str(irec)) // '.'
            call error(2106)
        endif
        
        readcellvecs = .true.
        periodiclat = .true.
        ncellsrep = (/1,1/)

    elseif (striccompare(recinput(iw(1):iw(2)),'n_sites')) then

        if (readnsites) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif

        if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
            nsites = str2int(recinput(iw(3):iw(4)))
            readnsites = .true.
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(2060)
        endif
            

    elseif (striccompare(recinput(iw(1):iw(2)),'max_coord')) then

        if (readmaxcoord) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif

        if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
            nmxcoord = str2int(recinput(iw(3):iw(4)))
            readmaxcoord = .true.
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(2061)
        endif

    elseif (striccompare(recinput(iw(1):iw(2)),'n_site_types')) then

        if (readnstypes) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif
        
        if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
            nsitetypes = str2int(recinput(iw(3):iw(4)))
            readnstypes = .true.        
        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(2033)
        endif
    
    elseif (striccompare(recinput(iw(1):iw(2)),'site_type_names')) then

        if (readstypenames) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif

        if (.not.readnstypes) then
            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
            call error(2034)
        endif
        
        if (nwords >= nsitetypes + 1) then
        
            if (nwords > nsitetypes + 1) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call warning(2001)
            endif
        
            allocate(sitetypenames(nsitetypes))
            do i = 1,nsitetypes
                sitetypenames(i) = recinput(iw(2*i+1):iw(2*i+2))
            enddo
            readstypenames = .true.

        else
            moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
            call error(2035)

        endif

    elseif (striccompare(recinput(iw(1):iw(2)),'lattice_structure')) then
        
        if (readlatstruct) then
            moreinfo = 'Keyword ' // recinput(iw(1):iw(2)) // ' has already been parsed; yet, it ' // &
            'appears again in line ' // trim(int2str(irec)) // ' in ' // trim(clatfname) // '.'
            call error(100)
        endif

        if (readnsites .and. readnstypes .and. readmaxcoord) then    
            allocate(xsite(nsites))
            allocate(ysite(nsites))
            allocate(sitetype(nsites))
            allocate(siteneighb1(nsites,0:nmxcoord))
        else
            call error(2062)    
        endif

        ! Initialize sitetype vector and nsitesread with zeros
        do i = 1,nsites
            sitetype(i) = 0
        enddo

        nsitesread = 0
                
        ! Start reading the information for each site 
        do while (io >= 0)

            call getrecord(ilatread,recinput,irec,io)

            call break_words(recinput,' ',remchar,nwords,iw)
            
            if (nwords == 0) then
                cycle
            
            elseif (striccompare(recinput(iw(1):iw(2)),'end_lattice_structure')) then
                exit

            elseif (str_check_expression(recinput,nwords,iw,'I4/2R8/AA/I4',.false.)) then
                
                nsitesread = nsitesread + 1
                
                isite = str2int(recinput(iw(1):iw(2)))
                
                if (isite > nsites .or. isite < 1) then
                    moreinfo = 'Site ' // trim(int2str(isite)) // ' out of bounds 1:' // trim(int2str(nsites)) // &
                               ' in line ' // trim(int2str(irec)) // '.'
                    call error(2063)       
                endif
                    
                xsite(isite) = str2dbl(recinput(iw(3):iw(4)))
                ysite(isite) = str2dbl(recinput(iw(5):iw(6)))                    
                    
                if (str_is_integer(recinput(iw(7):iw(8)),4)) then ! Site type specified by an integer
                
                    istype = str2int(recinput(iw(7):iw(8)))
                    if (istype < 1 .or. istype > nsitetypes) then
                        moreinfo = 'Site type ' // trim(int2str(istype)) // ' is out of bounds 1:' // &
                        trim(int2str(nsitetypes)) // ' in line ' // trim(int2str(irec)) // '.'
                        call error(2040)
                    endif
                
                else  ! Site type specified by the site type name
                    if (.not. allocated(sitetypenames)) then
                        call error(2064)
                    endif
                
                    istype = findinstrv(sitetypenames,1,1,nsitetypes,recinput(iw(7):iw(8)))
                    if (istype < 0) then
                        moreinfo = 'Invalid site type name "' // recinput(iw(7):iw(8)) // &
                                   '" in line ' // trim(int2str(irec)) // '.'
                        call error(2040)
                    endif
                endif  
                    
                if (sitetype(isite) == 0) then ! easy way to check for duplicate site definitions
                    sitetype(isite) = istype
                else
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(2065)
                endif
                    
                siteneighb1(isite,0) = str2int(recinput(iw(9):iw(10)))
                    
                if (siteneighb1(isite,0) > nmxcoord .or. siteneighb1(isite,0) < 0) then ! check if: 0 <= number of neihbors <= max_coord
                    moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                    call error(2063)
                endif
                
                if (nwords < siteneighb1(isite,0)+5) then ! check if all neighbors are listed
                    moreinfo = 'Expected ' // trim(int2str(siteneighb1(isite,0)+5)) // ' directives in expression "' // &
                               recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                    call error(2066)
                endif

                if (maxcoord < siteneighb1(isite,0)) maxcoord = siteneighb1(isite,0)
                                    
                do j = 1,siteneighb1(isite,0)
                    if (str_is_integer(recinput(iw(2*j+9):iw(2*j+10)),4)) then                    
                        jneigh = str2int(recinput(iw(2*j+9):iw(2*j+10)))
                        if (jneigh > nsites .or. isite < 1) then
                            moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                            call error(2063)       
                        endif
                        siteneighb1(isite,j) = jneigh
                    else
                        moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                        call error(2066)
                    endif
                enddo
                
            else
                
                moreinfo = 'Invalid expression "'// recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
                call error(2066)    
                
            endif

        enddo
        
        readlatstruct = .true.

    else    
        
        moreinfo = 'Unknown or invalid keyword ' // recinput(iw(1):iw(2)) // ' in line ' // trim(int2str(irec)) // '.'
        call error(2100)
    
    endif
    
    call getrecord(ilatread,recinput,irec,io)

enddo

! Check for completeness of input
if (nsitesread < nsites) call error(2067)

! Check self-consistency of lattice structure
call check_lattice()

return

end subroutine read_explicit_lattice

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_unit_cell()

use parser_module, only: int2str
use error_module, only: error, moreinfo, warning
implicit none

integer i, j, k, isite, jneigh

logical check

! Check that each neighbor is mentioned only once in each neighbor list

do i = 1,ncellsites

    do j = 1,cellsiteneighsf(i,0)
        do k = j+1,cellsiteneighsf(i,0)

            if (cellsiteneighsf(i,j) == cellsiteneighsf(i,k)) then
                moreinfo = 'Site ' // &
                    trim(int2str(cellsiteneighsf(i,j))) // &
                    ' appears more than once in self-cell neighbor list of site ' // &
                    trim(int2str(i))
                call error(2101)
                
            endif
        enddo
    enddo

    do j = 1,cellsiteneighnn(i,0)
        do k = j+1,cellsiteneighnn(i,0)

            if (cellsiteneighnn(i,j) == cellsiteneighnn(i,k)) then
                moreinfo = 'Site ' // &
                    trim(int2str(cellsiteneighnn(i,j))) // &
                    ' appears more than once in north neighbor list of site ' // &
                    trim(int2str(i))
                call error(2101)
                
            endif
        enddo
    enddo

    do j = 1,cellsiteneighne(i,0)
        do k = j+1,cellsiteneighne(i,0)

            if (cellsiteneighne(i,j) == cellsiteneighne(i,k)) then
                moreinfo = 'Site ' // &
                    trim(int2str(cellsiteneighne(i,j))) // &
                    ' appears more than once in north-east neighbor list of site ' // &
                    trim(int2str(i))
                call error(2101)
                
            endif
        enddo
    enddo

    do j = 1,cellsiteneighee(i,0)
        do k = j+1,cellsiteneighee(i,0)

            if (cellsiteneighee(i,j) == cellsiteneighee(i,k)) then
                moreinfo = 'Site ' // &
                    trim(int2str(cellsiteneighee(i,j))) // &
                    ' appears more than once in east neighbor list of site ' // &
                    trim(int2str(i))
                call error(2101)
                
            endif
        enddo
    enddo

    do j = 1,cellsiteneighse(i,0)
        do k = j+1,cellsiteneighse(i,0)

            if (cellsiteneighse(i,j) == cellsiteneighse(i,k)) then
                moreinfo = 'Site ' // &
                    trim(int2str(cellsiteneighse(i,j))) // &
                    ' appears more than once in south-east neighbor list of site ' // &
                    trim(int2str(i))
                call error(2101)
                
            endif
        enddo
    enddo

enddo

! Check for consistency of self-neighboring structure

do isite = 1,ncellsites ! for every site...
    do j = 1,cellsiteneighsf(isite,0) ! loop through neighbors
        
        jneigh = cellsiteneighsf(isite,j) ! for neighbor jneigh...
        check = .false.

        do k = 1,cellsiteneighsf(jneigh,0) ! try to find isite in jneigh's neighbor list
            if (isite == cellsiteneighsf(jneigh,k)) then
                check = .true. ! isite was found
                exit
            endif
        enddo

        if (.not.(check)) then
            moreinfo = 'Neighboring-site ' // trim(int2str(jneigh)) // &
                ' of site ' // trim(int2str(isite)) // ' is not a' // &
                ' neighbor of the latter.'
            call error(2102)
        endif

    enddo

enddo

return

end subroutine check_unit_cell

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_lattice()

use parser_module, only: int2str
use error_module, only: error, moreinfo, warning
implicit none

integer isite, j, jneigh, k, kneigh

logical foundisite

! Checking the lattice structure

! Check that the neighboring list has unique entries different than the original site
do isite = 1,nsites
    do j = 1,siteneighb1(isite,0)
        jneigh = siteneighb1(isite,j)
        
        if (jneigh == isite) then
            moreinfo = 'Site ' // trim(int2str(isite)) // &
                        ' neighbors with itself.'
            call warning(2103)    
        endif
                
        do k = j+1,siteneighb1(isite,0)
            kneigh = siteneighb1(isite,k)
            if (jneigh == kneigh) then
                moreinfo = 'Encountered neighboring site ' // trim(int2str(kneigh)) // &
                            ' twice for site ' // trim(int2str(isite)) // '.'
                call warning(2104)    
            endif
        enddo
    enddo
enddo        

! Check the parity between neighboring sites
do isite = 1,nsites
    do j = 1,siteneighb1(isite,0)
    
        jneigh = siteneighb1(isite,j)
        foundisite = .false. ! assume that isite does not appear in jneigh's list of neighbors
      
        do k = 1,siteneighb1(jneigh,0)
            kneigh = siteneighb1(jneigh,k)
            if (isite == kneigh) then
                foundisite = .true.
                exit
            endif
        enddo
        
        if (.not.foundisite) then
            moreinfo = 'Site ' // trim(int2str(isite)) // ' references ' // &
               'site ' // trim(int2str(jneigh)) // ' as a neighbor, but not the converse.'
            call error(2105)
        endif
        
    enddo
enddo

end subroutine check_lattice

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module lattice_setup_module
