! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module graph_functions_module

implicit none

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine find_dlevel_neighbors(neighbors1,d,nstart,startsites, &
                nneigh,indxs,neighlist)

implicit none

integer, intent(in) :: startsites(:)
integer, intent(in) :: neighbors1(:,0:)
integer, intent(in) :: d, nstart
integer, intent(out):: nneigh
integer, intent(out):: neighlist(:), indxs(0:d)

integer i, j, k, m, nneighpre, nneighcur, nsiteneighs
integer jneigh
logical newneigh

! Finds all up-to-d-level nearest neighbors of nstart sites contained in vector 
! startsites(:)

! Start from nstart sites stored in vector startsites
nneigh = nstart
do i = 1,nstart
    neighlist(i) = startsites(i)
enddo

nneighpre = 1
nneighcur = nneigh

indxs(0) = nneighcur

do i = 1,d
    
    ! For the sites of each level...
    do j = nneighpre,nneighcur
    
        nsiteneighs = neighbors1(neighlist(j),0)
                
    
        ! ... loop over all their neighbors
        do k = 1,nsiteneighs
            
            jneigh = neighbors1(neighlist(j),k)
            newneigh = .true. ! assume this neighbor is a new site,
                              ! not already included in the list
            
            do m = 1,nneigh
                ! if this site is already found in the list, flag and exit
                if (neighlist(m) == jneigh) then
                    newneigh = .false.
                    exit
                endif
            enddo
            
            if (newneigh) then
                nneigh = nneigh + 1
                neighlist(nneigh) = jneigh
            endif
            
        enddo         
        
    enddo
    
    nneighpre = nneighcur+1
    nneighcur = nneigh

    indxs(i) = nneighcur

enddo

end subroutine find_dlevel_neighbors

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine find_max_level(nsit,neighbors1,dmax,nstart,startsites)

implicit none

integer, intent(in)  :: nsit, nstart
integer, intent(in)  :: startsites(:)
integer, intent(in)  :: neighbors1(:,0:)
integer, intent(out) :: dmax

integer i, j, jneigh, k, m
integer nneigh, nneighpre, nneighcur, nsiteneighs
integer neighlist(nsit)

logical newneigh

! Finds the maximum level of the graph pattern contained in neighlist
! The maximum level is the maximum exponent dmax in the expression
!     v0 + A*v0 + A^2*v0 + ... + A^dmax*v0
! so that we get a v0 with non zero elements everywhere.
! A is the adjacency matrix of the graph and we assume that the graph is connected.
! Note that the level depends on the starting points contained in vector startsites.
! The number of starting points is stored in nstart.

nneigh = nstart
neighlist(:nstart) = startsites(:nstart)

nneighpre = 1
nneighcur = nneigh

dmax = 0
do while (nneigh < nsit)
    
    dmax = dmax + 1
    
    ! For the sites of each level...
    do j = nneighpre,nneighcur
    
        nsiteneighs = neighbors1(neighlist(j),0)
            
        ! ... loop over all their neighbors
        do k = 1,nsiteneighs
            
            jneigh = neighbors1(neighlist(j),k)
            newneigh = .true. ! assume this neighbor is a new site,
                              ! not already included in the list
            
            do m = 1,nneigh
                ! if this site is already found in the list, flag and exit
                if (neighlist(m) == jneigh) then
                    newneigh = .false.
                    exit
                endif
            enddo
            
            if (newneigh) then
                nneigh = nneigh + 1
                neighlist(nneigh) = jneigh
            endif
            
        enddo         
        
    enddo
    
    nneighpre = nneighcur+1
    nneighcur = nneigh

enddo

end subroutine find_max_level

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine create_next_permutation(total,chosen,isfixed,isspecific,isused,permutat, &
                                   nindx,lastindx,nlevel,speclevs,nonspeclevs, &
                                   checkflag,statusflag,completeflag)

logical, intent(in)    :: isspecific(:)
logical, intent(in)    :: checkflag
integer, intent(in)    :: total, chosen
integer, intent(out)   :: nlevel
logical, intent(inout) :: isfixed(:), isused(total)
logical, intent(inout) :: statusflag, completeflag
integer, intent(inout) :: permutat(:)
integer, intent(inout) :: nindx, lastindx
integer, intent(inout) :: speclevs(0:), nonspeclevs(0:)

integer i, j, prevlevel, nextobject

! Subroutine that generates permutations of "chosen" out of "total" objects.
! The permutations are returned up to level ncompl and are being checked by an 
! external subroutine to see if they satisfy some criterion. If the partially
! generated permutation satisfies the criterion, then the next level is filled.
! If not, the next object is considered for the current level.

! Input checkflag can have the values .true. or .false. depending on whether 
! the current partial permutation satisfies the criterion or not, respectively.
! Depending on its value, the next partial or complete permutation is generated.

! I/O variable statusflag takes the value .true. to indicate that more permutations
! can be generated or false to show that the last permutation has just been generated.

! Output completeflag is set to .true. if all levels are filled (complete permutation).

! In the first call of the subroutine statusflag must be set to .true. and .completeflag to .false.

! Example code: total = 10 and say that we want to choose 5 objects out of which the 1st and 3rd are fixed to 1 and 7.
! The following code illustrates how to call the subroutine to generate all permutations that satisfy a criterion.

!    integer total,chosen,permutat(5),speclevs(5),nonspeclevs(5),nlevel,nindx,lastindx,npermutat
!    logical isfixed(5),isspecific(5),isused(10),checkflag,statusflag,completeflag
!
!    total = 10
!    chosen = 5
!
!    permutat = (/1,0,7,0,0/)
!    isspecific = .true.
!    isspecific(2) = .false.
!    isspecific(4) = .false.
!
!    completeflag = .false.
!    statusflag = .true.
!    nindx = 0
!    npermutat = 0
!
!    do while (statusflag)
!
!        call create_next_permutation(total,chosen,isfixed,isspecific,isused,permutat, &
!                                       nindx,lastindx,nlevel,speclevs,nonspeclevs, &
!                                       checkflag,statusflag,completeflag)
!        if (.not.statusflag) exit
!        write(*,*) permutat
!        !checkflag = .true. ! always checks
!        read(*,*) checkflag ! the user checks the partial permutation
!        if (completeflag .and. checkflag) then
!            npermutat = npermutat + 1 ! count complete permutations
!            write(*,*) 'valid permutation found!'
!        endif
!        continue
!        
!    enddo
!    write(*,*) 'Found ',npermutat,' valid permutations!'
!
!    stop

! Array dimensions:
! fixed(1:chosen)
! permutat(1:chosen)
! isused(1:total)

if (.not.statusflag) return ! there are no more permutations to be generated
if (completeflag .and. nindx == 0) then
    statusflag = .false.
    return ! all levels are fixed and the permutation has been returned once
endif

completeflag = .false.

! *** Initialization if nindx == 0

if (nindx == 0) then

    isused(:total) = .false. 
    
    speclevs(0) = 0
    nonspeclevs(0) = 0
    do i = 1,chosen
    
        ! Note if level i is fixed and mark the object it contains as used
        if (permutat(i) == 0) then
            isfixed(i) = .false.
        else
            isfixed(i) = .true.
            isused(permutat(i)) = .true.
        endif        

        ! Populate the index vectors that give the specific and non-specific levels
        if (.not.isfixed(i)) then
            if (isspecific(i)) then
                speclevs(0) = speclevs(0) + 1
                speclevs(speclevs(0)) = i
            else
                nonspeclevs(0) = nonspeclevs(0) + 1
                nonspeclevs(nonspeclevs(0)) = i
            endif
        else
            if (.not.isspecific(i)) then
                ! This is invalid and means that there is a bug in the code, that is why it was 
                ! not put in the error handling subroutine as a legitimate error
                write(2222,*) 'FATAL ERROR from create_next_permutation subroutine of graph_functions_module:'
                write(2222,*) 'Cannot have a non-specific fixed level.'
                stop
            endif
        endif
        
    enddo

    lastindx = speclevs(0) + nonspeclevs(0)
    
    if (lastindx == 0) then
        completeflag = .true.
        return ! all levels are fixed
    endif
    
    if (speclevs(0) > 0) then
        ! populate the first non-fixed and specific level with the next available object
        nindx = 1
        nlevel = speclevs(nindx)
        do j = 1,total
            if (.not.isused(j)) then
                permutat(nlevel) = j
                isused(j) = .true.
                exit
            endif
        enddo
    
    else
        ! populate the first non-fixed and non-specific level with the next available object
        nindx = 1
        nlevel = nonspeclevs(nindx)
        do j = 1,total
            if (.not.isused(j)) then
                permutat(nlevel) = j
                isused(j) = .true.
                exit
            endif
        enddo    
    
    endif
        
    if (nindx == lastindx) completeflag = .true.
    return    
    
endif

! *** Propagation: make next permutation

if (checkflag) then ! The current permutation passed the test
    
    if (nindx < lastindx) then ! the current permutation is unfinished so we will continue with the next level
    
        nindx = nindx + 1
        
        ! Specify the next level to be populated
        if (nindx <= speclevs(0)) then
            nlevel = speclevs(nindx)
        else
            nlevel = nonspeclevs(nindx-speclevs(0))
        endif
        
        ! Populate the next level with the next available object
        do i = 1,total
            if (.not.isused(i)) then
                permutat(nlevel) = i
                isused(i) = .true.
                if (nindx == lastindx) completeflag = .true.
                return
            endif
        enddo

    else ! the current permutation is a finished one
    
        if (nonspeclevs(0) > 0) then ! if there are non-specific levels we are done: we don't want to iterate them more
            do i = 1,nonspeclevs(0)
                nlevel = nonspeclevs(i)
                isused(permutat(nlevel)) = .false.
                permutat(nlevel) = 0                
            enddo
            nindx = speclevs(0)
        endif

        nlevel = speclevs(nindx)
        
        ! If there are nonspecific levels, we have flushed them out and consider the last specific level.
        ! If there are only specific levels, we have reached the last one. 
        ! We will continue execution and try to populate the current level
        
        ! Check if there exists a next object to put on that level
        nextobject = 0
        do i = permutat(nlevel)+1,total
            if (.not.isused(i)) then
                nextobject = i
                exit
            endif
        enddo

        if (nextobject > 0) then
            
            ! Put the next object in the current level
            isused(permutat(nlevel)) = .false.
            permutat(nlevel) = nextobject
            isused(nextobject) = .true.
            if (nindx == lastindx) completeflag = .true.
            return
            
        endif
        
    endif
    
endif
       
! If we have reached this part of the code it means that 
! - either the (possibly unfinished) permutation did not pass the test, 
!   in which case we can consider the next object for that level
! - or that it did, but there are no more available objects for the current level,
!   in which case we have to rewind


if (.not. checkflag) then ! The current permutation did not pass the test

    ! Check if there exists a next object to put on that level
    nextobject = 0
    do i = permutat(nlevel)+1,total
        if (.not.isused(i)) then
            nextobject = i
            exit
        endif
    enddo

    if (nextobject > 0) then
        
        ! Put the next object in the current level
        isused(permutat(nlevel)) = .false.
        permutat(nlevel) = nextobject
        isused(nextobject) = .true.
        if (nindx == lastindx) completeflag = .true.
        return
        
    endif

endif
    
! At this point, there are no more available objects for the current level and 
! we have to rewind

prevlevel = 0
do
    ! Check if there exists a previous level
    if (nindx > 1) then

        nindx = nindx - 1

        if (nindx <= speclevs(0)) then
            prevlevel = speclevs(nindx)
        else
            prevlevel = nonspeclevs(nindx-speclevs(0))
        endif

    endif

    if (prevlevel > 0) then ! There exists a previous level
        
        isused(permutat(nlevel)) = .false. ! release the current level from the occupying object
        permutat(nlevel) = 0
        nlevel = prevlevel
        
        ! Check if there exists a next object to put in the previous level                
        nextobject = 0
        do i = permutat(nlevel)+1,total
            if (.not.isused(i)) then
                nextobject = i
                exit
            endif
        enddo
        
        if (nextobject > 0) then
            
            ! Put the next object in the previous level
            isused(permutat(nlevel)) = .false.
            permutat(nlevel) = i
            isused(i) = .true.
            exit
        
        else
            
            prevlevel = 0
            cycle
            
        endif

    else ! no previous level, last permutation was generated
    
        statusflag = .false.
        completeflag = .false. ! completeflag is set to false here since no new permutation was generated;
                               ! thus, the last permutation will not be double-counted by an external subroutine
        return            
    
    endif
    
enddo


end subroutine create_next_permutation

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!subroutine find_graph_isomorphs( &
!                nsit,neighbors1, &       ! total number of sites and neighboring structure
!                nobjsites,objsites, &    ! the number of candidate sites and the vector containing them
!                nspat,neighborspattrn, & ! number of sites in the pattern, neighboring structure
!                fixedsites, &
!                cntrpatrn,validpatrns)
!
!implicit none
!
!integer i, nlevel, lastlevel, cntrpatrn
!integer nsit, nspat, nobjsites
!integer sitemap(nspat), checksites(0:nspat)
!integer fixedsites(:), objsites(:), validpatrns(:,:)
!
!integer neighbors1(:,0:)
!integer neighborspattrn(:,0:)
!
!logical isfixed(nspat), isused(nobjsites), checkflag, statusflag, completeflag
!
!! Subroutine that finds isomorphisms.
!!
!! This works by finding the permutations of nobjsites that match the graph pattern
!! example:
!!                      1   2   3   4   5
!!       objsites = [ 874 322  43 102  12]
!!
!! the ith element of sitemap gives the objsites index of the site that corresponds/matches 
!! site i of the pattern, for instance if:
!!       sitemap  = [ 2 4 1 ]
!! Then: site 1 of the pattern corresponds to site 322 of the lattice,
!!       site 2 of the pattern corresponds to site 102 of the lattice,
!!       site 3 of the pattern corresponds to site 874 of the lattice
!!
!! Say now that we want to fix sites 1 and 3 of the pattern to objsites(4) and objsites(5), 
!! respectively. Then we would set
!!
!!       fixedsites = [ 4 0 5 ]
!
!! Example code:
!
!! integer i, j, nobjsites, nspat, cntrpatrn
!! integer, allocatable :: objsites(:)
!! integer, allocatable :: fixedsites(:)
!! integer, allocatable :: neighborspattrn(:,:)
!! integer, allocatable :: validpatrns(:,:)
!! 
!! ! Finds triangular patterns in a lattice. One has to generate a 
!! ! hexagonal lattice using the following syntax in lattice_input.dat:
!! !    lattice default_choice
!! !      hexagonal_periodic 1.25 8 8
!! !    end_lattice
!!
!! call read_lattice_setup()
!!
!! nobjsites = nsites
!! allocate(objsites(nobjsites))
!! do i = 1,nobjsites
!!     objsites(i) = i
!! enddo
!!
!! nspat = 3
!! allocate(neighborspattrn(nspat,0:nspat-1))
!! neighborspattrn(1,0:2) = (/2,2,3/)
!! neighborspattrn(2,0:2) = (/2,1,3/)
!! neighborspattrn(3,0:2) = (/2,1,2/)
!!
!! allocate(fixedsites(nspat))
!! fixedsites = (/56,0,71/)
!!
!! allocate(validpatrns(10000,nspat))
!!
!! call find_graph_isomorphs( &
!!             nsites,siteneighb1, & ! total number of sites and neighboring structure
!!             nobjsites,objsites, & ! the number of candidate sites and the vector containing them
!!             nspat,neighborspattrn, & ! number of sites in the pattern, neighboring structure
!!             fixedsites, &
!!             cntrpatrn,validpatrns)
!!
!! write(2222,*) cntrpatrn
!! do i = 1,cntrpatrn
!!     write(2222,*) (validpatrns(i,j), j = 1,nspat)
!! enddo
!! stop
!
!cntrpatrn = 0
!
!if (nobjsites < nspat) return
!
!! Populate sitemap and checksites vectors with the fixed sites
!checksites(0) = 0 ! The checksites array contains the sites that will be checked for 
!                  ! the isomorphism. This saves time, since every time we add a new
!                  ! site in the candidate pattern, we do not need to check all
!                  ! previously added sites.
!
!do i = 1,nspat
!    sitemap(i) = fixedsites(i);
!    if (fixedsites(i) /= 0) then
!        checksites(0) = checksites(0) + 1
!        checksites(checksites(0)) = i ! checksites numbering refers to the pattern
!    endif
!enddo
!
!! Check partial isomorphism for the fixed sites of this pattern (if more than one)
!if (checksites(0) > 1) then
!    call check_isomorph(nsit,neighbors1,nspat,neighborspattrn,sitemap,objsites,checksites,checkflag)
!else
!    checkflag = .true.
!endif
!
!if (.not.checkflag) return
!
!! Enter the loop in which permutations of sites are generated and the pattern is checked
!nlevel = 0
!statusflag = .true.
!completeflag = .false.
!call create_next_permutation(nobjsites,nspat,isfixed,isused,sitemap,nlevel, &
!                             lastlevel,checkflag,statusflag,completeflag)
!
!if (lastlevel == 0) then ! In this case all sites are fixed
!
!    if (checkflag) then ! If the pattern checks, add it to the valid patterns and return
!        cntrpatrn = cntrpatrn + 1
!        do i = 1,nspat
!            validpatrns(cntrpatrn,i) = objsites(sitemap(i))
!        enddo
!    endif
!    
!    return ! (return since there is nothing more to do)
!
!endif
!
!do while (statusflag)
!    
!    checksites(0) = 1
!    checksites(1) = nlevel
!
!    call check_isomorph(nsit,neighbors1,nspat,neighborspattrn,sitemap,objsites,checksites,checkflag)
!    
!    if (checkflag .and. completeflag) then
!        cntrpatrn = cntrpatrn + 1
!        do i = 1,nspat
!            validpatrns(cntrpatrn,i) = objsites(sitemap(i))
!        enddo
!    endif
!
!    call create_next_permutation(nobjsites,nspat,isfixed,isused,sitemap,nlevel, &
!                                 lastlevel,checkflag,statusflag,completeflag)
!
!enddo
!
!return
!
!end subroutine find_graph_isomorphs

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_isomorph(nsit,neighbors1, &
                nspat,neighborspattrn, &
                sitemap,objsites,checksites,checks)

implicit none

integer, intent(in)  :: nsit, nspat
integer, intent(in)  :: neighbors1(:,0:)
integer, intent(in)  :: neighborspattrn(:,0:)
integer, intent(in)  :: sitemap(nspat), objsites(:)
integer, intent(in)  :: checksites(0:)
logical, intent(out) :: checks

integer i, isite, j, jneigh, k

logical foundneigh


checks = .true.

do i = 1,checksites(0) ! For each of the sites to be checked

    isite = checksites(i)
    
    do j = 1,neighborspattrn(isite,0) ! find its neighbors in the pattern
    
        jneigh = neighborspattrn(isite,j)
        
        if (sitemap(jneigh) /= 0) then ! If that neighbor has been assigned 
                                       ! a site in the big lattice
            ! check that the same neighboring relation exists in the big lattice
            foundneigh = .false.
            do k = 1,neighbors1(objsites(sitemap(isite)),0)
                if (neighbors1(objsites(sitemap(isite)),k) == objsites(sitemap(jneigh))) then
                    foundneigh = .true.
                    exit
                endif
            enddo
            
            if (.not. foundneigh) then
                checks = .false.
                return
            endif
            
        endif
        
    enddo
    
enddo

return

end subroutine check_isomorph

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

real(8) function measureangle(x1,x2,x3,y1,y2,y3)

use constants_module, only: pi
implicit none

real(8), intent(in) :: x1, x2, x3, y1, y2, y3
real(8) u1, v1, u2, v2
real(8) innerprod, outerprod, costheta, sintheta, lengthprod, out1

! The range of this function is [-180.0,180.0]
! Be careful about the positive/negative convention.
! Angles are measured from vector 1-2 to 2-3:
!  3
!    \  angle measured as positive in the counter-clockwise direction
!     2 -- 1

! Use the following code to test if necessary:
!real(8) x1,x2,x3,y1,y2,y3, theta
!do theta = 0.d0,2*pi,2*pi/20.d0
!    x1 = 1.d0
!    y1 = 0.d0
!    x2 = 0.d0
!    y2 = 0.d0
!    x3 = dcos(theta)
!    y3 = dsin(theta)
!    write(*,*) measureangle(x1,x2,x3,y1,y2,y3)
!enddo
!pause
!stop

u1 = x1-x2;
v1 = y1-y2;

u2 = x3-x2;
v2 = y3-y2;

innerprod = u1*u2 + v1*v2;
outerprod = u1*v2 - v1*u2;
lengthprod = dsqrt(u1**2+v1**2)*dsqrt(u2**2+v2**2)

costheta = innerprod/lengthprod;
sintheta = outerprod/lengthprod;

sintheta = max(-1.d0,min(1.d0,sintheta)) ! numerical error may give values slightly out of [-1,1]
costheta = max(-1.d0,min(1.d0,costheta)) ! numerical error may give values slightly out of [-1,1]

if (dabs(sinTheta) < 1e-10 .and. dabs(cosTheta+1) < 1e-10) then
    measureangle = 180.d0
    return
elseif (dabs(sinTheta) < 1e-10 .and. dabs(cosTheta-1) < 1e-10) then
    measureangle = 0
    return
endif

out1 = dasin(sintheta)
if (sintheta < 0 .and. costheta < 0) then
    out1 = -pi-out1
elseif (sintheta > 0 .and. costheta < 0) then
    out1 = pi-out1
endif

measureangle = out1*180.d0/pi

return

end function measureangle

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module graph_functions_module
