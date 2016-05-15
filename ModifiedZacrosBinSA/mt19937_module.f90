! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

 module mt19937_module

! Mersenne Twister MT19937. Paper introducing the method:
! M. Matsumoto and T. Nishimura (1998) "Mersenne Twister: 
! "A 623-Dimensionally Equidistributed Uniform Pseudo-Random
! Number Generator", ACM Trans. Model. Comput. Simulat. 8(1): 3-30.

! Random number generator state vector
integer(4), allocatable :: rngstate(:)
! Random number generator period parameters
integer(4), parameter :: rngw = 32  ! word size
integer(4), parameter :: rngn = 624 ! degree of recurrence
integer(4), parameter :: rngm = 397 ! middle term
integer(4), parameter :: rngr = 31  ! separation point of one word
integer(4), parameter :: rnglo = 2147483647_4 ! B'01111111111111111111111111111111'
integer(4), parameter :: rngup = -rnglo-1 ! B'10000000000000000000000000000000'

! Random number generator tempering parameters
integer(4), parameter :: rnga = -1727483681_4 ! X'9908B0DF'
integer(4), parameter, dimension(0:1) :: rngaaugm = (/0,rnga/)
integer(4), parameter :: rngl = 18
integer(4), parameter :: rngu = 11
integer(4), parameter :: rngs = 7
integer(4), parameter :: rngb = -1658038656_4 ! X'9D2C5680'
integer(4), parameter :: rngt = 15
integer(4), parameter :: rngc = -272236544_4 ! X'EFC60000'

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine initialize_rng(iseed)

implicit none

! Use the generator from "The Art of Computer Programming: Seminumerical algorithms"
! by Donald Ervin Knuth Vol. 2, 2nd Ed., p. 102, line 25, table 1.

integer(4), intent(in) :: iseed
integer i

allocate(rngstate(0:rngn))

rngstate(0) = iand(iseed,-1)
rngstate(rngn) = 0

do i = 1,rngn-1
    rngstate(i) = int(mod(69069_8*rngstate(i-1),2_8**32),4)
enddo

return

end subroutine initialize_rng

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

real(8) function uniform_dev()

implicit none

integer k, y

if (rngstate(rngn) == 0) then ! generate the next rngn words

    do k = 0,rngn-1
        ! Computing x^u_(i)|x^l_(i+1)        
        y = ior(iand(rngstate(k),rngup),iand(rngstate(mod(k+1,rngn)),rnglo))
        ! Multiplying A
        rngstate(k) = ieor(ieor(rngstate(mod(k+rngm,rngn)),ishft(y,-1)),rngaaugm(iand(y,1)))
    enddo
    
    continue

endif

! Calculate x[i]T
y = rngstate(rngstate(rngn))
y = ieor(y,ishft(y,-rngu))
y = ieor(y,iand(ishft(y,rngs),rngb))
y = ieor(y,iand(ishft(y,rngt),rngc))
y = ieor(y,ishft(y,-rngl))

rngstate(rngn) = mod(rngstate(rngn)+1,rngn)

! Return a uniformly distributed random number
if (y < 0) then
    uniform_dev = (dble(y) + 2.d0**32)/(2.d0**32 - 1.d0)
else
    uniform_dev = dble(y)/(2.d0**32 - 1.d0)
endif

return

end function uniform_dev

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module mt19937_module
