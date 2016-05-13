! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module parser_module

use constants_module, only: iwrite

! the only constant used from constants_module is the generic input-output unit
! number iwrite for the case the errors of the parser_module need to be reported:

integer, parameter :: parser_module_geniounit = iwrite

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine getrecord(iunit,recinput,irec,io)

implicit none

integer iunit, irec, io, i

character(*) recinput

! Gets record from unit iunit.

read(iunit,'(a' // int2str(len(recinput)) // ')', iostat = io) recinput

if (io < 0) then ! cannot read past end of file
    recinput = ' ' ! nullify recinput
    return         ! and return WITHOUT advancing the irec record counter
endif

irec = irec + 1

! To increase compatibility between different systems (in which CR or LF or
! both characters are used as end-of-record characters), replace "weird"
! characters with space in the record just read.

do i = 1,len(recinput)
    if (recinput(i:i) == char(9)) then ! TAB character
        recinput(i:i) = ' '
    endif
    if (recinput(i:i) == char(10)) then ! LF character
        recinput(i:i) = ' '
    endif
    if (recinput(i:i) == char(11)) then ! HT character
        recinput(i:i) = ' '
    endif
    if (recinput(i:i) == char(12)) then ! FF character
        recinput(i:i) = ' '
    endif
    if (recinput(i:i) == char(13)) then ! CR character
        recinput(i:i) = ' '
    endif
enddo

return

end subroutine getrecord

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine break_words(stringin,delimiter,comment,nwords,iw)

implicit none

integer i, nwords, strln, iw(:)

logical inword

character(1) delimiter, comment
character(*) stringin

! Breaks a string into substrings (words) using the specified single character delimiter.
! Anything after the comment character is ignored. nwords is the number of words in the 
! original string. Array iw contains the indexes in the original string where the 
! initial and final letter of each word are.

nwords = 0
inword = .false.

if (stringin(1:1) == comment) then ! record is a comment
    return
endif

if (stringin(1:1) /= delimiter) then ! record begins with a word
    nwords = nwords + 1
    iw(1) = 1
    inword = .true.
endif

strln = len(stringin)

do i = 2,strln

    if (stringin(i:i) == comment) then ! rest of record is a comment
        if (inword) then
            iw(2*nwords) = i-1
        endif
        return
    endif

    if (stringin(i:i) == delimiter .and. &
        stringin(i-1:i-1) /= delimiter) then ! end of word
        inword = .false.
        iw(2*nwords) = i-1
    endif

    if (stringin(i:i) /= delimiter .and. &
        stringin(i-1:i-1) == delimiter) then ! beginning of word
        inword = .true.
        nwords = nwords + 1
        iw(2*nwords-1) = i
    endif

enddo

if (inword) then
    iw(2*nwords) = i-1
endif

return

end subroutine break_words

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

logical function str_is_integer(stringin,ntype)

integer i, ntype

character(*) stringin

logical isnum, numflag

! Check if the string represents an integer. We assume that all spaces (delimiters)
! have been stripped. 

select case (ntype)
    case (4)
        if (len(stringin) > 11) then
            isnum = .false.
        endif
    case (8)
        if (len(stringin) > 20) then
            isnum = .false.
        endif
    case default
        call parser_module_print_error('bad argument ntype in str_is_integer')
end select

numflag = .false.

do i = 1,len(stringin)

    select case(stringin(i:i))
        case('-','+')
            if (i /= 1) then
                isnum = .false.
                exit
            endif
        case('0':'9')
            numflag = .true.
            isnum = .true.
        case default
            isnum = .false.
            exit
    end select
    
enddo

if (numflag) then
    str_is_integer = isnum
else
    str_is_integer = .false.
endif    

return

end function str_is_integer

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

logical function str_is_real(stringin,ntype)

integer i, ntype

character(*) stringin
character(1) decs

logical isnum, decpflag, expflag, exp0flag, numflag

select case (ntype)
    case (4)
        decs = 'e'
    case (8)
        decs = 'd'
    case default
        call parser_module_print_error('bad argument ntype in str_is_real')
end select

decpflag = .false.
expflag = .false.
exp0flag = .false.
numflag = .false.

do i = 1,len(stringin)

    select case(stringin(i:i))
        case('-','+')
            if (.not.(i == 1 .or. exp0flag)) then ! signs allowed only in the beginning of the mantissa or the exponent
                isnum = .false.
                exit
            endif
            if (exp0flag) exp0flag = .false.
       case('.')
            if (decpflag .or. expflag) then ! only one decimal place allowed in the mantissa
                isnum = .false.
                exit
            endif
            decpflag = .true.
       case('E','e') ! only one E, e, D, or d exponent symbol allowed
            if (expflag) then
                isnum = .false.
                exit
            endif
! If one wants to be strict and use only E for real(4), and E or D for real(8)
! one should uncomment this if-block
!            if (decs == 'd') then
!                isnum = .false.
!                exit
!            endif           
            expflag = .true. ! flag that we are in the exponent
            exp0flag = .true. ! this flag will be set to false after the first exponent character is read
       case('D','d')
            if (expflag) then
                isnum = .false.
                exit
            endif
            if (decs == 'e') then
                isnum = .false.
                exit
            endif           
            expflag = .true. ! flag that we are in the exponent
            exp0flag = .true. ! this flag will be set to false after the first exponent character is read
        case('0':'9')
            if (exp0flag) exp0flag = .false.
            numflag = .true.
            isnum = .true.
        case default
            isnum = .false.
            exit
    end select
    
enddo

if ((decpflag .or. expflag) .and. (numflag)) then
    str_is_real = isnum
else
    str_is_real = .false.
endif

return

end function str_is_real

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

logical function str_check_expression(stringin,nwords,iw,patrn,strictnwords)

character(*) stringin, patrn

integer i, j, m, nwords, nstr, iw(:), nwordsp, iwp(200)

logical strictnwords

! Checks a string stringin which has already been parsed by the break_words subroutine
! (so nwords and iw are the number of words and indexes of stringin), to see whether
! it matches the pattern patrn. The latter is a string say '2AA/3I4/1R4' which specifies
! that stringin should contain 2 strings, 3 integers kind 4, and a real type 4. The specifiers
! are two letter strings: I4, I8, R4, R8, AA for integer(4), integer(8), real(4), real(8), and 
! string, respectively. Anything before these specifiers is the number of times the specifier
! should be repeated.
! Flag strictnwords specifies whether the stringin and the patrn are required to have the same 
! number of words.

str_check_expression = .true.

call break_words(patrn,'/','#',nwordsp,iwp)

! If strict word-number checking is required, calculate the number of words in the pattern
! and compare it to that of the input
if (strictnwords) then
    m = 0
    do i = 1,nwordsp
        if (iwp(2*i)-iwp(2*i-1) > 1) then ! parse number of repetitions
            if (.not.str_is_integer(patrn(iwp(2*i-1):iwp(2*i)-2),4)) then
                call parser_module_print_error('bad patrn argument in str_check_expression')
            else
                m = m + str2int(patrn(iwp(2*i-1):iwp(2*i)-2))            
            endif
        else
            m = m + 1             
        endif
    enddo
    if (m /= nwords) then
        str_check_expression = .false.
        return
    endif
endif    

nstr = 0
do i = 1,nwordsp
    
    m = 1
    if (iwp(2*i)-iwp(2*i-1) > 1) then ! parse number of repetitions
    
        if (.not.str_is_integer(patrn(iwp(2*i-1):iwp(2*i)-2),4)) then
            call parser_module_print_error('bad patrn argument in str_check_expression')
        endif
         
         m = str2int(patrn(iwp(2*i-1):iwp(2*i)-2))
         
    endif
    
    do j = 1,m
        
        nstr = nstr + 1
        
        if (nstr > nwords) then 
            str_check_expression = .false.
            return
        endif
        
        select case (patrn(iwp(2*i)-1:iwp(2*i)))
            case ('AA')
                continue
            case ('I4')
                if (.not.str_is_integer(stringin(iw(2*nstr-1):iw(2*nstr)),4)) then
                    str_check_expression = .false.
                    return
                endif
            case ('I8')
                if (.not.str_is_integer(stringin(iw(2*nstr-1):iw(2*nstr)),8)) then
                    str_check_expression = .false.
                    return
                endif
            case ('R4')
                if (.not.str_is_real(stringin(iw(2*nstr-1):iw(2*nstr)),4)) then
                    str_check_expression = .false.
                    return
                endif
            case ('R8')
                if (.not.str_is_real(stringin(iw(2*nstr-1):iw(2*nstr)),8)) then
                    str_check_expression = .false.
                    return
                endif
        end select                    
    
    enddo
    
enddo

return

end function str_check_expression

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine str_blank(stringin)

implicit none

integer i

character(*) stringin

do i = 1,len(stringin)
    stringin(i:i) = ' '
enddo

end subroutine str_blank

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine lowercase(stringin)

implicit none

integer i

character(*) stringin

do i = 1,len(stringin)
    
    if(stringin(i:i) == 'A') then
        stringin(i:i)='a'
    elseif (stringin(i:i) == 'B') then
        stringin(i:i)='b'
    elseif (stringin(i:i) == 'C') then
        stringin(i:i)='c'
    elseif (stringin(i:i) == 'D') then
        stringin(i:i)='d'
    elseif (stringin(i:i) == 'E') then
        stringin(i:i)='e'
    elseif (stringin(i:i) == 'F') then
        stringin(i:i)='f'
    elseif (stringin(i:i) == 'G') then
        stringin(i:i)='g'
    elseif (stringin(i:i) == 'H') then
        stringin(i:i)='h'
    elseif (stringin(i:i) == 'I') then
        stringin(i:i)='i'
    elseif (stringin(i:i) == 'J') then
        stringin(i:i)='j'
    elseif (stringin(i:i) == 'K') then
        stringin(i:i)='k'
    elseif (stringin(i:i) == 'L') then
        stringin(i:i)='l'
    elseif (stringin(i:i) == 'M') then
        stringin(i:i)='m'
    elseif (stringin(i:i) == 'N') then
        stringin(i:i)='n'
    elseif (stringin(i:i) == 'O') then
        stringin(i:i)='o'
    elseif (stringin(i:i) == 'P') then
        stringin(i:i)='p'
    elseif (stringin(i:i) == 'Q') then
        stringin(i:i)='q'
    elseif (stringin(i:i) == 'R') then
        stringin(i:i)='r'
    elseif (stringin(i:i) == 'S') then
        stringin(i:i)='s'
    elseif (stringin(i:i) == 'T') then
        stringin(i:i)='t'
    elseif (stringin(i:i) == 'U') then
        stringin(i:i)='u'
    elseif (stringin(i:i) == 'V') then
        stringin(i:i)='v'
    elseif (stringin(i:i) == 'W') then
        stringin(i:i)='w'
    elseif (stringin(i:i) == 'X') then
        stringin(i:i)='x'
    elseif (stringin(i:i) == 'Y') then
        stringin(i:i)='y'
    elseif (stringin(i:i) == 'Z') then
        stringin(i:i)='z'
    endif
    
enddo

return
end subroutine lowercase

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine strfind(stringin,token,ntimes,iw)

implicit none

integer i, ntimes, strln, tokln, iw(:)

character(*) stringin, token

strln = len(stringin)
tokln = len(token)

ntimes = 0

if (tokln > strln) then

    return

elseif (tokln == strln) then

    if (stringin == token) then
        ntimes = 1
        iw(1) = 1
    endif

    return

else

    do i = 1,strln-tokln+1
        if (stringin(i:tokln+i-1) == token) then
            ntimes = ntimes + 1
            iw(ntimes) = i
        endif
    enddo
        
endif    

end subroutine strfind

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

logical function striccompare(stringin1,stringin2)

implicit none

integer strln1, strln2

character(*) stringin1, stringin2

strln1 = len(stringin1)
strln2 = len(stringin2)

striccompare = .false.

call lowercase(stringin1)
call lowercase(stringin2)

if (strln1 == strln2) then
    if (stringin1 == stringin2) then
        striccompare = .true.
    endif
endif

return

end function striccompare

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

character(32) function int2str(i)

implicit none

integer i

! DEVELOPER COMMENT MICHSTAM: the formatted version of the following command
! is incompatible with fortran standards older than F2003 (for instanse it 
! will create a runtime error in pgf90). The unformatted version is always fine 
! but it may reduce the precision of the number converted to a string. Thus, 
! this construct should not be used when saving the state of a program in order 
! to resume the simulation at a later time.
write(int2str,*) i
!write(int2str,'(I32)') i

call remove_leading_blanks(int2str)

end function int2str

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

character(32) function int82str(i)

implicit none

integer(8) i

! DEVELOPER COMMENT MICHSTAM: the formatted version of the following command
! is incompatible with fortran standards older than F2003 (for instanse it 
! will create a runtime error in pgf90). The unformatted version is always fine 
! but it may reduce the precision of the number converted to a string. Thus, 
! this construct should not be used when saving the state of a program in order 
! to resume the simulation at a later time.
write(int82str,*) i
!write(int82str,'(I32)') i

call remove_leading_blanks(int82str)

end function int82str

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

character(32) function dbl2str(x)

implicit none

real(8) x

! DEVELOPER COMMENT MICHSTAM: the formatted version of the following command
! is incompatible with fortran standards older than F2003 (for instanse it 
! will create a runtime error in pgf90). The unformatted version is always fine 
! but it may reduce the precision of the number converted to a string. Thus, 
! this construct should not be used when saving the state of a program in order 
! to resume the simulation at a later time.
write(dbl2str,*) x
!write(dbl2str,'(ES32.16E3)') x

call remove_leading_blanks(dbl2str)

end function dbl2str

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

character(32) function real2str(x)

implicit none

real(4) x

! DEVELOPER COMMENT MICHSTAM: the formatted version of the following command
! is incompatible with fortran standards older than F2003 (for instanse it 
! will create a runtime error in pgf90). The unformatted version is always fine 
! but it may reduce the precision of the number converted to a string. Thus, 
! this construct should not be used when saving the state of a program in order 
! to resume the simulation at a later time.
write(real2str,*) x
!write(real2str,'(ES32.8E3)') x

call remove_leading_blanks(real2str)

end function real2str

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function str2int(stringin)

implicit none

character(*) stringin

read(stringin,*) str2int

end function str2int

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer(8) function str2int8(stringin)

implicit none

character(*) stringin

read(stringin,*) str2int8

end function str2int8

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

real(8) function str2dbl(stringin)

implicit none

character(*) stringin

read(stringin,*) str2dbl

end function str2dbl

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function findinstrv(stringin,m0,m1,m2,strcmp)

implicit none

integer i, m0, m1, m2
character(*) stringin(m0:)
character(*) strcmp

findinstrv = -1

do i = m1,m2
    if (trim(stringin(i)) == trim(strcmp)) then
        findinstrv = i
        return
    endif
enddo

end function findinstrv

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine remove_leading_blanks(stringin)

implicit none

integer i, j, strln

character(*) stringin

if (stringin(1:1) /= ' ') then
    return
endif

strln = len(stringin)

do i = 2,strln
    if (stringin(i:i) /= ' ') then
        exit
    endif
enddo

stringin(1:strln-i+1) = stringin(i:strln)

do j = strln-i+2,strln
    stringin(j:j) = ' '
enddo

end subroutine remove_leading_blanks

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine parser_module_print_error(errdesc)

integer iwrite

character(*) errdesc

iwrite = parser_module_geniounit

write(iwrite,'(/,a)') '***************'
write(iwrite,'(/,a)') 'Error message returned from parser_module: ' // errdesc //'.'
write(iwrite,'(/,a)') '***************'
stop

end subroutine parser_module_print_error

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module parser_module
