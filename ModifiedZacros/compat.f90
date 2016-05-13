! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

! The following is for compatibility with Compaq Fortran
! The following nonsense definition of function setvbuf3f 
! should be commented out when running in PGI
function setvbuf3f(lu, typ, size)
integer setvbuf3f, lu, typ, size
setvbuf3f=0
end