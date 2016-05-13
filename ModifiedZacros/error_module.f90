! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module error_module

use constants_module
use parser_module

implicit none      

integer, allocatable :: warncounters(:)

character(1024) moreinfo

contains

  subroutine cleanup_error_module
    ! Cleanup error globals
    if(allocated(warncounters)) deallocate(warncounters)
  end subroutine cleanup_error_module
      
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine error(ierror)

implicit none

integer ierror

write(iwrite,'(/,a)') '***************'

select case (ierror)
    case (0)
        continue
        
    ! *** General errors
    case (1)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': at least one non-optional input file is missing.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': restart file appears to exist ' // &
                              'but it cannot be opened or is corrupted.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)

    case (100)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid duplicate specification.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)

    ! *** Simulation setup parsing errors
    
    case (1001)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': random_seed keyword must be followed by a single integer in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1002)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': temperature keyword must be followed by a single real in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1003)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': pressure keyword must be followed by a single real in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1004)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': n_gas_specs keyword must be followed by a single integer in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1005)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': gas_specs_names keyword must be preceded by n_gas_specs keyword in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1006)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': gas_specs_names keyword must be followed by n_gas_specs strings in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1007)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': gas_molec_weights keyword must be preceded by n_gas_specs keyword in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1008)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': gas_molec_weights keyword must be followed by n_gas_specs reals in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1009)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': gas_molec_fracs keyword must be preceded by n_gas_specs keyword in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1010)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': gas_molec_fracs keyword must be followed by n_gas_specs reals in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1011)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': n_surf_specs keyword must be followed by a single integer in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1012)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': surf_specs_names keyword must be preceded by n_surf_specs keyword in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1013)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': surf_specs_names keyword must be followed by n_surf_specs strings in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1014)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': surf_specs_dent keyword must be preceded by n_surf_specs keyword in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1015)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': surf_specs_dent keyword must be followed by n_surf_specs integers in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1016)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': maxsteps keyword must be followed by a single integer or the word "infinity" in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1017)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': maxtime keyword must be followed by a single real or the word "infinity" in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1018)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': snapshots keyword can only be followed by: off; on event; on time [Dt(real8)]; in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1019)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': process_statistics keyword can only be followed by: off; on event; on time [Dt(real8)]; in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1020)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': species_numbers keyword can only be followed by: off; on event; on elemevent [i(int4)]; on time [Dt(real8)]; on logtime [t0(real8) c(real8)]; in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1021)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': event_report keyword must be followed by either "on" or "off" in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1022)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': wall_time keyword must be followed by a single integer in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1023)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': no_restart keyword cannot followed by any directive in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1024)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': debug_report_processes keyword cannot followed by any directive in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1025)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': debug_check_processes keyword cannot followed by any directive in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1026)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': debug_check_lattice keyword cannot followed by any directive in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1027)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': debug_newtons_method keyword cannot followed by any directive in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1028)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': debug_report_global_energetics keyword cannot followed by any directive in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1029)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': debug_check_global_energetics keyword cannot followed by any directive in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1030)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': gas_molec_weights keyword must be preceded by n_gas_specs keyword in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1031)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': gas_molec_weights keyword must be followed by n_gas_specs reals in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
        
        
    case (1100)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid simulation specification expression in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
!    case (1101)
!        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': random seed has not been specified in ' // trim(csimfname) // '.'
    case (1102)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid random seed equal to zero has been specified in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1103)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': temperature has not been specified in ' // trim(csimfname) // '.'
    case (1104)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid non-positive temperature has been specified in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1105)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': pressure has not been specified in ' // trim(csimfname) // '.'
    case (1106)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid negative pressure has been specified in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1107)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': number of gas species has not been specified in ' // trim(csimfname) // '.'
    case (1108)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid negative number of gas species has been specified in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1109)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': gas species molecular weights have not been specified in ' // trim(csimfname) // '.'
    case (1110)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid non-positive gas species molecular weight has been specified in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1111)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': gas species molar fractions have not been specified in ' // trim(csimfname) // '.'
    case (1112)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid gas species molar fraction, outside the range of [0,1] has been specified in ' // &
                              trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1113)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': number of surface species has not been specified in ' // trim(csimfname) // '.'
    case (1114)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid non-positive number of surface species has been specified in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1115)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': surface species dentations have not been specified in ' // trim(csimfname) // '.'
    case (1116)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid non-positive surface species dentation has been specified in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1117)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': n_gauss_pts keyword must be followed by a positive integer in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1118)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid non-positive number of points for Gauss quadrature has been specified in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1119)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': k_newton_max keyword must be followed by a positive integer in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1120)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid non-positive number of maximum Newton iterations has been specified in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1121)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': tol_dx_newton keyword must be followed by a (small) positive real in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1122)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid non-positive number of Newton Dx tolerance has been specified in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1123)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': tol_rhs_newton keyword must be followed by a (small) positive real in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1124)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid non-positive number of Newton RHS tolerance has been specified in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1125)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': gas molar fractions sum up to a number greater than 1 in ' // trim(csimfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)

    ! *** Lattice setup parsing errors

    case (2001)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': file ' // trim(clatfname) // ' must start with the keyword lattice followed by ' // &
                              'one of the following directives: default_choice, periodic, or explicit.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2002)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': lattice has already been specified in ' // trim(clatfname) // '. ' // &
                              'Expected end_lattice keyword.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2003)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid or absent default lattice specification in ' // trim(clatfname) // '. ' // &
                              'A valid specification consists of one of the following keywords: triangular_periodic, rectangular_periodic, or hexagonal_periodic; ' // &
                              'followed by a real and two integers.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2030)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': cell_vectors keyword cannot be followed by any directive in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2031)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': two lines following the cell_vectors keyword must contain two reals each, specifying ' // &
                              'the unit cell vectors in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2032)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': repeat_cell keyword must be followed by two integers in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2033)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': n_site_types keyword must be followed by a single integer in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2034)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site_type_names keyword must be preceded by n_site_types keyword in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2035)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site_type_names keyword must be followed by n_site_types strings in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2036)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': n_cell_sites keyword must be followed by a single integer in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2037)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site_types keyword must be preceded by n_cell_sites keyword in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2038)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site_types keyword must be preceded by n_site_types keyword in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2039)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site_types keyword must be followed by n_cell_sites integers or site type names ' // &
                              '(if already defined) in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2040)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid site type or site type name following site_types keyword in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)        
    case (2041)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site_coordinates keyword must be preceded by n_cell_sites keyword in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)        
    case (2042)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': n_cell_sites lines following the cell_vectors keyword must contain two reals each, specifying ' // &
                              'each site''s coordinates in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2043)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': neighboring_structure keyword must be preceded by n_cell_sites keyword in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)        
    case (2044)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid neighboring structure specification in ' // trim(clatfname) // '.'
        write(iwrite,'(a)') 'Check spacing, make sure that site numbers must be integers, and keep in mind that only the self, north, northeast, east, and ' // &
                            'southeast neighbors have to be specified. Anything else is invalid.'
        write(iwrite,'(a)') '  Right examples: 1-2 self  or  2-1 southeast'
        write(iwrite,'(a)') '  Wrong examples: 1 - 2 self  or  1- 2 self'
        write(iwrite,'(/,a)') 'More information:'
        write(iwrite,'(a)') trim(moreinfo)
    case (2045)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': incomplete lattice specification in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)        
    case (2060)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': n_sites keyword must be followed by a single integer in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2061)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': max_coord keyword must be followed by a single integer in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2062)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': lattice_structure keyword must be preceded by n_sites, max_coord, n_site_types ' // &
                              '(and possibly site_type_names) in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2063)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site number out of bounds in ' // trim(clatfname) // '.' // &
                              'Check numbering of site and/or neighbors, with respect to max_coord and n_sites.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2064)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site type names have not been defined; yet, the site types are specified by strings in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2065)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': duplicate site definition in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2066)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': for explicit lattice definition, each line in the lattice_structure block in ' // trim(clatfname) // ' must contain: '
        write(iwrite,'(a)')   '  site_number(int4) xcoord(real8) ycoord(real8) site_type(int4 or string) nneigbors(int4) neigh_1(int4) neigh_2(int4) ... neigh_nneigbors(int4)'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2067)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': missing information for at least one site in ' // trim(clatfname) // '.'
    case (2100)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid lattice specification expression in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2101)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': multiple encounters of a site in neighboring arrays in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2102)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid parity of neighboring arrays in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2105)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': absence of site-site neighboring parity in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2106)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': encountered end-of-file prematurely in ' // trim(clatfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)

    ! *** Mechanism setup parsing errors
    
    case (3001)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': file ' // trim(cmechfname) // ' must start with the mechanism keyword ' // &
                              'followed by blocks specifying elementary steps.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3002)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': step keyword can be followed by a single string in which no spaces are allowed.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3003)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid specification of gas products in ' // trim(cmechfname) // ': ' // &
                              'keyword gas_reacs_prods must be followed by pairs of a gas species name and the stoichiometric coefficient.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3004)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': reference to an unknown species in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3005)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': sites keyword must be followed by a single integer in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3006)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': initial/final keywords must be preceded by sites keyword in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3007)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': the lines after initial keyword in ' // trim(cmechfname) // ' specify the state ' // &
                              'of each site and must contain: '
        write(iwrite,'(a)')   '  entity_number(int4) species(string) dentate(int4)'        
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3008)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': entity-species information mismatch in state specification in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3009)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': error in elementary step neighboring list in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3010)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': illegal self-neighboring in elementary step specification in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3011)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': duplicate pair in elementary step neighboring list in ' // trim(cmechfname) // &
                              '. Remember, the order does not play role, so for instance 1-2 is the same as 2-1.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3012)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': at least one site in an elementary step specification has no neighbors in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3013)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site_types keyword must be followed by integers or strings specifying the types for all sites ' // &
                              'in an elementary step specification in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3014)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': reference to an unknown site type in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3015)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': pre_expon keyword must be followed by either a single real or seven reals ' // &
                              '(for temperature dependent pre-exponentials) in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3016)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': pe_ratio keyword is invalid for irreversible steps in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3017)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': pe_ratio keyword must be followed by a real in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3018)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': activ_eng keyword must be followed by a real in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3022)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': an irreversible step block must close with keyword ' // &
                              'end_step and a reversible one with end_reversible_step.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3023)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': opening elementary step section before closing the previous one in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3024)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': incomplete information given for an elementary step in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3025)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid numbering of participating entities in an elementary step in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3026)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': incomplete dentate assignment of participating entities in an elementary step in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3027)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': if any one of the following keywords ' // &
                              '[site_types, pre_expon, pe_ratio, activ_eng] ' // &
                              'is encountered outside of a variant block, one can no longer specify variants ' // &
                              'of an elementary step in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3028)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid dentation specification in initial or final state of an elementary step in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3031)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': prox_factor keyword must be followed by a real in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)    
    case (3032)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': supplied proximity factor value is out of range 0.0 - 1.0 in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)    
    case (3033)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': once variants have been specified for an elementary step, the use of keywords ' // &
                              'site_types, pre_expon, pe_ratio, prox_factor, activ_eng, is invalid outside a variant block in file ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)    
    case (3034)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword angles must be preceded by keyword sites in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3035)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword angles must be preceded by keyword neighboring in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3036)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword angles must be followed by expressions each of which specifies the angle ' // &
                              'between three sites in ' // trim(cmechfname) // '. These expressions must follow the syntax Site1-Site2-Site3:AngleValue.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3037)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid sequence of sites in angle specification in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3038)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': neighboring structure in an angle specification sequence of sites is incompatible ' // &
                              'with that elementary step''s defined neighboring structure in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)    
    case (3040)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword absl_orientation must be preceded by keyword sites in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3041)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword absl_orientation must be preceded by keyword neighboring in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3042)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword absl_orientation must be followed by an expression that specifies the angle ' // &
                              'between an edge of the graph and the x-axis unit vector in ' // trim(cmechfname) // '. This expression must follow the syntax Site1-Site2:AngleValue.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3043)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid edge in absolute orientation specification in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3044)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': duplicate site assignment in initial/final ' // &
                              'state of an elementary step in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    
    case (3100)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid mechanism specification expression in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (3101)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': encountered end-of-file prematurely in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
        
        
    ! *** Energetics setup parsing errors
    
    case (8001)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': file ' // trim(cenergfname) // ' must start with the energetics keyword ' // &
                              'followed by blocks specifying clusters.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8002)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': cluster keyword can be followed by a single string in which no spaces are allowed.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8004)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': reference to an unknown species in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8005)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': sites keyword must be followed by a single integer in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8006)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': lattice_state keyword must be preceded by sites keyword in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8007)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': the lines after lattice_state keyword in ' // trim(cenergfname) // ' specify the state ' // &
                              'of each site and must contain: '
        write(iwrite,'(a)')   '  entity_number(int4) species(string) dentate(int4)'        
        write(iwrite,'(a)')   'If you want to declare an unknown/unspecific site state, use the "&" symbol in all three fields.'        
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8008)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': entity-species information mismatch in state specification in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8009)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': error in cluster neighboring list in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8010)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': illegal self-neighboring in cluster specification in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8011)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': duplicate pair in cluster neighboring list in ' // trim(cenergfname) // &
                              '. Remember, the order does not play role, so for instance 1-2 is the same as 2-1.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8012)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': at least one site in an cluster specification has no neighbors in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8013)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site_types keyword must be followed by integers or strings specifying the types for all sites ' // &
                              'in a cluster specification in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8014)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': reference to an unknown site type in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8018)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': cluster_eng keyword must be followed by a real in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8023)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': opening cluster section before closing the previous one in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8024)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': incomplete information given for a cluster in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8025)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid numbering of participating entities in a cluster in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8026)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': incomplete dentate assignment of participating entities in a cluster in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8027)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': if any one of the following keywords ' // &
                              '[site_types, cluster_eng] ' // &
                              'is encountered outside of a variant block, one can no longer specify variants ' // &
                              'of a cluster in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8028)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid dentation specification in lattice state of a cluster in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8029)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': duplicate site assignment in lattice state of a cluster in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8031)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword site_types must be preceded by keyword sites in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8032)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword angles must be preceded by keyword sites in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8033)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword angles must be preceded by keyword neighboring in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8034)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword angles must be followed by expressions each of which specifies the angle ' // &
                              'between three sites in ' // trim(cenergfname) // '. These expressions must follow the syntax Site1-Site2-Site3:AngleValue.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8035)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid sequence of sites in angle specification in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8036)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': neighboring structure in an angle specification sequence of sites is incompatible ' // &
                              'with that cluster''s defined neighboring structure in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
   
    case (8050)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword absl_orientation must be preceded by keyword sites in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8051)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword absl_orientation must be preceded by keyword neighboring in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8052)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': keyword absl_orientation must be followed by an expression that specifies the angle ' // &
                              'between an edge of the graph and the x-axis unit vector in ' // trim(cenergfname) // '. This expression must follow the syntax Site1-Site2:AngleValue.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8053)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid edge in absolute orientation specification in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8054)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': once variants have been specified for a cluster, the use of keywords site_types, ' // &
                              'graph_multiplicity, cluster_eng, angles, no_mirror_images, absl_orientation, is invalid outside a variant block in file ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)    

    case (8100)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid energetics specification expression in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (8101)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': encountered end-of-file prematurely in ' // trim(cenergfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)



    ! *** Initial state setup parsing errors
    
    case (4001)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': file ' // trim(cstatefname) // ' must start with the initial_state keyword ' // &
                              'followed by blocks or single lines containing species seeding instructions.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4002)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': reference to an unknown species in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4003)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site number out of bounds in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4004)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': seed_on_sites keyword must be followed by a species name (string) and ' // &
                              'nspecdent integers in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4005)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': seed_multiple keyword must be followed by a species name (string) and ' // &
                              'a positive integer in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4006)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': neighboring keyword is invalid for species with a single dentate in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4007)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': error in seed_multiple instruction neighboring list in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4008)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': illegal self-neighboring in elementary step specification in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4009)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': duplicate pair in elementary step neighboring list in ' // trim(cstatefname) // &
                              '. Remember, the order does not play role, so for instance 1-2 is the same as 2-1.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4010)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': site_types keyword must be followed by integers or strings specifying the types for all sites ' // &
                              'in an seed_multiple instruction in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4011)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': reference to an unknown site type in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4012)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': opening species seeding instruction section before closing the previous one in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4013)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': incomplete information given for a seeding instruction in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)


    case (4100)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid species seeding instruction in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (4101)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': dentate neighboring has not been specified in ' // trim(cstatefname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)

    case (4501)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': seeding of multiple adsorbates on the surface cannot continue.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
        
    ! *** KMC simulation runtime errors

    case (5001)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid parity of adsorbspecposi and latticestate arrays.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5002)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': cannot remove entity until all processes and clusters in which ' // &
                              'it is involved are removed first.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5003)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': a process involves an adsorbate with an incorrect mapping ' // &
                               'between elementary step sites and lattice sites.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5004)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': a process involves an adsorbate which does not list that process ' // &
                              'in its participation list.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5005)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': an adsorbate lists a process that contains an incorrect mapping ' // &
                              'between elementary step sites and lattice sites.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5006)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': an adsorbate lists a process which does not correctly point back to that adsorbate.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5007)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': an adsorbate lists a process which lists that adsorbate as a different species.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5008)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': an adsorbate lists a process which lists that adsorbate with wrong dentation.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5009)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': invalid entity has been specified for removal.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5010)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': attempted to add an entity on an invalid lattice site.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5011)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': Newton''s method failed to find the time for the occurrence of an elementary event.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)

    case (5013)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': a cluster involves an adsorbate with an incorrect mapping ' // &
                               'between elementary step sites and lattice sites.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5014)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': a cluster involves an adsorbate which does not list that process ' // &
                              'in its participation list.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5015)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': an adsorbate lists a cluster that contains an incorrect mapping ' // &
                              'between elementary step sites and lattice sites.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5016)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': an adsorbate lists a cluster which does not correctly point back to that adsorbate.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5017)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': an adsorbate lists a cluster which lists that adsorbate as a different species.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5018)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': an adsorbate lists a cluster which lists that adsorbate with wrong dentation.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (5020)
        write(iwrite,'(/,a)') 'Error code ' // trim(int2str(ierror)) // ': upon executing a process the global energy after the process ' // &
                              'was not equal to energy before plus the delta-energy of that process.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
        

    case default
        write(iwrite,'(/,a)') 'Error - Unspecified error code.'      
       
end select

write(iwrite,'(/,a)') '***************'

write(iwrite,'(/,a)') '> ABNORMAL TERMINATION DUE TO FATAL ERROR <'
    
close(iread)
close(ilatread)
close(imechread)
close(iwrite)
close(ihistory)
close(iprocdbg)

stop
    
end subroutine error

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine initialize_warning_counters()

implicit none

allocate(warncounters(1:10))

warncounters = 0

return

end subroutine initialize_warning_counters

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine warning(iwarn)

implicit none

integer iwarn

select case (iwarn)
    case (0)
        continue
        
    ! *** Warnings
    case (1)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': lost track of clock time ' // &
                              '(processor clock possibly messed-up). Restart data will be written now and execution will terminate.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)

    case (2)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': no more processes can occur ' // &
                              '(this may indicate that the surface is poisoned). The program will terminate.'
    case (3)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': restart information could not be written.' // &
                              ' You will not be able to resume this simulation from the point it stopped.'
    case (1001)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': more than n_gas_specs gas species names appear in ' // &
                              trim(csimfname) // '. Only the first n_gas_specs will be read.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1002)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': more than n_gas_specs gas species molecular weights' // &
                              ' appear in ' // trim(csimfname) // '. Only the first n_gas_specs will be read.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1003)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': more than n_gas_specs gas species molar fractions' // &
                              ' appear in ' // trim(csimfname) // '. Only the first n_gas_specs will be read.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1004)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': more than n_surf_specs surface species names appear in ' // &
                              trim(csimfname) // '. Only the first n_surf_specs will be read.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1005)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': more than n_surf_specs surface species dentations ' // &
                              'appear in ' // trim(csimfname) // '. Only the first n_surf_specs will be read.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (1100)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': read past the end of file ' // trim(csimfname) // '. ' // &
                              'It is advisable to end the simulation input with the finish keyword.'

    case (2001)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': more than n_site_types site type names names appear in ' // &
                              trim(clatfname) // '. Only the first n_site_types will be read.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2100)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': read past the end of file ' // trim(clatfname) // '. ' // &
                              'It is advisable to end the lattice input with the end_lattice keyword.'
        write(iwrite,'(/,a)') '***************'
    case (2103)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': self-neighboring in ' // trim(clatfname) // '. This can be normal for small lattices.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
    case (2104)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': possibly invalid lattice neighbor list in ' // trim(clatfname) // '. This can be normal for small lattices.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)

    case (3029)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': negative activation energy specified for an elementary step in ' // trim(cmechfname) // '.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)

    case (3100)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': read past the end of file ' // trim(cmechfname) // '. ' // &
                              'It is advisable to end the mechanism input with the end_mechanism keyword.'

    case (4100)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': read past the end of file ' // trim(cstatefname) // '. ' // &
                              'It is advisable to end the initial state input with the end_initial_state keyword.'

    case (8100)
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': read past the end of file ' // trim(cenergfname) // '. ' // &
                              'It is advisable to end the energetics input with the end_energetics keyword.'

    case (10006)
        
        if (warncounters(1) == maxrepwarnings) then
            return
        endif
    
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': Newton-Raphson may have not converged in subroutine' // &
                              ' calculate_elemstep_rate. The subroutine tried to find the time when the next reaction event will' // &
                              ' occur for the specified temperature ramp.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
        
        warncounters(1) = warncounters(1) + 1
        
        if (warncounters(1) == maxrepwarnings) then
            write(iwrite,'(/,a)') 'Further warning code ' // trim(int2str(iwarn)) // ' messages will be suppressed during this calculation...'
        endif

    case (10007)
        
        if (warncounters(2) == maxrepwarnings) then
            return
        endif
        
        write(iwrite,'(/,a)') '***************'
        write(iwrite,'(/,a)') 'Warning code ' // trim(int2str(iwarn)) // ': Newton-Raphson has detected that a propensity tends to zero' // &
                              ' for long times and yet a solution was not found in subroutine calculate_elemstep_rate.' // &
                              ' The next reaction event time will be set to a huge value. This behavior is to be expected for simulated' // &
                              ' annealing calculations, otherwise the issue needs further investigation.'
        write(iwrite,'(/,a)') 'More information: '
        write(iwrite,'(a)') trim(moreinfo)
     
        warncounters(2) = warncounters(2) + 1
        
        if (warncounters(2) == maxrepwarnings) then
            write(iwrite,'(/,a)') 'Further warning code ' // trim(int2str(iwarn)) // ' messages will be suppressed during this calculation...'
        endif

end select    

write(iwrite,'(/,a)') '***************'

return

end subroutine warning

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module error_module
