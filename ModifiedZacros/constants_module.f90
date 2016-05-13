! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module constants_module

implicit none

! every how many seconds a restart file will be written
integer, parameter :: dtwrestart = 3600

! for repetitive warnings: the maximum number of message that will be issued for each such warning
integer, parameter :: maxrepwarnings = 10

! *** Mathematical constants

! pi value
real(8), parameter :: pi = 3.141592653589793d0
! e value
real(8), parameter :: e = 2.718281828459046d0

! *** Physical constants (in SI or SI-derived units)

! ~ Speed of light
real(8), parameter :: clight = 2.99792458d+08 ! m/s

! ~ Avogadro's number
real(8), parameter :: navog = 6.02214179d+23 ! 1/mol

! ~ Gas constant
real(8), parameter :: rgas = 8.314472d+00   ! J/K/mol

! ~ Plank's constant
real(8), parameter :: hplank = 6.62606896d-34 ! J*s

! ~ Boltzmann's constant
real(8), parameter :: kboltz = rgas/navog ! J/K

! *** Conversion factors

! ~ Energy: 1 J = engconv User_units
!real(8), parameter :: enrgconv = 1.d0 ! J
!real(8), parameter :: enrgconv = 1.d-3*navog ! kJ/mol
real(8), parameter :: enrgconv = 6.24150974d+18 ! eV
!real(8), parameter :: enrgconv = 1.d+07 ! erg
!real(8), parameter :: enrgconv = 2.3901d-4 ! kcal
!real(8), parameter :: enrgconv = 2.3901d-4*navog ! kcal/mol
!real(8), parameter :: enrgconv = 2.3901d-1 ! cal

! ~ Distanse: 1 m = engconv User_units
!real(8), parameter :: distconv = 1.d0 ! m
!real(8), parameter :: distconv = 1.d+02 ! cm
real(8), parameter :: distconv = 1.d+10 ! Angstrom

! ~ Time: 1 s = engconv User_units
real(8), parameter :: timeconv = 1.d0 ! s
!real(8), parameter :: timeconv = 2.7777777777777778d-4 ! h
!real(8), parameter :: timeconv = 1.6666666666666667d-2 ! min
!real(8), parameter :: timeconv = 1.d+09 ! ns
!real(8), parameter :: timeconv = 1.d+12 ! ps

! ~ Mass:
real(8), parameter :: massconv = 1.d+0 ! kg

! ~ Pressure:
!real(8), parameter :: presconv = 1.d+0 ! Pa
real(8), parameter :: presconv = 1.d-5 ! bar
!real(8), parameter :: presconv = 9.8692d-06 ! atm

! *** KMC code parameters

real(8), parameter :: angletol = 1e-3

! commenting character
character(1), parameter :: remchar = '#'

! double the maximum number of words separated by a delimiter that can be read 
integer, parameter :: maxwords2 = 6000

! record length (in number of characters) that can be read 
integer, parameter :: lengthrecinp = 8192

! allowed length (in characters) of the species', site-types' and mechanism-steps' names
integer, parameter :: nnam0 = 64

! *** File names and unit numbers

! simulation input file
character(30), parameter :: csimfname = 'simulation_input.dat'
integer, parameter :: iread = 101

! lattice specification input file name and unit number
character(30), parameter :: clatfname = 'lattice_input.dat'
integer, parameter :: ilatread = 102

! mechanism specification input file name and unit number
character(30), parameter :: cmechfname = 'mechanism_input.dat'
integer, parameter :: imechread = 103

! energetics specification input file name and unit number
character(30), parameter :: cenergfname = 'energetics_input.dat'
integer, parameter :: ienergread = 103

! initial state specification input file name and unit number
character(30), parameter :: cstatefname = 'state_input.dat'
integer, parameter :: istateread = 104

! Contains the cluster occupancy matrix
character(30), parameter :: clusoccfname = 'clusterocc.bin'
integer, parameter :: clusteroccwrite = 105

! Contains energies
character(30), parameter :: Efname = 'E.bin'
integer, parameter :: Ewrite = 106

! Contains energies
character(30), parameter :: Histfname = 'Hist.bin'
integer, parameter :: Histwrite = 107
! restart file
character(30), parameter :: crestartfname = 'restart.inf'
integer, parameter :: irestart = 201

! general output file
character(30), parameter :: cgenoutfname = 'general_output.txt'
integer, parameter :: iwrite = 202

! lattice specification output file
character(30), parameter :: clatoutfname = 'lattice_output.txt'
integer, parameter :: ilatwrite = 203

! lattice state snapshots (history) output file
character(30), parameter :: chistoryfname = 'history_output.txt'
integer, parameter :: ihistory = 204

! process statistics output file
character(30), parameter :: cprocstatfname = 'procstat_output.txt'
integer, parameter :: iprocstat = 205

! species numbers (transiently) output file
character(30), parameter :: cspecnumfname = 'specnum_output.txt'
integer, parameter :: ispecnum = 206

! output file for properties we want to do SA on
character(30), parameter :: Propfname = 'Prop_output.bin'
integer, parameter :: Propfnum = 208

!Edit: Taylor Robie
!02/22/2016
! Running counter of reaction propensity * dt
character(30), parameter :: PropCountfname = 'PropCounter_output.bin'
integer, parameter :: PropCountfnum = 209

! process debug data output file
character(30), parameter :: cprocdbgfname = 'process_debug.txt'
integer, parameter :: iprocdbg = 301

! global energetics debug data output file
character(30), parameter :: cglbenerdbgfname = 'globalenerg_debug.txt'
integer, parameter :: iglbenergdbg = 302

! Newton's method debug data output file
character(30), parameter :: cnewtndbgfname = 'newton_debug.txt'
integer, parameter :: inewtndbg = 303

! process propensities debug data output file
character(30), parameter :: cprocpropensdbgfname = 'propens_debug.txt'
integer, parameter :: iprocpropensdbg = 304

! process propensities initial catalogue debug data output file
character(30), parameter :: cprocinicatpropensdbgfname = 'propens_inicatalogue.txt'
integer, parameter :: iprocinicatpropensdbg = 305

! cluster contributions debug data output file
character(30), parameter :: cclusterdbgfname = 'cluster_debug.txt'
integer, parameter :: iclusterdbg = 306

! cluster contributions initial catalogue debug data output file
character(30), parameter :: cclusterinidbgfname = 'cluster_inicatalogue.txt'
integer, parameter :: iclusterinidbg = 307

! process propensities debug data state-input generated file
character(30), parameter :: cprocpropensstatedbgfname = 'state_input_inikmc_debug.txt'
integer, parameter :: iprocpropensstatedbg = 308

end module constants_module
