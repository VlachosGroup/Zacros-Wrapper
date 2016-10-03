! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module output_handle_module

use lattice_setup_module
use simulation_setup_module
use mechanism_setup_module
use energetics_setup_module, only: nclusters, clusterOcc, clustergraphmultipl

use kmc_simulation_handle_module
use energetics_handle_module, only: globalenergy
use lattice_handle_module
use sampling_handle_module

implicit none

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine write_history_header()

implicit none

integer i

if (.not.output_snapshots) return

write(ihistory,'(' // int2str(ngasspecs+1) // '(a,2x))') &
            'Gas_Species:        ', (trim(gasspecsnames(i)),i=1,ngasspecs)

write(ihistory,'(' // int2str(nsurfspecs+1) // '(a,2x))') &
            'Surface_Species:    ', (trim(surfspecsnames(i)),i=1,nsurfspecs)

if (allocated(cellsiteneighsf)) then ! If periodic lattice has been specified save the simulation box
    write(ihistory,'(a)') 'Simulation_Box:         '
    do i = 1,2
        write(ihistory,'(2x,2(ES32.16E3,2x))') ncellsrep(i)*vcell(1,i), ncellsrep(i)*vcell(2,i)
    enddo
endif

write(ihistory,'(' // int2str(nsitetypes+1) // '(a,2x))') &
            'Site_Types:         ', (trim(sitetypenames(i)),i=1,nsitetypes)

end subroutine write_history_header

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine save_snaphshot()

implicit none

integer i, j

if (.not.output_snapshots) return

if (snap_on_event) then

    if (mod(curstep-1_8,dkeventsnap) == 0_8) then

        snapshnum = snapshnum + 1_8

        ! *** Write sample to history file
        write(ihistory,'(a,I20,1x,I20,1x,ES30.16,1x,ES30.16,1x,ES30.16)') 'configuration ', snapshnum, curstep-1_8, prevtime, temp+tramp*curtime, globalenergy
        do i = 1,nsites
            write(ihistory,'(4(I10,1x))') i,(latticestate(i,j), j = 1,3)
        enddo
        write(ihistory,'('//int2str(ngasspecs)//'(I20,1x))') (gasspecsnums(j), j = 1,ngasspecs)
    
    endif
        
else
    
    if (snap_on_logtime) then
        
        do while (snaptime <= min(curtime,maxtime) + 2*tiny(dtsnap))
        
            snapshnum = snapshnum + 1_8

            ! *** Write sample to history file
            write(ihistory,'(a,I20,1x,I20,1x,ES30.16,1x,ES30.16,1x,ES30.16)') 'configuration ', snapshnum, curstep-1_8, snaptime, temp+tramp*snaptime, globalenergy
            do i = 1,nsites
                write(ihistory,'(4(I10,1x))') i,(latticestate(i,j), j = 1,3)
            enddo
            write(ihistory,'('//int2str(ngasspecs)//'(I20,1x))') (gasspecsnums(j), j = 1,ngasspecs)
            snaptime = snaptime*dtsnap
            
        enddo
            
    else
    
        do while (snaptime <= min(curtime,maxtime) + 2*tiny(dtsnap))
        
            snapshnum = snapshnum + 1_8

            ! *** Write sample to history file
            write(ihistory,'(a,I20,1x,I20,1x,ES30.16,1x,ES30.16,1x,ES30.16)') 'configuration ', snapshnum, curstep-1_8, snaptime, temp+tramp*snaptime, globalenergy
            do i = 1,nsites
                !write(ihistory,'(4(I10,1x))') i,(latticestate(i,j), j = 1,3)
				!
				!write(Histwrite) i
				!write(Histwrite) latticestate(i,1)
				!write(Histwrite) latticestate(i,2)
				!write(Histwrite) latticestate(i,3)
				
            enddo
            write(ihistory,'('//int2str(ngasspecs)//'(I20,1x))') (gasspecsnums(j), j = 1,ngasspecs)
            snaptime = snaptime + dtsnap
            
        enddo
    
    endif
    
endif    

return

end subroutine save_snaphshot

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine write_procstat_header()

implicit none

integer i

if (.not.output_procstat) return

! If processes statistics saving requested write the names of the elementary 
! steps in the header of the corresponding file
write(iprocstat,'(' // int2str(nelemsteps+1) // '(a30,1x))') 'Overall', &
                                      (trim(elemstepnames(i)),i=1,nelemsteps)

end subroutine write_procstat_header

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine save_procstatistics()

implicit none

integer i

if (.not.output_procstat) return

if (procstat_on_event) then

    if (mod(curstep-1_8,dkeventprocstat) == 0_8) then
        
        procstatnum = procstatnum + 1_8

        ! *** Write process statistics info
        write(iprocstat,'(a,I20,1x,I20,1x,ES30.16)') 'configuration ', procstatnum, curstep-1_8, prevtime
    
        write(iprocstat,'(' // int2str(nelemsteps+1) // '(ES30.16,1x))') (elemstep_avgtime(i), i = 0,nelemsteps)
        write(iprocstat,'(' // int2str(nelemsteps+1) // '(I20,1x))') (elemstep_noccur(i), i = 0,nelemsteps)
        
    endif
    
else
    
    if (procstat_on_logtime) then
    
        do while (procstattime <= min(curtime,maxtime) + 2*tiny(dtprocstat))
        
            procstatnum = procstatnum + 1_8

            ! *** Write process statistics info
            write(iprocstat,'(a,I20,1x,I20,1x,ES30.16)') 'configuration ', procstatnum, curstep-1_8, procstattime
            
            write(iprocstat,'(' // int2str(nelemsteps+1) // '(ES30.16,1x))') (elemstep_avgtime(i), i = 0,nelemsteps)
            write(iprocstat,'(' // int2str(nelemsteps+1) // '(I20,1x))') (elemstep_noccur(i), i = 0,nelemsteps)
            procstattime = procstattime*dtprocstat

        enddo

    else
    
        do while (procstattime <= min(curtime,maxtime) + 2*tiny(dtprocstat))
        
            procstatnum = procstatnum + 1_8

            ! *** Write process statistics info
            write(iprocstat,'(a,I20,1x,I20,1x,ES30.16)') 'configuration ', procstatnum, curstep-1_8, procstattime
            
            write(iprocstat,'(' // int2str(nelemsteps+1) // '(ES30.16,1x))') (elemstep_avgtime(i), i = 0,nelemsteps)
            write(iprocstat,'(' // int2str(nelemsteps+1) // '(I30,1x))') (elemstep_noccur(i), i = 0,nelemsteps)
            procstattime = procstattime + dtprocstat

        enddo
    
    endif
    
endif    

return

end subroutine save_procstatistics

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine write_specnums_header()

implicit none

integer i

if (.not.output_specnum) return

! If species number reporting requested write the names of the surface and
! gas species in the header of the corresponding file
write(ispecnum,'(a20,1x,a20,1x,a30,1x,a30,1x,a30,' // trim(int2str(nsurfspecs+ngasspecs)) // '(a20,1x))') &
          'Entry','Nevents','Time','Temperature','Energy',(trim(surfspecsnames(i)),i=1,nsurfspecs), &
          (trim(gasspecsnames(i)),i=1,ngasspecs)

end subroutine write_specnums_header

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine save_specnums(mproc)

use state_setup_module, only: nadsorb

implicit none

integer i, mproc, i2

if (.not.output_specnum) return

if (specnum_on_event) then

    if (specnum_on_eleventoccur) then
        
        if (proctypesites(mproc,0) == ieleventoccur) then
        
            specnumnum = specnumnum + 1_8

            ! *** Write species numbers info
                write(ispecnum,'(I20,1x,I20,1x,ES30.16,1x,ES30.16,1x,ES30.16,' // trim(int2str(nsurfspecs+ngasspecs)) // '(I20,1x))') &
                  specnumnum, curstep-1_8, prevtime,temp+tramp*prevtime,globalenergy, &
                  (sum(adsorbspecposi(1:nadsorb,0),mask = adsorbspecposi(1:nadsorb,0) == i)/i,i=1,nsurfspecs), &
                  (gasspecsnums(i),i=1,ngasspecs)
				  
			! Extra binary output
			!write(Specfnum) (sum(adsorbspecposi(1:nadsorb,0),mask = adsorbspecposi(1:nadsorb,0) == i)/i,i=1,nsurfspecs)	! Record species numbers in binary file
			!write(clusteroccwrite) (clusterOcc(i) / clustergraphmultipl(i), i=1, nclusters)
			!write(Ewrite) specnumtime	
			!write(Ewrite) globalenergy
			!write(Propfnum) (propvec(i), i=1, nSAparams)
			write(PropCountfnum) (propCountvec(i) + propvec(i) * (specnumtime - prevtime), i=1, nSAparams)
			!write(procstatfnum) (elemstep_noccur(i), i = 1, nelemsteps)																		! Write process statistics info			
			write(SAfnum) (elemstep_noccur(i) - ( propCountvec(i) + propvec(i) * (specnumtime - prevtime) ), i=1, nSAparams)				! Record W for sensitivity analysis 
			
        endif
        
    else
        
        if (mod(curstep-1_8,dkeventspecnum) == 0_8) then
        
            specnumnum = specnumnum + 1_8

            ! *** Write species numbers info
                write(ispecnum,'(I20,1x,I20,1x,ES30.16,1x,ES30.16,1x,ES30.16,' // trim(int2str(nsurfspecs+ngasspecs)) // '(I20,1x))') &
                  specnumnum, curstep-1_8, prevtime,temp+tramp*prevtime,globalenergy, &
                  (sum(adsorbspecposi(1:nadsorb,0),mask = adsorbspecposi(1:nadsorb,0) == i)/i,i=1,nsurfspecs), &
                  (gasspecsnums(i),i=1,ngasspecs)        
        
			! Extra binary output
			!write(Specfnum) (sum(adsorbspecposi(1:nadsorb,0),mask = adsorbspecposi(1:nadsorb,0) == i)/i,i=1,nsurfspecs)	! Record species numbers in binary file
			!write(clusteroccwrite) (clusterOcc(i) / clustergraphmultipl(i), i=1, nclusters)
			!write(Ewrite) specnumtime	
			!write(Ewrite) globalenergy
			!write(Propfnum) (propvec(i), i=1, nSAparams)
			write(PropCountfnum) (propCountvec(i) + propvec(i) * (specnumtime - prevtime), i=1, nSAparams)
			!write(procstatfnum) (elemstep_noccur(i), i = 1, nelemsteps)																		! Write process statistics info			
			write(SAfnum) (elemstep_noccur(i) - ( propCountvec(i) + propvec(i) * (specnumtime - prevtime) ), i=1, nSAparams)				! Record W for sensitivity analysis 
		
        endif
        
    endif
    
else
    
    if (specnum_on_logtime) then
    
        do while (specnumtime <= min(curtime,maxtime) + 2*tiny(dtspecnum))
        
            specnumnum = specnumnum + 1_8

            ! *** Write species numbers info
            write(ispecnum,'(I20,1x,I20,1x,ES30.16,1x,ES30.16,1x,ES30.16,' // trim(int2str(nsurfspecs+ngasspecs)) // '(I20,1x))') &
                  specnumnum, curstep-1_8, specnumtime,temp+tramp*specnumtime,globalenergy, &
                  (sum(adsorbspecposi(1:nadsorb,0),mask = adsorbspecposi(1:nadsorb,0) == i)/i,i=1,nsurfspecs), &
                  (gasspecsnums(i),i=1,ngasspecs)     
				  
			! Extra binary output
			!write(Specfnum) (sum(adsorbspecposi(1:nadsorb,0),mask = adsorbspecposi(1:nadsorb,0) == i)/i,i=1,nsurfspecs)	! Record species numbers in binary file
			!write(clusteroccwrite) (clusterOcc(i) / clustergraphmultipl(i), i=1, nclusters)
			!write(Ewrite) specnumtime
			!write(Ewrite) globalenergy
			!write(Propfnum) (propvec(i), i=1, nSAparams)
			write(PropCountfnum) (propCountvec(i) + propvec(i) * (specnumtime - prevtime), i=1, nSAparams)
			!write(procstatfnum) (elemstep_noccur(i), i = 1, nelemsteps)																		! Write process statistics info			
			write(SAfnum) (elemstep_noccur(i) - ( propCountvec(i) + propvec(i) * (specnumtime - prevtime) ), i=1, nSAparams)				! Record W for sensitivity analysis 
			
            specnumtime = specnumtime * dtspecnum

        enddo

    else

        do while (specnumtime <= min(curtime,maxtime) + 2*tiny(dtspecnum))
        
            specnumnum = specnumnum + 1_8

            ! *** Write species numbers info
            write(ispecnum,'(I20,1x,I20,1x,ES30.16,1x,ES30.16,1x,ES30.16,' // trim(int2str(nsurfspecs+ngasspecs)) // '(I20,1x))') &
                  specnumnum, curstep-1_8, specnumtime,temp+tramp*specnumtime,globalenergy, &
                  (sum(adsorbspecposi(1:nadsorb,0),mask = adsorbspecposi(1:nadsorb,0) == i)/i,i=1,nsurfspecs), &
                  (gasspecsnums(i),i=1,ngasspecs)
				  
			! Extra binary output
			!write(Specfnum) (sum(adsorbspecposi(1:nadsorb,0),mask = adsorbspecposi(1:nadsorb,0) == i)/i,i=1,nsurfspecs)	! Record species numbers in binary file
			!write(clusteroccwrite) (clusterOcc(i) / clustergraphmultipl(i), i=1, nclusters)
			!write(Ewrite) specnumtime	
			!write(Ewrite) globalenergy
			!write(Propfnum) (propvec(i), i=1, nSAparams)
			!write(PropCountfnum) (propCountvec(i) + propvec(i) * (specnumtime - prevtime), i=1, nSAparams) 	! Include truncation term
			write(PropCountfnum) (propCountvec(i), i=1, nSAparams)
			!write(procstatfnum) (elemstep_noccur(i), i = 1, nelemsteps)																		! Write process statistics info			
			write(SAfnum) (elemstep_noccur(i) - propCountvec(i), i=1, nSAparams)
			!write(SAfnum) (elemstep_noccur(i) - ( propCountvec(i) + propvec(i) * (specnumtime - prevtime) ), i=1, nSAparams)				! Record W for sensitivity analysis, Include truncation term
				  
            specnumtime = specnumtime + dtspecnum

        enddo

    endif

endif    

do i = 1,nelemsteps
	if (dtPrior > 0.d0) then
		propCountvec(i) = propCountvec(i) + propvec(i) * dtPrior
	endif
end do

end subroutine save_specnums

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine dump_propensity_tables()

use state_setup_module, only: nadsorb

implicit none

integer i, iproc, j, si, sj, k, m, tempint, cntr
integer, allocatable :: indxk(:)

real(8) tempvar !, activenrg, preexpfac, eventpropensity0
real(8), allocatable :: propensarray(:)

logical flagsitedif, takesprecedensejoveri

! Caution: this subroutine must be ran as the last procedure before terminating
! i.e. *after* writing the restart file

do cntr = 1,2
    
    if (cntr == 1) then
    
        allocate(indxk(nprocesses))
        allocate(propensarray(nprocesses))
        do iproc = 1,nprocesses
            indxk(iproc) = iproc    
            propensarray(iproc) = procpropenst0(iproc)
        enddo
    
    elseif (cntr == 2) then

        deallocate(event_times_heap)
        deallocate(event_times_labels)
        deallocate(event_times_indexes)
        deallocate(proctypesites)
        deallocate(procpropenst0)
        deallocate(procdeltaenrg)
        deallocate(nadsorprocparticip)
        deallocate(adsorprocparticip)
        deallocate(indxk)
        deallocate(propensarray)
        nprocesses = 0
        curtime = 0.d0
        debug_report_processes = .false.
        call catalogue_all_processes()
        allocate(indxk(nprocesses))
        allocate(propensarray(nprocesses))
        do iproc = 1,nprocesses
            indxk(iproc) = iproc    
            propensarray(iproc) = procpropenst0(iproc)           
        enddo
    
    endif
    
    ! In the following loops we are just sorting the propensities based on the:
    !  (i)  magnitude
    !  (ii) the sites at which the corresponding process takes place
    ! The propensities with higher magnitude take precedence. If two propensities correspond to the
    ! same elementary step, the one happening at the site with the higher number takes precedence (and
    ! of course for multi-site processes, each site is checked). It may happen that two propensities 
    ! are equal but correspond to processes with different site types, in this case the one with the 
    ! more sites participating takes precedence.

    do i = 1,nprocesses
        do j = i+1,nprocesses
        
            takesprecedensejoveri = .false.
            
            if (propensarray(i) < propensarray(j)) then
                
                takesprecedensejoveri = .true.
                
            elseif ( (dabs(propensarray(i)) < 1.d-135 .and. dabs(propensarray(j)) < 1.d-135) .or. &
                 (dabs(propensarray(i)/propensarray(j)-1.d0) < 5.d-16) ) then
                
    !            if (dabs(propensarray(i)) < 1.d-15) then
    !            continue
    !            endif
                
                if (elemstepnsites(proctypesites(indxk(i),0)) < elemstepnsites(proctypesites(indxk(j),0))) then
                    
                    takesprecedensejoveri = .true.
                    
                elseif (elemstepnsites(proctypesites(indxk(i),0)) == elemstepnsites(proctypesites(indxk(j),0))) then
                    
                    flagsitedif = .false.
                    
                    do k = 1,elemstepnsites(proctypesites(indxk(i),0))
                        si = proctypesites(indxk(i),k)
                        sj = proctypesites(indxk(j),k)
                        if (si < sj) then
                            takesprecedensejoveri = .true.
                            flagsitedif = .true.
                            exit
                        elseif (si > sj) then
                            flagsitedif = .true.
                            exit                
                        endif
                    enddo
                    
                    if (.not.flagsitedif .and. proctypesites(indxk(i),0) < proctypesites(indxk(j),0)) then
                        takesprecedensejoveri = .true.
                    endif
                    
                endif
            
            endif
            
            if (takesprecedensejoveri) then
                tempvar = propensarray(i)
                propensarray(i) = propensarray(j)
                propensarray(j) = tempvar
                tempint = indxk(i)
                indxk(i) = indxk(j)
                indxk(j) = tempint
            endif
            
        enddo
    enddo
    
    if (cntr == 1) then
    
        open(unit=iprocpropensdbg,file=trim(cprocpropensdbgfname),status='unknown')

        do i = 1,nprocesses
            k = indxk(i)
            write(iprocpropensdbg,'(ES32.16E3,' // int2str(elemstepnsites(proctypesites(k,0))+1) // '(1x,i10))') &
                propensarray(i), proctypesites(k,0), &
                (proctypesites(k,m), m = 1,elemstepnsites(proctypesites(k,0)))
        enddo

        close(iprocpropensdbg)
        
    elseif (cntr == 2) then
    
        open(unit=iprocinicatpropensdbg,file=trim(cprocinicatpropensdbgfname),status='unknown')

        do i = 1,nprocesses
            k = indxk(i)
            write(iprocinicatpropensdbg,'(ES32.16E3,' // int2str(elemstepnsites(proctypesites(k,0))+1) // '(1x,i10))') &
                propensarray(i), proctypesites(k,0), &
                (proctypesites(k,m), m = 1,elemstepnsites(proctypesites(k,0)))
        enddo

        close(iprocinicatpropensdbg)
            
    endif

enddo

open(unit=iprocpropensstatedbg,file=trim(cprocpropensstatedbgfname),status='unknown')

write(iprocpropensstatedbg,'(a)') 'initial_state'

do i = 1,nadsorb
    if (adsorbspecposi(i,0) > 0) then
        write(iprocpropensstatedbg,'(a,' // int2str(surfspecsdent(adsorbspecposi(i,0))) // '(1x,i10))') &
            '   seed_on_sites ' // trim(surfspecsnames(adsorbspecposi(i,0))), &
            (adsorbspecposi(i,k), k = 1,surfspecsdent(adsorbspecposi(i,0)))
    endif
enddo

write(iprocpropensstatedbg,'(a)') 'end_initial_state'

close(iprocpropensstatedbg)

end subroutine dump_propensity_tables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine dump_cluster_contribution_tables()

use state_setup_module, only: nadsorb
use energetics_setup_module, only: clustergraphmultipl, clusterenrg,           &
                                   clusternsites

implicit none

integer i, iclust, j, jcluster, si, sj, k, m, tempint, cntr
integer, allocatable :: indxk(:)

real(8) tempvar
real(8), allocatable :: energarray(:)

logical flagsitedif, takesprecedensejoveri

! Caution: this subroutine must be ran as the last procedure before terminating
! i.e. *after* writing the restart file

do cntr = 1,2
    
    if (cntr == 1) then
    
        allocate(indxk(nglobclust))
        allocate(energarray(nglobclust))
        do iclust = 1,nglobclust
            indxk(iclust) = iclust   
            jcluster = globclustypesites(iclust,0)
            energarray(iclust) = clusterenrg(jcluster)/clustergraphmultipl(jcluster)
        enddo
    
    elseif (cntr == 2) then

        deallocate(globclustypesites)
        deallocate(nadsorglobclusparticip)
        deallocate(adsorglobclusparticip)
        deallocate(indxk)
        deallocate(energarray)
        nglobclust = 0
        debug_report_globenerg = .false.
        call initialize_energetics()
        allocate(indxk(nglobclust))
        allocate(energarray(nglobclust))
        do iclust = 1,nglobclust
            indxk(iclust) = iclust    
            jcluster = globclustypesites(iclust,0)
            energarray(iclust) = clusterenrg(jcluster)/clustergraphmultipl(jcluster)
        enddo
    
    endif
    
    ! In the following loops we are just sorting the energy contributions based on the:
    !  (i)  magnitude
    !  (ii) the sites of the corresponding cluster
    ! The clusters with higher energies take precedence. If two global clusters correspond to the
    ! same pattern, the one involving a site with the higher number takes precedence (and
    ! of course for multi-site clusters, each site is checked). It may happen that two global cluster energies 
    ! are equal but correspond to clusters with different site types, in this case the one with the 
    ! more sites participating takes precedence.

    do i = 1,nglobclust
        do j = i+1,nglobclust
        
            takesprecedensejoveri = .false.
            
            if (energarray(i) < energarray(j)) then
                
                takesprecedensejoveri = .true.
                
            elseif ( (dabs(energarray(i)) < 1.d-35 .and. dabs(energarray(j)) < 1.d-35) .or. &
                 (dabs(energarray(i)/energarray(j)-1.d0) < 5.d-16) ) then
                
    !            if (dabs(energarray(i)) < 1.d-15) then
    !            continue
    !            endif
                
                if (clusternsites(globclustypesites(indxk(i),0)) < clusternsites(globclustypesites(indxk(j),0))) then
                    
                    takesprecedensejoveri = .true.
                    
                elseif (clusternsites(globclustypesites(indxk(i),0)) == clusternsites(globclustypesites(indxk(j),0))) then
                    
                    flagsitedif = .false.
                    
                    do k = 1,clusternsites(globclustypesites(indxk(i),0))
                        si = globclustypesites(indxk(i),k)
                        sj = globclustypesites(indxk(j),k)
                        if (si < sj) then
                            takesprecedensejoveri = .true.
                            flagsitedif = .true.
                            exit
                        elseif (si > sj) then
                            flagsitedif = .true.
                            exit                
                        endif
                    enddo
                    
                    if (.not.flagsitedif .and. globclustypesites(indxk(i),0) < globclustypesites(indxk(j),0)) then
                        takesprecedensejoveri = .true.
                    endif
                    
                endif
            
            endif
            
            if (takesprecedensejoveri) then
                tempvar = energarray(i)
                energarray(i) = energarray(j)
                energarray(j) = tempvar
                tempint = indxk(i)
                indxk(i) = indxk(j)
                indxk(j) = tempint
            endif
            
        enddo
    enddo
    
    if (cntr == 1) then
    
        open(unit=iclusterdbg,file=trim(cclusterdbgfname),status='unknown')

        do i = 1,nglobclust
            k = indxk(i)
            write(iclusterdbg,'(ES32.16E3,' // int2str(clusternsites(globclustypesites(k,0))+1) // '(1x,i10))') &
                energarray(i), globclustypesites(k,0), &
                (globclustypesites(k,m), m = 1,clusternsites(globclustypesites(k,0)))
        enddo

        close(iclusterdbg)
        
    elseif (cntr == 2) then
    
        open(unit=iclusterinidbg,file=trim(cclusterinidbgfname),status='unknown')

        do i = 1,nglobclust
            k = indxk(i)
            write(iclusterinidbg,'(ES32.16E3,' // int2str(clusternsites(globclustypesites(k,0))+1) // '(1x,i10))') &
                energarray(i), globclustypesites(k,0), &
                (globclustypesites(k,m), m = 1,clusternsites(globclustypesites(k,0)))
        enddo

        close(iclusterinidbg)
            
    endif

enddo

open(unit=iprocpropensstatedbg,file=trim(cprocpropensstatedbgfname),status='unknown')

write(iprocpropensstatedbg,'(a)') 'initial_state'

do i = 1,nadsorb
    if (adsorbspecposi(i,0) > 0) then
        write(iprocpropensstatedbg,'(a,' // int2str(surfspecsdent(adsorbspecposi(i,0))) // '(1x,i10))') &
            '   seed_on_sites ' // trim(surfspecsnames(adsorbspecposi(i,0))), &
            (adsorbspecposi(i,k), k = 1,surfspecsdent(adsorbspecposi(i,0)))
    endif
enddo

write(iprocpropensstatedbg,'(a)') 'end_initial_state'

close(iprocpropensstatedbg)

end subroutine dump_cluster_contribution_tables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module output_handle_module
