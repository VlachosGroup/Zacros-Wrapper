! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module rates_handle_module

use constants_module
use mechanism_setup_module, only: calculate_preex_factor, tempfun, nSAparams

implicit none

! Declares a type which is singularly used to hold work arrays.
! The issue is that some functions are called over and over. We do not want to
! have to allocate and reallocate those arrays.  Currently, only holds arrays
! previously 
type work_array_type
  integer, allocatable :: adsexclude(:)
  integer, allocatable :: clustercontrib(:)
  integer, allocatable :: reactadsnums(:)
  integer, allocatable :: neighlist(:)
  integer, allocatable :: indxs(:)
  integer, allocatable :: validpatrns(:,:)
  integer, allocatable :: adsorbspecposi(:, :)
end type 

! A private work array.
! When/if we move to multi-processing, we may need more than one. 
type (work_array_type), private, target, allocatable :: work_prv(:)
private :: work_array_type
private :: get_work_array
contains
  
  subroutine prealloc_rhm(numthreads)
    use energetics_setup_module, only: clusterlevels, nclusters,               &
                                       clusternmxsites
    use mechanism_setup_module, only: elemstepnmxsites, elemstepreactnts
    use lattice_setup_module, only: maxcoord, nsites
    use simulation_setup_module, only: maxdent
    ! Does memory pre-allocation. 
    ! Should be called after setup by only one thread
    ! It should be called only once.
    integer, intent(in) :: numthreads
    integer :: i

    if( .not. allocated(work_prv)) then
      allocate(work_prv(numthreads))
      do i = 1, numthreads
        allocate(work_prv(i)%adsexclude(0:elemstepnmxsites))
        allocate(work_prv(i)%clustercontrib(0:10*nclusters*clusternmxsites))
        allocate(work_prv(i)%reactadsnums(maxval(elemstepreactnts(:,0))))
        allocate(work_prv(i)%validpatrns(100*maxcoord,clusternmxsites))
        allocate(work_prv(i)%neighlist(100*maxcoord))
        allocate(work_prv(i)%indxs(0:maxval(clusterlevels)))
        allocate(work_prv(i)%adsorbspecposi(nsites + elemstepnmxsites, 0:maxdent))
      enddo
    endif

  end subroutine prealloc_rhm

  subroutine cleanup_rhm
    ! Deletes work arrays
    integer :: i

    if(allocated(work_prv)) then
      do i = 1, size(work_prv) 
        deallocate(work_prv(i)%adsexclude)
        deallocate(work_prv(i)%clustercontrib)
        deallocate(work_prv(i)%reactadsnums)
        deallocate(work_prv(i)%validpatrns)
        deallocate(work_prv(i)%neighlist)
        deallocate(work_prv(i)%indxs)
        deallocate(work_prv(i)%adsorbspecposi)
      enddo
      deallocate(work_prv)
    endif

  end subroutine cleanup_rhm 

  function get_work_array()
      ! This function returns a pointer to an allocated work array.
      ! In multiprocessing, we may need to more than one work array. We will only
      ! need modify this function and add a locking mechanism.
      ! Should be called in parallel section only
      
!$    use omp_lib, only: omp_in_parallel, omp_get_thread_num
      
      type(work_array_type), pointer :: get_work_array
      
!$    if(omp_in_parallel()) then
!$        get_work_array => work_prv(omp_get_thread_num()+1)
!$    else
          get_work_array => work_prv(1)
!$    endif
  end function get_work_array
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine calculate_elemstep_rate( jelemstep, latsitesmap, activenrg,         &
                                    activenrg0, deltaenrg, deltaenrg0,         &
                                    deltalattenerg, eventpropensity0, der,          &
                                    eventtime, randnum)
use random_deviates_module, only: expon_dev_from_uniform
use simulation_setup_module, only: inewtonstats, rnewtonstats, tpdsim, curtime,&
                                   etol1, etol2, inewtndbg, ngp, kmaxnewton,   &
                                   temp, tramp, gp, gw
use mechanism_setup_module, only: debug_newtons_method, elemstepnames, tempfun, timeoftemp0
use constants_module, only: enrgconv, kboltz, timeconv
use error_module, only: moreinfo, error, warning, warncounters
use parser_module, only: dbl2str, int2str

implicit none

integer, intent(in) :: jelemstep
integer, intent(in) :: latsitesmap(:)
real(8), intent(in) :: randnum

real(8), intent(out) :: activenrg, activenrg0
real(8), intent(out) :: deltaenrg, deltaenrg0
real(8), intent(out) :: deltalattenerg
real(8), intent(out) :: eventpropensity0
real(8), dimension(nSAparams), intent(out) :: der							! the derivative of propensity with respect to some parameter
real(8), intent(out) :: eventtime

integer i, k, i2

real(8) r1, lninvphi, tgp, integral_eventpropensity, err1, err2, timetemp0
real(8) incremeventtime, incremeventtimeprev, preexpfac, eventpropensity

call calculate_activ_energy( jelemstep, latsitesmap, activenrg,                &
                             activenrg0, deltaenrg, deltaenrg0,                &
                             deltalattenerg )                             
                                                                          
call calculate_preex_factor(jelemstep,temp,preexpfac)         
eventpropensity0 = (preexpfac*timeconv)                                        &
                   * dexp(-(activenrg/enrgconv)/kboltz/temp)            

if (.not.tpdsim) then

	! Compute the derivative for sensitivity analysis, key place where analytical derivative information is included
	do i2 = 1,nSAparams
      IF (i2 == jelemstep) THEN
        der(i2) = eventpropensity0				! da/dlnk = a
      ELSE
        der(i2) = 0								! 0 if the rate constant and reaction type do not match
      ENDIF
    end do

    eventpropensity = eventpropensity0
    
    eventtime = curtime + expon_dev_from_uniform(eventpropensity, randnum)

else ! temperature ramp has been defined

    call calculate_preex_factor( jelemstep,                                    &
                                 tempfun(temp,tramp,curtime,1),                &
                                 preexpfac )

    if (dabs(preexpfac) < 1.D-150) then
       eventtime = huge(1.d0)
       return
    endif

    r1 = randnum
    lninvphi = -dlog(r1)
    
    ! If the simulation time has reached infinity (no more events can happen)
    ! just set eventtime to infinity This is done to allow the code to terminate
    ! properly
    if (curtime >= huge(1.d0)) then
        eventtime = huge(1.d0)
        return
    endif
    
    ! Time at which the absolute temperature reaches zero (useful if T is dropping)
    timetemp0 = timeoftemp0(temp,tramp,1)
    
    ! Initial guess for the time needed for the next event occurrence 
    ! Obtain initial guess assuming no change in the event propensity
    eventpropensity = (preexpfac*timeconv)                                     &
                      * dexp( -(activenrg/enrgconv) / kboltz                   &
                               / tempfun(temp,tramp,curtime,1) )
    incremeventtimeprev = 1.d0/eventpropensity*lninvphi
    
    !! Refine initial guess
    !if (curtime < 0.9*maxtime) then
    !    incremeventtimeprev = (maxtime-curtime)/2.d0
    !else
    !    incremeventtimeprev = 1.d-15
    !endif

    if (debug_newtons_method) then
        !$OMP CRITICAL
        write(inewtndbg,'(a)') 'curtime=' // trim(dbl2str(curtime))            &
                         // 'incremeventtimeprev = '                           &
                         // trim(dbl2str(incremeventtimeprev))
        write(inewtndbg,'(6A20)') 'iter', 'incremeventtimeprev',               &
                                  'incremeventtime', 'err1', 'err2'
        !$OMP END CRITICAL
    endif
    
    k = 0

    do while (.true.)

        ! *** Calculate the propensity function and its integral at the event
        ! time
        call calculate_preex_factor( jelemstep,                                &                                     
                           tempfun(temp,tramp,curtime+incremeventtimeprev,1),  & 
                                     preexpfac)
        eventpropensity = (preexpfac*timeconv)                                 &
                          * dexp( -(activenrg/enrgconv)/kboltz                 &
                        / tempfun(temp,tramp,curtime+incremeventtimeprev,1) )

        integral_eventpropensity = 0.d0
        
        ! MSTAM 10-AUG-2013: the block below assumes that the propensity is 0 at 
        ! Temperature = 0 K and thus, we only need to calculate the integral up to
        ! the time that T = 0 K, if this time is finite and less than 
        ! curtime+incremeventtimeprev, of course. These apply to simmulated annealing 
        ! calculations where we quench the system. Note that this is OK for thermal 
        ! rates. If we want to simulate quantum tunnelling we will need to add
        ! one more line that computes the rest of the integral. This should be easy,
        ! as the propensity for times > t-where-T=0K will remain constant.
        do i = 1,ngp
            ! Gauss points are from 0 to 1: make the transformation from
            ! absolute t. The transformation is xi = tau/t => t*dxi = dtau
            tgp = curtime + gp(i)*min(incremeventtimeprev,max(timetemp0-curtime,0.d0))

            call calculate_preex_factor( jelemstep,                            &
                                         tempfun(temp,tramp,tgp,1), preexpfac )
          
            integral_eventpropensity = integral_eventpropensity                &
                                       + gw(i)*incremeventtimeprev             &
                                         * ( (preexpfac*timeconv)              &
                                         * dexp( -(activenrg/enrgconv)         &
                                                 / kboltz                      &
                                                 /tempfun(temp,tramp,tgp,1) ) )
            
            ! continue
        enddo
    
        ! *** Perform next Newton iteration
        k = k + 1
    
        ! The equation is:
        !  integral_eventpropensity - ln(1/r2) = 0
        ! and we are solving:
        !  ln(integral_eventpropensity) - ln(lninvphi) = 0
        ! The logarithmization improves convergence since integral_eventpropensity is a steep 
        ! function of time.
    
        incremeventtime = incremeventtimeprev                                  &
                          - (dlog(integral_eventpropensity/lninvphi))          &
                            /(eventpropensity/integral_eventpropensity)

        if (incremeventtime < 0.d0) then
            incremeventtime = 1.d-15
        endif
        
        ! *** Check norms and proceed accordingly
        err1 = dabs(incremeventtime - incremeventtimeprev)                     &
               / dabs(incremeventtime)
        err2 = dabs(integral_eventpropensity/lninvphi-1.d0)    
        
        if (debug_newtons_method) then
           !$OMP CRITICAL
            write(inewtndbg,'(I20,5ES20.7E3)') k, incremeventtimeprev,         &
                                               incremeventtime, err1, err2, r1
           !$OMP END CRITICAL
        endif
        
        if (integral_eventpropensity > tiny(1.d0) .and.                        &
            eventpropensity < tiny(1.d0) .and.                                 &
            dabs(integral_eventpropensity/lninvphi-1.d0) > etol2) then
            
            if (debug_newtons_method) then
                !$OMP CRITICAL
                write(inewtndbg,'(a)') 'It seems there is no solution. ' //    &
                'The time will be set to a huge number.'
                write(inewtndbg,'(a)') '~~~~~~~~~~~~~~~~~~~~~~'//              &
                '~~~~~~~~~~~~~~~~~~~~~~~~~'
                !$OMP END CRITICAL
            endif
            
            if (warncounters(2) == maxrepwarnings) then
                !$OMP CRITICAL
                moreinfo = 'Elem. step ' // trim(int2str(jelemstep))           &
                           // ', ' // trim(elemstepnames(jelemstep))           &
                           // '; RelErrDx = '&
                           // trim(dbl2str(err1)) // '. RelErrRHS = '          &
                           // trim(dbl2str(err2)) // char(13) // char(10)      &
                           // 'curtime = ' // trim(dbl2str(curtime))           &
                           // '; incremeventtime = '                           &
                           // trim(dbl2str(incremeventtime))                   &
                           // '; incremeventtimeprev = '                       &
                           // trim(dbl2str(incremeventtimeprev))               &
                           // '; eventpropensity = '                           &
                           // trim(dbl2str(eventpropensity))                   &
                           // char(13) // char(10)                             &
                           // 'integral_eventpropensity = '                    &
                           // trim(dbl2str(integral_eventpropensity))          &
                           // '; lninvphi = ' // trim(dbl2str(lninvphi))
                call warning(10007)
                !$OMP END CRITICAL
            endif
            
            eventtime = huge(1.d0)
            return
            
        endif        
        
        ! this evaluates to .true. if eventtime == NaN
        if (.not.(incremeventtime == incremeventtime)) then 
            !$OMP CRITICAL
            moreinfo = 'Elem. step ' // trim(int2str(jelemstep))               &
                       // ', ' // trim(elemstepnames(jelemstep))               &
                       // char(13) // char(10) // 'curtime = '                 &
                       // trim(dbl2str(curtime)) // '; incremeventtime = '     &
                       // trim(dbl2str(incremeventtime))                       &
                       // '; incremeventtimeprev = '                           &
                       // trim(dbl2str(incremeventtimeprev)) // char(13)       &
                       // char(10) // 'eventpropensity = '                     &
                       // trim(dbl2str(eventpropensity))                       &
                       // '; integral_eventpropensity = '                      &
                       // trim(dbl2str(integral_eventpropensity))              &
                       // '; lninvphi = ' // trim(dbl2str(lninvphi))
            call error(5011)
            !$OMP END CRITICAL
        endif
    
        if (((err1 < etol1) .and. (err2 < etol2)) .or. (k == kmaxnewton)) exit
    
        incremeventtimeprev = incremeventtime

    enddo

    eventtime = curtime + incremeventtime

    if (k == kmaxnewton) then
        
        if (warncounters(1) == maxrepwarnings) then
            !$OMP CRITICAL
            moreinfo = 'Elem. step ' // trim(int2str(jelemstep))                &
                       // ', ' // trim(elemstepnames(jelemstep))                &
                       // '; RelErrDx = '                                       &
                       // trim(dbl2str(err1)) // '. RelErrRHS = '               &
                       // trim(dbl2str(err2)) // char(13) // char(10)           &
                       // 'curtime = ' // trim(dbl2str(curtime))                &
                       // '; incremeventtime = '                                &
                       // trim(dbl2str(incremeventtime))                        &
                       // '; incremeventtimeprev = '                            &
                       // trim(dbl2str(incremeventtimeprev))                    &
                       // '; eventpropensity = '                                &
                       // trim(dbl2str(eventpropensity))                        &
                       // char(13) // char(10)                                  &
                       // 'integral_eventpropensity = '                         &
                       // trim(dbl2str(integral_eventpropensity))               &
                       // '; lninvphi = ' // trim(dbl2str(lninvphi))            
            call warning(10006)
            !$OMP END CRITICAL
        endif
        
        inewtonstats(1) = inewtonstats(1) + 1
        
    endif

    if (debug_newtons_method) then
        !$OMP CRITICAL
        write(inewtndbg,'(a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        !$OMP END CRITICAL
    endif

    inewtonstats(0) = inewtonstats(0) + 1
    rnewtonstats(0) = (rnewtonstats(0)*(inewtonstats(0)-1) + k)                &
                      / dble(inewtonstats(0))
    rnewtonstats(1) = max(rnewtonstats(1),err1)
    rnewtonstats(2) = max(rnewtonstats(2),err2)

endif

return

end subroutine calculate_elemstep_rate

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine calculate_activ_energy( jelemstep, latsitesmap, activenrg,          &
                                   activenrg0, deltaenrg, deltaenrg0,          &
                                   deltalattenerg)
use mechanism_setup_module, only: acteng, omega, reverselemstep

implicit none

integer, intent(in) :: jelemstep
integer, intent(in) :: latsitesmap(:)
real(8), intent(out) :: activenrg, activenrg0, deltaenrg, deltaenrg0,          &
                        deltalattenerg

call calculate_derxn( jelemstep, latsitesmap, deltaenrg, deltaenrg0,           &
                      deltalattenerg )

if (reverselemstep(jelemstep) == 0) then  
! we are looking at an irreversible reaction

    activenrg0 = acteng(jelemstep)
    activenrg = activenrg0 + omega(jelemstep)*(deltaenrg-deltaenrg0)

elseif (reverselemstep(jelemstep) > jelemstep) then 
! we are looking at the forward step of a reversible reaction
    
    activenrg0 = acteng(jelemstep)
    activenrg = activenrg0 + omega(jelemstep)*(deltaenrg-deltaenrg0)

    activenrg0 = max(0.d0,activenrg0,deltaenrg0)
    activenrg = max(0.d0,activenrg,deltaenrg)
        
else  
! we are looking at the reverse step of a reversible reaction

    ! The deltaenrg that we have calculated is for the reverse elementary step
    ! that we are interested in ...
    
    ! that's why the opposite sign here ...
    activenrg0 = acteng(reverselemstep(jelemstep)) + deltaenrg0 
    activenrg = activenrg0 + (1.d0-omega(reverselemstep(jelemstep)))           &
                             * (deltaenrg-deltaenrg0) ! ... and here.

    activenrg0 = max(0.d0,activenrg0,deltaenrg0)
    activenrg = max(0.d0,activenrg,deltaenrg)
    
endif

return

end subroutine calculate_activ_energy

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine calculate_derxn(jelemstep,latsitesmap,derxn,derxn0,deltalattenerg)

use energetics_setup_module, only: clusterangles, clusternsites,               &
                                   specclusterparticip, clusteradsorbpos,      &
                                   clusternanglessequenc, clusterneigh,        &
                                   clusternmxangles, clusternmxcoord,          &
                                   clusternomirrorimgs, clusterabslorientat,   &
                                   clusterorientationedge,                     &
                                   clusterorientationangles, clusternmxsites,  &
                                   clusterlatticestate, clusterspecs 
use state_setup_module, only: nadsorb_global => nadsorb
use energetics_handle_module, only: adsorglobclusparticip, globclustypesites,  &
                                    nadsorglobclusparticip
use energetics_setup_module, only: clusterenrg, clustergraphmultipl,           &
                                   clusterstype, clusterlevels,                &
                                   nspecclusterparticip
use lattice_handle_module, only: adsorbspecposi, get_latticestate,             &
                                 check_lattice_state, find_pattern_match
use lattice_setup_module, only: siteneighb1, sitetype
use mechanism_setup_module, only: debug_check_lattice, elemstepadsorbposin,    &
                                  elemstepreactnts, elemstepproducts,          &
                                  elemstepadsorbposfn, elemstepgases
use simulation_setup_module, only: surfspecsdent, gasenergies
use graph_functions_module, only: find_dlevel_neighbors

implicit none
integer, intent(in) :: jelemstep
integer, intent(in) :: latsitesmap(:)
real(8), intent(out) :: derxn, derxn0, deltalattenerg 

integer i, ik, ispec, dmax
integer j, jsitep, jsiteelstep, jsitelattic, jglobclust, jadsorb, jcluster,    &
        k, kadsorb, kspec, ksites, m, mmolec, nneigh, nstart, cntrpatrn,       &
        nadsexclude, nloops, ipatrnsite, ilattcsite, threadnum, nadsorb

real(8) gasenerg, enrginistate, enrginistate0, enrgfinstate, enrgfinstate0

type(work_array_type), pointer :: work
integer, pointer, dimension(:, :) :: latticestate
work => get_work_array()
latticestate => get_latticestate()

derxn = 0.d0
derxn0 = 0.d0

work%adsexclude(0) = elemstepreactnts(jelemstep,0)
work%clustercontrib(0) = 0

do i = 1, elemstepreactnts(jelemstep,0) ! for each of the reactant species

    jsiteelstep = elemstepadsorbposin(jelemstep,i,1)
    jsitelattic = latsitesmap(jsiteelstep)
    jadsorb = latticestate(jsitelattic,1) ! get the adsorbate/entity number 
    work%adsexclude(i) = jadsorb

enddo

do i = 1, elemstepreactnts(jelemstep,0) ! for each of the reactant species

    jadsorb = work%adsexclude(i)

    ! for each of the global clusters this entity participates        
    gloclustloop: do j = 1,nadsorglobclusparticip(jadsorb) 

        jglobclust = adsorglobclusparticip(jadsorb,j,1)
        mmolec = adsorglobclusparticip(jadsorb,j,2)
        jcluster = globclustypesites(jglobclust,0) ! get the cluster type
        
        ! Check if that cluster has already been encountered. For this
        ! we note that if any of the excluded adsorbates participates in
        ! the cluster as an entity with order less than mmolec, then the
        ! cluster has already been counted.
        do m = 1,mmolec-1
            ! the site that the molecule occupies on the lattice
            jsitep = globclustypesites( jglobclust,                            &
                                        clusteradsorbpos(jcluster,m,1) )
            kadsorb = latticestate(jsitep,1)
        
            if( any(work%adsexclude(1:work%adsexclude(0)) == kadsorb) )        &
               cycle gloclustloop
        enddo

        work%clustercontrib(0) = work%clustercontrib(0) + 1
        work%clustercontrib(work%clustercontrib(0)) = jglobclust
    enddo gloclustloop
    
enddo

! Calculate the initial state energy contributions from the reactants
! configuration as well as at the zero coverage limit
enrginistate = 0.d0
enrginistate0 = 0.d0

clustloopeng0: do j = 1,work%clustercontrib(0)
    jglobclust = work%clustercontrib(j)
    jcluster = globclustypesites(jglobclust,0) ! get the cluster type

    enrginistate = enrginistate                                                &
                   + clusterenrg(jcluster)/clustergraphmultipl(jcluster)

    do k = 1,clusternsites(jcluster)
      if( .not. any( latticestate(globclustypesites(jglobclust,k),1)           &
                     .eq. work%adsexclude(1:work%adsexclude(0))) )             &
          cycle clustloopeng0
    enddo
    
    enrginistate0 = enrginistate0                                              &
                    + clusterenrg(jcluster)/clustergraphmultipl(jcluster)
    
enddo clustloopeng0

! nadorb used to be modified then unmodified. That's really bad mojo.
! Now use local version.
nadsorb = nadsorb_global
! Check lattice structure (for debugging)
if (debug_check_lattice) then
   call check_lattice_state(1,nadsorb,adsorbspecposi,latticestate)
endif

! Make a temporary change in the state of the lattice 

! Save entity numbers for the reactants at the initial state
do i = 1,elemstepreactnts(jelemstep,0)
    work%reactadsnums(i)                                                       &
      = latticestate( latsitesmap(elemstepadsorbposin(jelemstep,i,1)), 1)
enddo

work%adsexclude(0) = 0

! We add the product adsorbates temporarily, numbered as
! nadsorb+1:nadsorb+elemstepproducts(jelemstep,0)
do i = 1, elemstepproducts(jelemstep,0)
    
    nadsorb = nadsorb + 1

    ispec = elemstepproducts(jelemstep,i)

    work%adsexclude(0) = work%adsexclude(0) + 1
    work%adsexclude(work%adsexclude(0)) = nadsorb
    
    work%adsorbspecposi(nadsorb, 0) = ispec
    
    do j = 1,surfspecsdent(ispec)
    
        jsiteelstep = elemstepadsorbposfn(jelemstep,i,j)
        jsitelattic = latsitesmap(jsiteelstep)
        latticestate(jsitelattic,1) = nadsorb
        latticestate(jsitelattic,2) = ispec
        latticestate(jsitelattic,3) = j
        
        work%adsorbspecposi(nadsorb, j) = jsitelattic
            
    enddo
    
enddo

! Check lattice for (partial) consistency
if (debug_check_lattice) then
   call check_lattice_state(nadsorb-elemstepproducts(jelemstep,0)+1,nadsorb, work%adsorbspecposi,latticestate)
endif

threadnum = 1
nadsexclude = work%adsexclude(0) 
enrgfinstate = 0.d0
enrgfinstate0 = 0.d0
do i = 1, nadsexclude
    
    ik = work%adsexclude(i)
    kspec = work%adsorbspecposi(ik, 0)

    clusterloop: do j = 1, nspecclusterparticip(kspec)

        mmolec = specclusterparticip(kspec, j, 2)
        if(mmolec < 1 .or. mmolec > clusternmxsites) cycle clusterloop

        jcluster = specclusterparticip(kspec, j, 1)
        ! *** INITIAL CHECKS ***
               
            
        ! We will immediately check to see whether molecule i occupies on
        ! the lattice the site types that the pattern jcluster requires        
        do k = 1, surfspecsdent(kspec)
            ipatrnsite = clusteradsorbpos(jcluster, mmolec, k) 
            ilattcsite = work%adsorbspecposi(ik, k)
            if ( clusterstype(jcluster, ipatrnsite) > 0                    &
                 .and. ( clusterstype(jcluster, ipatrnsite)                &
                         /= sitetype(ilattcsite) ) )                       &
                cycle clusterloop
        enddo

        ! *** CANDIDATE SITES ***
        
        ! Make a list of candidate sites of the lattice. These sites could
        ! be participating in a process pertaining to elementary step
        ! jelemstep.
        ! dmax = max level of elementary step i starting from
        ! sites occupied by molecule i
        nstart = surfspecsdent(kspec)
        dmax = clusterlevels(jcluster, mmolec) 
        
        call find_dlevel_neighbors( siteneighb1, dmax, nstart,             &
                                    work%adsorbspecposi(ik, 1:nstart),      &
                                    nneigh, work%indxs, work%neighlist)
                                        
        ! *** FINDING VALID PATTERNS ***
        call find_pattern_match( jcluster, mmolec, clusternmxcoord,        &
                                 clusternmxangles, clusternsites,          &
                                 clusterspecs, clusterstype, clusterneigh, &
                                 clusterlatticestate, clusteradsorbpos,    &
                                 clusternanglessequenc, clusterangles,     &    
                                 clusternomirrorimgs,                      &
                                 clusterabslorientat, &
                                   clusterorientationedge, &
                                     clusterorientationangles, &
                                 ! the number of candidate sites and the
                                 ! vector containing them
                                 nneigh, work%neighlist, work%adsexclude,  &
                                 cntrpatrn, work%validpatrns)
        
        ! *** ADDING IDENTIFIED ENERGETIC CONTRIBUTIONS ***
        enrgfinstate = cntrpatrn * clusterenrg(jcluster)                   &
                       / clustergraphmultipl(jcluster)                     &
                       + enrgfinstate
        
        enrgfinstate0 = count_finstates(cntrpatrn, clusternsites(jcluster),&
                                        work%validpatrns,                  &
                                        latticestate(:, 1),                &
                                        work%adsexclude(1:nadsexclude) )   &
                        * clusterenrg(jcluster)                            &
                        / clustergraphmultipl(jcluster)                    &
                        + enrgfinstate0
        
    enddo clusterloop
enddo 
                
derxn = enrgfinstate - enrginistate
derxn0 = enrgfinstate0 - enrginistate0

! Save the lattice energy change due to the reaction
deltalattenerg = derxn

! Add the contributions from the gas species to the reaction energy at the given
! coverage and at the zero coverage limit
do i = 1,elemstepgases(jelemstep,0)
    gasenerg = gasenergies(elemstepgases(jelemstep,2*i-1))                     &
               * elemstepgases(jelemstep,2*i)
    derxn0 = derxn0 + gasenerg
    derxn = derxn + gasenerg
enddo

! Now revert the change
do i = 1,elemstepreactnts(jelemstep,0)
    
    ispec = elemstepreactnts(jelemstep,i)

    do j = 1,surfspecsdent(ispec)
    
        jsiteelstep = elemstepadsorbposin(jelemstep,i,j)
        jsitelattic = latsitesmap(jsiteelstep)
        
        latticestate(jsitelattic,1) = work%reactadsnums(i)
        latticestate(jsitelattic,2) = elemstepreactnts(jelemstep,i)
        latticestate(jsitelattic,3) = j

    enddo
    
enddo

! Check lattice structure (for debugging)
if (debug_check_lattice) then
    call check_lattice_state(1,nadsorb_global,adsorbspecposi,latticestate)
endif

return

contains
! Helper subroutines needed only here (for now)

  integer function count_finstates(cntrpatrn, nsites, validpatrns,             &
                                   latticestate, adsexclude)
    ! Counts number of valid states for enrgfinstate0
    !
    ! As far as I can see, this whole algorithm performs set intersections
    ! between adsexclude and non-trivial segments of the latticestate array.
    integer, intent(in) :: cntrpatrn
    integer, intent(in) :: nsites
    integer, intent(in) :: validpatrns(:, :)
    integer, intent(in) :: latticestate(:)
    integer, intent(in) :: adsexclude(:)

    integer :: w, k, ksites


    count_finstates = 0
    patternloop: do k = 1, cntrpatrn ! For every pattern...
        
        do ksites = 1, nsites
          w = latticestate(validpatrns(k, ksites))
          if(.not. any(adsexclude .eq. w)) cycle patternloop
        enddo
        
        count_finstates = count_finstates + 1
    
    enddo patternloop
  end function count_finstates

end subroutine calculate_derxn



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


end module rates_handle_module
