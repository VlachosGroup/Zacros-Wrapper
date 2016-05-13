! This Software is property of University College London. The Licensee shall reproduce a copyright notice on every
! copy of the Software (including partial copies) and on any accompanying manuals and documentation in the form 
! "Copyright (c) University College London, All rights reserved". Trademark and other proprietary notices must also 
! be reproduced but the Licensee has no other right to use the name, arms, trademark, logo or other designation
! of University College London.

module heap_functions_module

use constants_module

! This is an implementation of a binary heap.
! Heap structure setup:
!    heap_elements(1:size)    the elements to be stored in the binary tree
!    heap_labels(1:size)        the labels of the trees elements. These may 
!                            not necessarily be consecutive numbers from 
!                            1 to size, they may be any integers
!    array_indexes(:)        mapping that takes as input an element label
!                            and returns the heap index.

implicit none

integer heapcapacity0

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!subroutine construct_heap(values,n_process, &
!                event_times_heap,event_times_labels,event_times_indexes)
!
!implicit none
!
!include 'z_arraybounds.f'
!
!integer size, itemp, i, j, n_process
!real(8) event_times_heap(heapcapacity0)
!integer event_times_labels(heapcapacity0)
!integer event_times_indexes(heapcapacity0)
!
!Real(8) values(N_cod0)
!
!do i = 1,n_process
!    event_times_heap(i) = values(i)
!    event_times_labels(i) = i
!    event_times_indexes(i) = i
!enddo
!
!call heap_fullsort(event_times_heap,event_times_labels, &
!        event_times_indexes,n_process)
!
!return
!end

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_fullsort(heap_elements,heap_labels,array_indexes,size)

implicit none

integer size, itemp, i, j
real(8) dtemp
real(8) heap_elements(:)
integer heap_labels(:)
integer array_indexes(:)

do i = 1,size
    do j = i+1,size
        if ( .not.(has_priority(heap_elements(i),heap_elements(j))) ) then

            dtemp = heap_elements(i)
            heap_elements(i) = heap_elements(j)
            heap_elements(j) = dtemp

            itemp = heap_labels(i)
            heap_labels(i) = heap_labels(j)
            heap_labels(j) = itemp

            array_indexes(heap_labels(i)) = i
            array_indexes(heap_labels(j)) = j

        endif
    enddo
enddo

return
end subroutine heap_fullsort

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_relabel(heap_elements,heap_labels,array_indexes,size, &
                        cur_element_label,new_element_label)

implicit none

integer new_element_label, cur_element_label, size
real(8) heap_elements(:)
integer heap_labels(:)
integer array_indexes(:)

! This subroutine is intended for filling holes in array_indexes after 
! removal of an element from the heap. It will re-label heap item
! array_indexes(cur_element_label) to new_element_label
! The subroutine will return an error if new_element_label is in use

if (array_indexes(new_element_label) /= 0) then

    write(iwrite,*) 'ERROR: re-labeling not allowed: new index already in use'
    return

else

    heap_labels(array_indexes(cur_element_label)) = new_element_label
    array_indexes(new_element_label) = array_indexes(cur_element_label)
    array_indexes(cur_element_label) = 0
    
endif

end subroutine heap_relabel

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_update_element(heap_elements,heap_labels,array_indexes,size, &
                       updated_element_index,updated_element_value)

implicit none

integer index, newindex, size
integer updated_element_index
real(8) updated_element_value
real(8) heap_elements(heapcapacity0)
integer heap_labels(:)
integer array_indexes(:)
logical boolv1, boolv2, boolv3, boolv4

! This subroutine updates item with label updated_element_index
! and resorts the heap

index = array_indexes(updated_element_index)

if ( (index <= size) .and. (index >= 1) ) then

    ! boolv1 evaluates to .true. iff
    ! (we are updating an internal element) .AND. 
    !    (the updated element has priority over 
    !     the parent of the updated element)
    boolv1 = (elem_parent(index) >= 1)
    if (boolv1) then
        boolv1 = has_priority(updated_element_value,heap_elements(elem_parent(index)))
    endif

    if (boolv1) then ! we need to perform upheap
    
        do while ( elem_parent(index) >= 1 )

            if ( has_priority(updated_element_value,heap_elements(elem_parent(index))) ) then

                heap_elements(index) = heap_elements(elem_parent(index))
                heap_labels(index) = heap_labels(elem_parent(index))

                array_indexes(heap_labels(elem_parent(index))) = index

                index = elem_parent(index)

            else

                exit

            endif

        enddo
    
        heap_elements(index) = updated_element_value
        heap_labels(index) = updated_element_index

        array_indexes(updated_element_index) = index

    else  ! we need to perform downheap

        do
        
            ! boolv2 evaluates to .true. iff
            ! (there exists a left child) .AND. 
            !    (the left child of the updated element has priority 
            !    over the updated element of the heap)
            boolv2 = (elem_left_child(index) <= size)
            if (boolv2) then
                boolv2 = has_priority(heap_elements(elem_left_child(index)),updated_element_value)
            endif

            ! boolv3 evaluates to .true. iff
            ! (there exists a right child) .AND. 
            !    (the right child of the updated element has priority 
            !    over the updated element of the heap)
            boolv3 = (elem_right_child(index) <= size)
            if (boolv3) then
                boolv3 = has_priority(heap_elements(elem_right_child(index)),updated_element_value)
            endif

            if (.not.(boolv2 .or. boolv3)) then
                exit
            endif
    
            ! boolv4 evaluates to .true. iff
            ! (there does not exist a right child) .OR. 
            !    (the left child of the removed element has priority 
            !    over the right child of the removed element)
            boolv4 = (elem_right_child(index) > size)
            if (.not.(boolv4)) then
                boolv4 = has_priority(heap_elements(elem_left_child(index)),heap_elements(elem_right_child(index)))
            endif
                            
            if (boolv4) then
                newindex = elem_left_child(index)
            else
                newindex = elem_right_child(index)
            endif
            
            heap_elements(index) = heap_elements(newindex)
            heap_labels(index) = heap_labels(newindex)

            array_indexes(heap_labels(newindex)) = index

            index = newIndex

        enddo
        
        heap_elements(index) = updated_element_value
        heap_labels(index) = updated_element_index

        array_indexes(updated_element_index) = index
    
    endif
    
else 
    write(iwrite,*) 'ERROR: element not in heap'
    stop
endif

return
end subroutine heap_update_element

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_remove_element(heap_elements,heap_labels,array_indexes,size, &
                       removed_element_label)

implicit none

integer index, newindex, removed_element_label, size
real(8) heap_elements(:)
integer heap_labels(:)
integer array_indexes(:)
logical boolv1, boolv2, boolv3, boolv4

! This subroutine removes item with label removed_element_label
! from the heap

index = array_indexes(removed_element_label)

if ( (index <= size) .and. (index >= 1) ) then

    ! boolv1 evaluates to .true. iff
    ! (we are removing an internal element) .AND. 
    !    (the last element of the heap has priority over 
    !     the parent of the removed element)
    boolv1 = (elem_parent(index) >= 1)
    if (boolv1) then
        boolv1 = has_priority(heap_elements(size),heap_elements(elem_parent(index)))
    endif

    if (boolv1) then ! we need to perform upheap
    
        do while ( elem_parent(index) >= 1 )

            if ( has_priority(heap_elements(size),heap_elements(elem_parent(index))) ) then

                heap_elements(index) = heap_elements(elem_parent(index))
                heap_labels(index) = heap_labels(elem_parent(index))

                array_indexes(heap_labels(elem_parent(index))) = index

                index = elem_parent(index)

            else

                exit

            endif

        enddo
    
        heap_elements(index) = heap_elements(size)
        heap_labels(index) = heap_labels(size)

        array_indexes(heap_labels(size)) = index

    else  ! we need to perform downheap

        do
        
            ! boolv2 evaluates to .true. iff
            ! (there exists a left child) .AND. 
            !    (the left child of the removed element has priority 
            !    over the last element of the heap)
            boolv2 = (elem_left_child(index) <= size)
            if (boolv2) then
                boolv2 = has_priority(heap_elements(elem_left_child(index)),heap_elements(size))
            endif

            ! boolv3 evaluates to .true. iff
            ! (there exists a right child) .AND. 
            !    (the right child of the removed element has priority 
            !    over the last element of the heap)
            boolv3 = (elem_right_child(index) <= size)
            if (boolv3) then
                boolv3 = has_priority(heap_elements(elem_right_child(index)),heap_elements(size))
            endif

            if (.not.(boolv2 .or. boolv3)) then
                exit
            endif
    
            ! boolv4 evaluates to .true. iff
            ! (there does not exist a right child) .OR. 
            !    (the left child of the removed element has priority 
            !    over the right child of the removed element)
            boolv4 = (elem_right_child(index) > size)
            if (.not.(boolv4)) then
                boolv4 = has_priority(heap_elements(elem_left_child(index)),heap_elements(elem_right_child(index)))
            endif
                            
            if (boolv4) then
                newindex = elem_left_child(index)
            else
                newindex = elem_right_child(index)
            endif
            
            heap_elements(index) = heap_elements(newindex)
            heap_labels(index) = heap_labels(newindex)

            array_indexes(heap_labels(newindex)) = index

            index = newIndex

        enddo
        
        if (index < size) then

            heap_elements(index) = heap_elements(size)
            heap_labels(index) = heap_labels(size)

            array_indexes(heap_labels(size)) = index

        endif
    
    endif
    
    array_indexes(removed_element_label) = 0
    heap_elements(size) = 0
    heap_labels(size) = 0

    size = size - 1

else 
    write(iwrite,*) 'ERROR: element not in heap'
    stop
endif

return
end subroutine heap_remove_element

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_insert_element(heap_elements,heap_labels, &
                               array_indexes,size,inserted_element,inserted_index)

implicit none

integer size, index, inserted_index
real(8) inserted_element
real(8) heap_elements(:)
integer heap_labels(:)
integer array_indexes(:)

index = size + 1

if (index > heapcapacity0) then
    write(iwrite,*) 'ERROR: heap maxed out'
    stop
endif

do while (index > 1)
    
    if (has_priority(inserted_element,heap_elements(elem_parent(index)))) then

        heap_elements(index) = heap_elements(elem_parent(index))
        heap_labels(index) = heap_labels(elem_parent(index))

        array_indexes(heap_labels(elem_parent(index))) = index

        index = elem_parent(index)

    else

        exit

    endif

enddo

heap_elements(index) = inserted_element
heap_labels(index) = inserted_index
array_indexes(heap_labels(index)) = index

size = size + 1

end subroutine heap_insert_element

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_check_status(heap_elements,heap_labels, &
                             array_indexes,size)
implicit none

integer size, index
real(8) heap_elements(:)
integer heap_labels(:)
integer array_indexes(:)

do index = 1,size
    
    if (elem_left_child(index) <= size) then

        if ( has_priority(heap_elements(elem_left_child(index)),heap_elements(index)) ) then
            write(iwrite,*) 'Problem in heap node',index,': incorrect sorting - left child has priority'
            stop
        endif

    endif

    if (elem_right_child(index) <= size) then

        if ( has_priority(heap_elements(elem_right_child(index)),heap_elements(index)) ) then
            write(iwrite,*) 'Problem in heap node',index,': incorrect sorting - right child has priority'
            stop
        endif

    endif

    if (array_indexes(heap_labels(index)) /= index) then
        write(iwrite,*) 'Problem in array index ',index,' - heap node ', &
                   array_indexes(index),':incorrect pointer'
        stop
    endif
    
enddo

end subroutine heap_check_status

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function elem_parent(index)

implicit none

integer index

elem_parent = int(floor(dble(index)/2.D0))

return
end function elem_parent

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function elem_left_child(index)

implicit none

integer index

elem_left_child = 2*index

return
end function elem_left_child

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function elem_right_child(index)

implicit none

integer index

elem_right_child = 2*index + 1

return
end function elem_right_child

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

logical function has_priority(element1,element2)

implicit none

real(8) element1, element2

if (element1 < element2) then
    has_priority = .true.
else
    has_priority = .false.
endif

return
end function has_priority

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module heap_functions_module
