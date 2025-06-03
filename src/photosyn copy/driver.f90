program main

    use photosyn, only : can_photosyn
    
    real(8),allocatable :: gpp(:)
    real(8):: gppx(20000)
    
    call can_photosyn(gpp)
    
    !gppx(1:365*24) = gpp
!    gppx(:) = gpp
    
!    print*, gppx(3000:3100)
!    print*, gppx(19990:20000)
    print*,'success!'

end program main