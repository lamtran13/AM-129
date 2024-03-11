! File: fput.f90
! Author: Lam Tran
! Purpose: This is going to be your main driver routine which calls the
! appropriate routines from the remaining files.

program fput
    use utility
    use setup_module
    use leapfrog_module
    use output_module

    implicit none

     call setup_init() 

      allocate(x(nMasses))
      allocate(v(nMasses))
    
     call set_ics()

     print *, 'initial conditions:'
     print *, (x(i), i = 1,nMasses)

     call leapfrog()
  
     print *, 'leapfrog result:'
     print *, (x(i), i = 1,nMasses)
    
     call output_write()

      deallocate(x)
      deallocate(v)

end program fput