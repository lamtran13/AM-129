! leapfrog_module.f90
! Author: Lam Tran
! This holds and updates the solution depending
! on the physical parameters and initial conditions

module leapfrog_module
    use utility
    use setup_module
    implicit none
    
    public :: leapfrog

contains

    subroutine leapfrog()
        implicit none
        
        allocate(x0(nMasses),x1(nMasses))
        allocate(result(nMasses,nSteps)) 

        x0 = x
        x1 = x

        do n = 1,nSteps
         do i = 2,nMasses-1
           x(i) = 2*x1(i)-x0(i)+K*dt**2*(x1(i+1)-2*x1(i)+x1(i-1))*(1+alpha*(x1(i+1)-x1(i-1)))
         end do

          x0 = x1
          x1 = x

         result(:,n) = x
        end do 
    
        deallocate(x0,x1)
    end subroutine leapfrog

end module leapfrog_module