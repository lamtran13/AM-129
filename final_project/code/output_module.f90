! File: output_module.f90
! This writes the solution out to disc as needed.
! Author: Lam Tran

module output_module
    use setup_module
    use utility
    implicit none

 contains

    subroutine output_write()
        implicit none
        character(len=35) :: ofile

    !! File name for ascii output
        ofile = 'sol_'//trim(runName)//'.dat'

    !! File open
        open(unit=20, file=ofile, status='unknown')

    !! Write into a file:
        do i=1, nSteps
            write(20,*) i, (result(j,i), j = 1,nMasses)
        end do

    !! Output format specifiers
    !920 format(1x, i5, f16.8, 1x, f16.8, 1x, f16.12)

    !! File close
        close(20)

  end subroutine output_write

end module output_module
