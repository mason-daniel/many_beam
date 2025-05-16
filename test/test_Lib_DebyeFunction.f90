
!*      A simple test program to show functioning of DebyeFunction
!*      gfortran -ffree-line-length=256 Lib_ColouredTerminal.f90 Lib_DebyeFunction.f90 test_Lib_DebyeFunction.f90 -o test_Lib_DebyeFunction

    program test_Lib_DebyeFunction
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      Daniel Mason
!*      (c) UKAEA Jan 2025
        use Lib_ColouredTerminal
        use Lib_DebyeFunction
        use iso_fortran_env
        implicit none

        ! real(kind=real64)       ::      x
        ! integer                 ::      n

        real(kind=real64),dimension(3)          ::      testResult  
        real(kind=real64),dimension(3)          ::      groundTruth = (/0.77750463411233972d0,0.49308264399060786d0,-0.19552792724292167d0/)
        real(kind=real64)                       ::      tol = 1d-12
        
        !integer                 ::      ii

    !---    test Debye function 
        testResult(1) = DebyeFunction(1,1.0d0)
        print *,"DebyeFunction(1,1.0d0) ",testResult(1),numericalApprox(1,1.0d0)
        testResult(2) = DebyeFunction(2,2.0d0)
        print *,"DebyeFunction(2,2.0d0) ",testResult(2),numericalApprox(2,2.0d0)
          
        testResult(3) = dDebyeFunction(1,1.0d0)        
        print *,"dDebyeFunction(1,1.0d0) ",testResult(3),( DebyeFunction(1,1.0d0+1.0d-6) - DebyeFunction(1,1.0d0-1.0d-6) )/2.0d-6


        if (compareResult(testResult,groundTruth,tol)) then
            print *,colour(GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if


    contains
!---^^^^^^^^

        pure real(kind=real64) function numericalApprox(n,x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                  ::      n
            real(kind=real64),intent(in)        ::      x
            real(kind=real64)           ::      yy
            integer,parameter           ::      NSTEPS = 10000000
             
            integer                     ::      ii
            real(kind=real64)           ::      dy,zz

            dy = x/NSTEPS
            numericalApprox = 0.0d0
            do ii = 0,NSTEPS
                yy = ii*dy
                if (ii==0) then
                    zz = 0.0d0
                    if (n==1) zz = 1.0d0
                else
                    zz = yy**n/(exp(yy)-1)
                end if
                if ( (ii==0).or.(ii==NSTEPS) ) then
                    numericalApprox = numericalApprox + zz 
                else if (mod(ii,2)==1) then
                    numericalApprox = numericalApprox + 4*zz 
                else
                    numericalApprox = numericalApprox + 2*zz 
                end if
            end do
            numericalApprox = numericalApprox * n * x**(-n) * dy / 3
            return
        end function numericalApprox

    end program test_Lib_DebyeFunction