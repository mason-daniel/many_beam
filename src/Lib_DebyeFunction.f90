
!*  A module to compute the Debye function 
!*  
!*      Dn(x) = n/x^n   int_t=0^x     t^n/(e^t - 1) dt
!*
!*  https://en.wikipedia.org/wiki/Debye_function
!*  
!*      in range |x| < 2 pi , n>= 1


    module Lib_DebyeFunction
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*      Daniel Mason 
!*      (c) UKAEA Jan 2025
        use iso_fortran_env
        implicit none
        private

        include "Bernoulli_numbers.h"
        real(kind=real64),private,parameter     ::      PI  = 3.141592653589790d0

        public      ::      DebyeFunction
        public      ::      dDebyeFunction
        public      ::      BernoulliNumber
        public      ::      BernoulliFunction
        public      ::      DebyeWallerFactor_approx


    contains
!---^^^^^^^^

        pure real(kind=real64) function BernoulliFunction(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute Bernooulli function 
    !*      f(x) = x/(exp[x]-1)
            real(kind=real64),intent(in)        ::      x
            real(kind=real64),parameter         ::      MAXX = 690.0d0          !   I want to compute exp(x), which must be < 1d308.
            if (x>MAXX) then 
                BernoulliFunction = 0.0d0
            else if (x<-MAXX) then
                BernoulliFunction = -x
            else
                BernoulliFunction = 1/(exp(x) - 1)
            end if
            return
        end function BernoulliFunction


        pure real(kind=real64) function BernoulliNumber(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute nth Bernoulli number
    !*      note that B_n = 0 for n odd if n>1
            integer,intent(in)              ::      n
            if (n<0) then
                BernoulliNumber = 0.0d0
            else if (n>ubound(Bernoulli_numerator,dim=1)) then
                if (mod(n,2)==1) then
                    BernoulliNumber = 0.0d0
                else
                    BernoulliNumber = 1.0d0
                end if
            else
                BernoulliNumber = Bernoulli_numerator(n) / Bernoulli_denominator(n)
            end if
            return
        end function BernoulliNumber

        pure real(kind=real64) function DebyeFunction(n,x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the Debye function D_n(x)
            integer,intent(in)              ::      n
            real(kind=real64),intent(in)    ::      x
            
            integer             ::      kk
            real(kind=real64)   ::      x2kon2kfact,term


            if (n<1) then
                DebyeFunction = 0
            else if (x>=2*PI) then
                DebyeFunction = n*gamma( n+1.0d0 )*RiemannZetaFunction( n+1 )*x**(-n)
            else if (x<=-2*PI) then
                DebyeFunction = 0.0d0
            else                  
                DebyeFunction = 1 - n*x/(2*n+2)
                x2kon2kfact = 1.0d0
                do kk = 1,ubound(Bernoulli_numerator,dim=1)
                    x2kon2kfact = x2kon2kfact * x*x / ( 2*kk*(2*kk-1) )
                    term = BernoulliNumber(2*kk)*x2kon2kfact / (2*kk+n)
                    if (abs(term) < 1.0d-20*DebyeFunction) exit
                    DebyeFunction = DebyeFunction + n*term
                end do

            !   smoothly transition between analytic approximations at the end of the convergence region
                if (x>6) then
                    term = (2*PI - x)/(2*PI - 6)
                    DebyeFunction = DebyeFunction*term + n*gamma( n+1.0d0 )*RiemannZetaFunction( n+1 )*x**(-n)*(1-term)
                end if
                
            end if

            return
        end function DebyeFunction


        pure real(kind=real64) function dDebyeFunction(n,x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the derivative of the Debye function D_n(x)
            integer,intent(in)              ::      n
            real(kind=real64),intent(in)    ::      x
            
            if (abs(x)>0) then
                dDebyeFunction = n*( BernoulliFunction(x) - DebyeFunction(n,x) )/x
            else
                dDebyeFunction = 0.0d0
            end if

            return
        end function dDebyeFunction


        pure real(kind=real64) function RiemannZetaFunction(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      zeta function for small integer values of n
            integer,intent(in)          ::      n
            real(kind=real64),dimension(2:10),parameter         ::      ZETA = (/       &
                                        1.6449340668d0,                            &
                                        1.2020569032d0,                            &
                                        1.0823232337d0,                            &
                                        1.0369277551d0,                            &
                                        1.0173430620d0,                            &
                                        1.0083492774d0,                            &
                                        1.0040773562d0,                            &
                                        1.0020083928d0,                            &
                                        1.0009945751d0         /)
            
            
            if (n<=1) then
                RiemannZetaFunction = huge(1.0d0)
            else if (n>ubound(ZETA,dim=1)) then
                RiemannZetaFunction = 1.0d0
            else
                RiemannZetaFunction = ZETA(n)
            end if
            return
        end function RiemannZetaFunction
        
        pure real(kind=real64) function DebyeWallerFactor_approx(T,A,ThetaD)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the Willis & Pryor (1975) approximation for the DWF
    !*      given the temperature, Debye-Waller Temp and the atomic mass A  (in mass units, C12 = 12)
    !*      return answer in A^2
            real(kind=real64),intent(in)            ::      T       !   temperature
            real(kind=real64),intent(in)            ::      A       !   atomic mass
            real(kind=real64),intent(in)            ::      ThetaD  !   dbeye temp

            real(kind=real64)       ::      xx,zz
            
            xx = T/ThetaD
            if (xx*2*PI <= 1) then
                !   T is very low. So can't use full debye func approx.                
                zz = gamma( 2.0d0 )*RiemannZetaFunction( 2 )
            else
                zz = xx*DebyeFunction(1,1/xx)
            end if

            DebyeWallerFactor_approx = ( 11492/(A*ThetaD) ) * ( zz + 0.25d0 )
            

            return
        end function DebyeWallerFactor_approx

    end module Lib_DebyeFunction