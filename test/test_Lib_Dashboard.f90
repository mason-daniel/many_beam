!   simple test code for dashboard

    program test_Lib_Dashboard
!---^^^^^^^^^^^^^^^^^^^^^^^^^^
        use Lib_Dashboard
        implicit none

        call setState(COLUMNAR)
        call reportDashboard()
        print *,""

        call unsetState(COLUMNAR)
        call setState(LOSSY)
        call reportDashboard()
        print *,""
        

        print *,""
        print *,"PASS"
        print *,""
    end program test_Lib_Dashboard