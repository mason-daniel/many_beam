
!   A dashboard file for the many_beam code
!

    integer,public,parameter        ::      PROGRAM_STATE = int32
    integer,private,parameter       ::      FLAG_NAME_LEN = 32
    integer,private,parameter       ::      NUMBER_FLAGS = 8

    character(len=FLAG_NAME_LEN),dimension(NUMBER_FLAGS),private,parameter  ::      FLAG_NAME =     &
                (/  "columnar approx                 ",                                                     &
                    "thermal diffuse scatter         ",                                                     &
                    "two beam approx                 ",                                                     &
                    "gauge transform                 ",                                                     &
                    "output phase field              ",                                                     &
                    "output .png                     ",                                                     &
                    "output diffract pattern         ",                                                     &
                    "Fibonacci spiral precession     "   /)

    integer(kind=PROGRAM_STATE),public,parameter    ::      COLUMNAR        = int( b'00000001',kind=PROGRAM_STATE )
    integer(kind=PROGRAM_STATE),public,parameter    ::      LOSSY           = int( b'00000010',kind=PROGRAM_STATE )
    integer(kind=PROGRAM_STATE),public,parameter    ::      TWOBEAM         = int( b'00000100',kind=PROGRAM_STATE )
    integer(kind=PROGRAM_STATE),public,parameter    ::      GAUGETRANSFORM  = int( b'00001000',kind=PROGRAM_STATE )
    integer(kind=PROGRAM_STATE),public,parameter    ::      OPXYZ           = int( b'00010000',kind=PROGRAM_STATE )
    integer(kind=PROGRAM_STATE),public,parameter    ::      OPPNG           = int( b'00100000',kind=PROGRAM_STATE )
    integer(kind=PROGRAM_STATE),public,parameter    ::      OPDIFF          = int( b'01000000',kind=PROGRAM_STATE )
    integer(kind=PROGRAM_STATE),public,parameter    ::      FIBONACCI_PREC  = int( b'10000000',kind=PROGRAM_STATE )
                