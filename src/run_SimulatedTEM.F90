
    program run_SimulatedTEM
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*      read in a file, compute the simulated TEM image
!*      and output as png
!*
!*      Daniel Mason
!*      (c) UKAEA October 2024
!*      
!*      version history
!*          v.0.0.1         Oct 24      First working version
!*          v.0.0.2         Feb 26      bright field vs dark field imaging permitted. 
!*          v.0.0.3         Feb 26      Fix to selection rules - thanks to Thomas Jourdan (CEA) for spotting mod vs modulo bug
!*          v.0.0.4         Mar 26      Option to output image from all diffracted beams
!*          v.0.0.5         Apr 26      Option to find centre-of-mass in diffraction pattern
!*          v.0.1.0         May 26      thanks to Aslak Fellman (Helsinki) spotted bug in handling of unknown element type. Change to default op filename for .xyz files. dump .xyz file in tilted position ( serial run only )
!*          v.0.1.1         Oct 26      tweaks to output of diffraction patterns
!*          v.1.0.0         Nov 26      Major version number 

       
        use iso_fortran_env
        use Lib_Lattices
        use Lib_Elements
        use Lib_XYZFiles
        use Lib_Filenames
        use Lib_Png
        use Lib_IntegrateManyBeams
        use Lib_ImagingSpace
        use Lib_AtomSpace
        use Lib_Quaternions
        use Lib_Gvectors
        use Lib_Callipers
        use Lib_CommandLineArguments
        use Lib_ComputePhaseFactor
        use Lib_RotationMatrices
        use Lib_ManyBeam
        use Lib_Dashboard
#ifdef MPI
        use mpi_f08
#endif
        implicit none


        character(len=8),parameter          ::      VERSION = "1.0.0"


        real(kind=real64),parameter         ::      PI = 3.14159265390d0
        integer                             ::      rank = 0,nProcs = 1

        integer,parameter                   ::      CLA_FILE   = 1
        integer,parameter                   ::      CLA_MICRO  = 2
        integer,parameter                   ::      CLA_DIFF   = 3
        integer,parameter                   ::      CLA_CALC   = 4
        integer,parameter                   ::      CLA_OUTPUT = 5

    !---    command line parameter input        
        type(CommandLineArguments)          ::      cla
        
        character(len=256)                  ::      filename = ""           !   input filename
        character(len=256)                  ::      outfile = ""            !   output filename prefix
        character(len=8)                    ::      latticename = UNKNOWN_LATTICE
        real(kind=real64)                   ::      T = 300.0d0             !   temperature (K)
        real(kind=real64)                   ::      V = 200.0d0             !   accelerator voltage (kV)
        integer                             ::      M = 4                   !   image subdivisions per lattice param
        real(kind=real64),dimension(3)      ::      a0_in = 0.0             !   indicative lattice parameter (A)
        integer                             ::      nPrec = 1               !   precession angles
        real(kind=real64)                   ::      precAngle = 0.003d0     !   precession angle (rad)
        logical                             ::      attset = .false.        
        integer,dimension(:),allocatable    ::      hkl                     !   the g-vector for the aperture. allocatable because its either hkl or hjkl 
        logical                             ::      yflip = .true.          !   Right-hand rule for png output - y=0 at bottom.
        logical                             ::      pngBlur = .true.        !   blurring radius sigma/a
        logical                             ::      pngNorm = .true.        !   normalise output intensity
        real(kind=real64)                   ::      n_g = 1                 !   which reflection to make bright
        logical                             ::      df = .true. , bf =.false.       !   find dark field - note bf not used except in command line interface -bright
        logical                             ::      all_beams = .false.     !   construct image from all beams 
        logical                             ::      com_beams = .false.     !   construct image from centre of mass of beams
        real(kind=real64)                   ::      tilt_max = 5.0          !   maximum stage tilt (deg)
        integer                             ::      ntilt = 1000            !   number of tilt angles searched for good diffraction pattern
        real(kind=real64)                   ::      tweak = 0.01d0          !   tweak in tilt angle permitted to optimise diffraction pattern
        real(kind=real64)                   ::      imin = 0.001d0          !   intensity threshold for "important" beam
        integer                             ::      zlayers = 50            !   number of z layers to compute rho & x* grad x
        logical                             ::      dbg = .false.   
        real(kind=real64),dimension(9)      ::      Rmat = (/ 1,0,0,0,1,0,0,0,1 /)
        real(kind=real64),dimension(9)      ::      Umat = (/ 1,0,0,0,1,0,0,0,1 /)
        logical                             ::      specifiedRmat = .false.
        logical                             ::      specifiedUmat = .false.
        integer                             ::      select_grain = LIB_CLA_NODEFAULT_I
        logical                             ::      tiltonly = .false.      !   option to quit after finding goo tilt angle

    !---    physical variables
        type(IntegrateManyBeams)            ::      imb
        integer                             ::      Nx,Ny           !   image size
        real(kind=real64),dimension(:,:),allocatable        ::      img         !   greyscale image of intensities
        real(kind=real64),dimension(:,:,:),allocatable      ::      com_img     !   coloured image of centre of mass
        integer,dimension(4)                ::      dat
        character(len=256)                  ::      dummy
        integer,dimension(:,:),allocatable  ::      hkl_important
        type(AtomSpace)                     ::      as

    !---    dummy variables        
        integer             ::      ierror
        integer             ::      ii,ix,iy,iyp
        real(kind=real64)   ::      dd
        character(len=5)    ::      aa
        logical             ::      ok
        integer,dimension(2,3)  ::  bb
        real(kind=real64),dimension(3,3)    ::  a_super
        real(kind=real64),dimension(3)      ::  rr

    !---    timing
        integer,parameter               ::      TIMER = kind(1_int32)
        integer(kind=TIMER),parameter   ::      T_INIT      = 1_TIMER
        integer(kind=TIMER),parameter   ::      T_TILT      = 2_TIMER
        integer(kind=TIMER),parameter   ::      T_OUTPUT    = 3_TIMER
        integer(kind=TIMER),parameter   ::      T_INTEGRATE = 4_TIMER
        integer(kind=TIMER),parameter   ::      T_TOTAL     = 6_TIMER
        type(Callipers),dimension(T_TOTAL)      ::      tt




!******************************************************************************    
!*    
!*      START
!*    
!*    
!******************************************************************************    
    

        tt(T_TOTAL) = Callipers_ctor() 
        tt(T_INIT) = Callipers_ctor()
#ifdef MPI        
        call MPI_INIT(ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
#endif              



!******************************************************************************    
!*    
!*      MAJOR RUNTIME MODES
!*    
!*    
!******************************************************************************    
            
        call setState(COLUMNAR,.true.)
        call setState(LOSSY,.false.)
        call setState(TWOBEAM,.false.)
        call setState(GAUGETRANSFORM,.false.)
        call setState(OPXYZ,.false.)
        call setState(OPXYZPLUS,.false.)
        call setState(OPPNG,.true.)
        call setState(OPDIFF,.false.)


!******************************************************************************    
!*    
!*      COMMAND LINE ARGUMENTS
!*    
!*    
!******************************************************************************    
    
        
        
    !---    read command line arguments
        cla = CommandLineArguments_ctor(40)  
        call setCategories(cla,(/   "file handling         ",       &
                                    "microscope settings   ",       &
                                    "diffraction conditions",       &
                                    "calculation           ",       &
                                    "output                " /))
        call setProgramDescription( cla, "run_SimulatedTEM" )
        call setProgramVersion( cla, VERSION )   
          
    
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"            input filename",CLA_FILE )           
        ii=0 ; call get( cla,"a0",a0_in,ii ,LIB_CLA_OPTIONAL,"         lattice parameter(s)",CLA_FILE )
        if (ii==1) a0_in(2:3) = a0_in(1)
        call get( cla,"lattice",latticename ,LIB_CLA_OPTIONAL,"    lattice type",CLA_FILE )            
        outfile = filename                                       
        call get( cla,"o",outfile ,LIB_CLA_OPTIONAL,"          output filename prefix",CLA_FILE )           
        if ( (trim(getSuffix(outfile))=="xyz").or.(trim(getSuffix(outfile))=="lammps").or.(trim(getSuffix(outfile))=="lmp").or.(trim(getSuffix(outfile))=="dat").or.(trim(getSuffix(outfile))=="png") ) then
            outfile = trim(removeSuffix(outfile))
        end if
        
        
        call get( cla,"T",T ,LIB_CLA_OPTIONAL,"          temperature (K)",CLA_MICRO )                                                        
        call get( cla,"V",V ,LIB_CLA_OPTIONAL,"          accelerating voltage (kV)",CLA_MICRO )                                                        
        call get( cla,"theta",tilt_max,LIB_CLA_OPTIONAL,"      maximum foil tilt angle (deg)",CLA_MICRO )           
        call get( cla,"ntilt",ntilt,LIB_CLA_OPTIONAL,"      number of foil tilt angles considered",CLA_MICRO )           
        call get( cla,"tweak",tweak,LIB_CLA_OPTIONAL,"      fine-tune tweak angle as fraction of theta",CLA_MICRO )           
        ii=9; call get( cla,"R",Rmat,ii,LIB_CLA_OPTIONAL,"          specify tilt matrix",CLA_MICRO )           
        ii=9; call get( cla,"U",Umat,ii,LIB_CLA_OPTIONAL,"          specify crystal orientation in input file - if x axis is [uvw] then U=ux,uy,uz,vx...wz",CLA_MICRO )           
        call get( cla,"grain",select_grain,LIB_CLA_OPTIONAL,"      assuming column 4 of .xyz is grain number, 5:13 are deformation gradient, choose orientation of one grain",CLA_MICRO )           


        ii=0; dat = LIB_CLA_NODEFAULT_I; call get( cla,"g",dat,ii ,LIB_CLA_REQUIRED,"            g-vector reflection to output",CLA_DIFF )
        if ( (ii==3).or.(ii==4) ) then
            allocate(hkl(ii))
            hkl(1:ii) = dat(1:ii)       
        end if
        call get( cla,"ng",n_g ,LIB_CLA_OPTIONAL,"         tilt foil until this reflection bright",CLA_DIFF )            
        ok = isState(FIBONACCI_PREC)
        call get( cla,"Fibonacci",ok ,LIB_CLA_OPTIONAL,"  use Fibonacci spiral to generate precession angles",CLA_DIFF )                                                        
        if (hasArgument(cla,"Fibonacci")) call setState(FIBONACCI_PREC,ok) 
        call get( cla,"nPrecAngle",nPrec ,LIB_CLA_OPTIONAL," number of precession images to take",CLA_DIFF )                                                        
        call get( cla,"precAngle",precAngle ,LIB_CLA_OPTIONAL,"  precession angle (rad)",CLA_DIFF )             
        !df = (n_g /= 1.0d0)
        call get( cla,"dark",df ,LIB_CLA_OPTIONAL,"       output dark field image",CLA_DIFF )                 
        call get( cla,"bright",bf ,LIB_CLA_OPTIONAL,"     output bright field image (set -bf or -df not both!)",CLA_DIFF )                 
        if (hasArgument(cla,"bright")) df = .not. bf        
        call get( cla,"all",all_beams ,LIB_CLA_OPTIONAL,"        output image using all beams except [000] ( requires many beam )",CLA_DIFF )      
        call get( cla,"com",com_beams ,LIB_CLA_OPTIONAL,"        output image using com of all beams except [000] ( requires many beam )",CLA_DIFF )      
        


    !---    set major runtime attributes        
        attset = isState(COLUMNAR)
        call get( cla,"col",attset ,LIB_CLA_OPTIONAL,"        use columnar approximation",CLA_CALC )                 
        if (hasArgument(cla,"col")) call setState(COLUMNAR,attset)

        attset = isState(LOSSY)
        call get( cla,"loss",attset ,LIB_CLA_OPTIONAL,"       use complex crystal structure factor for inelastic electron scattering",CLA_CALC )    
        if (hasArgument(cla,"loss")) call setState(LOSSY,attset)

        attset = isState(TWOBEAM)
        call get( cla,"2",attset ,LIB_CLA_OPTIONAL,"          compute in two beam approximation",CLA_CALC )                 
        if (hasArgument(cla,"2")) call setState(TWOBEAM,attset)
 
                    


        call get( cla,"imin",imin,LIB_CLA_OPTIONAL,"       intensity threshold for beam selection",CLA_CALC )         
        call get( cla,"zlayers",zlayers,LIB_CLA_OPTIONAL,"    number of z-layers in rho calculation block - reduce to lower memory footprint",CLA_CALC )         
        call get( cla,"M",M,LIB_CLA_OPTIONAL,"          number of image subdivisions per lattice parameter",CLA_CALC )         

        
        attset = isState(GAUGETRANSFORM)
        call get( cla,"gauge",attset ,LIB_CLA_OPTIONAL,"      use gauge transformation (deprecated)",CLA_CALC )                 
        if (hasArgument(cla,"gauge")) call setState(GAUGETRANSFORM,attset)

        call get( cla,"tilt",tiltonly ,LIB_CLA_OPTIONAL,"       tilt foil then end - no image generated",CLA_CALC )                 
        
 
        attset = isState(OPXYZPLUS)
        call get( cla,"opxyz+",attset ,LIB_CLA_OPTIONAL,"     output phase field",CLA_OUTPUT )                 
        if (hasArgument(cla,"opxyz+")) call setState(OPXYZPLUS,attset)

        attset = isState(OPXYZ)
        call get( cla,"opxyz",attset ,LIB_CLA_OPTIONAL,"      output atom positions after tilt (serial only)",CLA_OUTPUT )                 
        if (hasArgument(cla,"opxyz")) call setState(OPXYZ,attset)
        if (isState(OPXYZ).and.(nProcs>1)) then
            if (rank==0) print *,"run_SimulatedTEM warning - asked for -opxyz but >1 procs. Ignoring"
            call setState(OPXYZ,.false.)
        end if
            

        attset = isState(OPPNG)
        call get( cla,"png",attset ,LIB_CLA_OPTIONAL,"        output result as .png file",CLA_OUTPUT )                 
        if (hasArgument(cla,"png")) call setState(OPPNG,attset)

        attset = isState(OPDIFF)
        call get( cla,"opDP",attset ,LIB_CLA_OPTIONAL,"       output diffraction pattern",CLA_OUTPUT )                 
        if (hasArgument(cla,"opDP")) call setState(OPDIFF,attset)

        call get( cla,"blur",pngBlur ,LIB_CLA_OPTIONAL,"       blur .png file by smoothing width sigma",CLA_OUTPUT )                                                        
        call get( cla,"norm",pngNorm ,LIB_CLA_OPTIONAL,"       normalise .png file intensity",CLA_OUTPUT )                                                        
        call get( cla,"yflip",yflip ,LIB_CLA_OPTIONAL,"      output png with y=0 at bottom, to match ovito orientation z out of screen",CLA_OUTPUT )                 
        call get( cla,"dbg",dbg ,LIB_CLA_OPTIONAL,"        additional debug output",CLA_OUTPUT )                 
        if (rank==0) call report(cla)
        if (hasHelpArgument(cla)) call errorExit("done")
        if (.not. allRequiredArgumentsSet(cla)) call errorExit("error - required arguments unset")



    !---    check for input sanity        
        if (hasArgument(cla,"imin") .and. isState(TWOBEAM)) then
            if (rank==0) then
                print *,"run_SimulatedTEM warning - argument conflict - attempt to set -2 and -imin parameters"
                print *,"run_SimulatedTEM warning - overriding two beam mode"
            end if
            call setState(TWOBEAM,.false.)
        end if
        if ((rank==0) .and. .not. hasArgument(cla,"f")) print *,"error - expected argument '-f <filename>'"
        if ((rank==0) .and. .not. allocated(hkl)) print *,"error - expected argument '-g h,k,l' or '-g h,i,k,l'"


        if (hasArgument(cla,"bright") .and. hasArgument(cla,"dark")) then
            if (rank==0) print *,"warning - used both -bf and -df"
        end if



        zlayers = max(2,zlayers)
        if (mod(zlayers,2) == 1) zlayers = zlayers+1



        if (all_beams .and.  isState(TWOBEAM)) then
            if (rank==0) print *,"run_SimulatedTEM warning - asked for 'all' beams in two beam mode. Ignoring '-all'"
            all_beams = .false.
        end if

        
        if (com_beams .and.  isState(TWOBEAM)) then
            if (rank==0) print *,"run_SimulatedTEM warning - asked for 'com' image in two beam mode. Ignoring '-all'"
            com_beams = .false.
        end if

        if (.not. (com_beams.or.all_beams) ) then
            specifiedRmat = hasArgument(cla,"R")
        else if (.not. specifiedRmat) then
            Rmat = (/1,0,0,0,1,0,0,0,1/)
            specifiedRmat = .true.
            if (rank==0) print *,"run_SimulatedTEM info - asked for 'all beams' or 'com image', specifying tilt -R 1,0,0,0,1,0,0,0,1 "
        end if

        if (.not. hasArgument(cla,"o")) then

            if (isState(COLUMNAR)) outfile = trim(outfile)//".c"
            if (.not. isState(LOSSY)) outfile = trim(outfile)//".e"
            if (isState(TWOBEAM))  outfile = trim(outfile)//".2"
            if (nPrec>1) then
                if (isState(FIBONACCI_PREC)) then
                    outfile = trim(outfile)//".f"
                else
                    outfile = trim(outfile)//".p"
                end if
            end if
            if (com_beams) outfile = trim(outfile)//".com"
            if (all_beams) outfile = trim(outfile)//".all"
            if (.not. df) outfile = trim(outfile)//".bf"
            if (specifiedRmat) then
                outfile = trim(outfile)//".fixR"
            else
                outfile = trim(outfile)//".autoR"
            end if

            if (.not. (com_beams.or.all_beams)) then
                write(aa,fmt='(i2)') hkl(1)
                outfile = trim(outfile)//".g"//trim(adjustl(aa))
                write(aa,fmt='(i2)') hkl(2)
                outfile = trim(outfile)//trim(adjustl(aa))
                write(aa,fmt='(i2)') hkl(3)
                outfile = trim(outfile)//trim(adjustl(aa))
                if (size(hkl)==4) then
                    write(aa,fmt='(i2)') hkl(4)
                    outfile = trim(outfile)//trim(adjustl(aa))
                end if
                write(aa,fmt='(f4.2)') n_g
                outfile = trim(outfile)//".ng"//trim(adjustl(aa))
            end if
        

        
        end if
        
        opxyz_file_prefix = outfile
         
        

        specifiedUmat = hasArgument(cla,"U")



        call delete(cla)
        

!******************************************************************************    
!*    
!*      WELCOME
!*    
!*    
!******************************************************************************    
        
        if (rank==0) then
            print *,"run_SimulatedTEM v."//trim(adjustl(VERSION))
            print *,"^^^^^^^^^^^^^^^^^^^"//repeat("^",len_trim(adjustl(VERSION)))
            call get_environment_variable("HOSTNAME",dummy)
            if (len_trim(dummy)==0) then
                call get_environment_variable("HOST",dummy)
                if (len_trim(dummy)==0) dummy="localhost"           !   ??
            end if
            print *,"   running on ",nProcs," processors on "//trim(dummy)
            print *,""
            print *,"   file handling"
            print *,"   ^^^^^^^^^^^^^"
            print *,"       input atom positions    """//trim(filename)//""""
            print *,"       output file prefix      """//trim(outfile)//""""
            if (trim(latticename)==UNKNOWN_LATTICE) then
                print *,"       lattice type            TBD"
            else
                print *,"       lattice type            """//trim(latticename)//""""
            end if
            if (any(a0_in<=0)) then
                print *,"       indicative latt param   TBD"
            else
                write(*,fmt='(a,3f12.8)') "        indicative latt param   ",a0_in
            end if
            !if (specifiedUmat) then
                print *,"       specified crystal orientation"
                write(*,fmt='(a,3f12.8)') "          ",Umat(1:7:3)
                write(*,fmt='(a,3f12.8)') "          ",Umat(2:8:3)
                write(*,fmt='(a,3f12.8)') "          ",Umat(3:9:3)
            !end if
            if (select_grain/=LIB_CLA_NODEFAULT_I) &
                print *,"       select grain orientation ",select_grain

            print *,""

            print *,"   microscope settings"
            print *,"   ^^^^^^^^^^^^^^^^^^^"
            write(*,fmt='(a,f12.3,a)') "        accelerator voltage     ",V," (kV)"
            write(*,fmt='(a,f12.3,a)') "        temperature             ",T," (K)"
            if (specifiedRmat) then
                print *,"       specified tilt"
                write(*,fmt='(a,3f12.8)') "       ",Rmat(1:7:3)
                write(*,fmt='(a,3f12.8)') "       ",Rmat(2:8:3)
                write(*,fmt='(a,3f12.8)') "       ",Rmat(3:9:3)
            else
                write(*,fmt='(a,f12.3,a)') "        max foil tilt angle     ",tilt_max," (deg)"
                print *,"       number of foil tilts    ",ntilt
                write(*,fmt='(a,f12.6,a)') "        tweak foil tilt         ",tweak*tilt_max," (deg)"
            end if
            print *,""

            print *,"   diffraction condition"
            print *,"   ^^^^^^^^^^^^^^^^^^^^^"
            write(*,fmt='(a,4i3,a)') "        output g-vector         ",hkl
            if (isState(TWOBEAM)) print *,"       two beam approximation"
            print *,"       number precession imgs  ",nPrec
            if (nPrec>1) &
            write(*,fmt='(a,f12.6,a)') "        precession angle        ",precAngle," (rad)"
            write(*,fmt='(a,f12.3,a)') "        bright reflection n_g = ",n_g
            if (df) then
                print *,"       output dark field image"
            else
                print *,"       output bright field image (hkl = [000])"
            end if
            print *,""

            print *,"   calculation settings"
            print *,"   ^^^^^^^^^^^^^^^^^^^^"
            call reportDashboard(o=8)
            print *,"       tilt then quit          ",tiltonly
            if (.not. isState(TWOBEAM)) &
            print *,"       beam intensity thresh   ",imin   
            print *,"       size of rho calc block  ",zlayers
            print *,""

            print *,"   output settings"
            print *,"   ^^^^^^^^^^^^^^^"
            if (isState(OPPNG)) &
            print *,"       output .png  blur       ",pngBlur," normalised? ",pngNorm
            print *,"       yflip .png              ",yflip
            print *,"       debug mode              ",dbg
            print *,""  

        end if





!******************************************************************************    
!*    
!*      CONSTRUCTOR
!*    
!*    
!******************************************************************************    

        if (rank==0) then
            print *,""
            print *,"initialising"
            print *,"^^^^^^^^^^^^"
        end if

        if (dbg) then
            LIB_IMB_DBG = dbg
            Lib_ComputePhaseFactor_DBG = dbg
            LIB_MANYBEAM_DBG = dbg
            Lib_ImagingSpace_DBG = dbg
        end if

        call Lib_IntegrateManyBeams_init_MPI()
        !call pause(tt(T_INIT))


    !---    constructor - read input file and determine parameters of calculation
        !tt(T_CTOR) = Callipers_ctor()
        if (specifiedUmat) then
            imb = IntegrateManyBeams_ctor( latticename,a0_in,filename,T,V,M, nPrec = nPrec, precAngle = precAngle, rho_in = max(4.0d0,norm2(n_g*hkl)+4), zlayers = zlayers, dfg_bar = reshape( Umat,(/3,3/))  )
        else if (select_grain/=LIB_CLA_NODEFAULT_I) then
            imb = IntegrateManyBeams_ctor( latticename,a0_in,filename,T,V,M, nPrec = nPrec, precAngle = precAngle, rho_in = max(4.0d0,norm2(n_g*hkl)+4), zlayers = zlayers, select_grain = select_grain )
        else
            imb = IntegrateManyBeams_ctor( latticename,a0_in,filename,T,V,M, nPrec = nPrec, precAngle = precAngle, rho_in = max(4.0d0,norm2(n_g*hkl)+4), zlayers = zlayers )
        end if
        call pause(tt(T_INIT))


        
!******************************************************************************    
!*    
!*      FIND GOOD DIFFRACTION CONDITIONS
!*    
!*    
!******************************************************************************    

        
        if ((rank==0) .and. isState(OPDIFF)) then
            call outputDiffractionPattern( imb,trim(outfile)//"_before",1024,maxIntensity = .false.,yflip=yflip )
        end if

    !---    tilt the foil to "good" diffraction conditions
        
        if (rank==0) then
            print *,""
            print *,"finding diffraction condition"
            print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
        end if
        tt(T_TILT) = Callipers_ctor()        
        call tiltFoil(imb,tilt_max*PI/180,hkl,n_g,df)
        call pause(tt(T_TILT))


        tt(T_OUTPUT) = Callipers_ctor()
        call pause(tt(T_OUTPUT))
        if ((rank==0) .and. isState(OPDIFF)) then
            call start(tt(T_OUTPUT))
            call outputDiffractionPattern( imb,trim(outfile)//"_ready",1024,maxIntensity = .false.,yflip=yflip )
            call outputDiffractionPattern( imb,trim(outfile)//"_max",1024,maxIntensity = .true.,yflip=yflip )
            call pause(tt(T_OUTPUT))
        end if
        
        call start(tt(T_TILT))
        if (isState(TWOBEAM)) then
            call changeGvectors( imb , hkl_in = hkl ) 
        else
            call selectImportantBeams( imb, imin = imin, hkl = hkl_important )
            call changeGvectors( imb , rule = LIB_IMB_BLOCK_HIGH_SG , hkl_in = hkl_important )
        end if

        call setImagingSpace(imb)                        
        call pause(tt(T_TILT))


    ! !---    compute the phase field x = exp[ - ig.u(r) ]        
    !     tt(T_PHASE) = Callipers_ctor()
        
    !     call pause(tt(T_PHASE))



        if (rank==0) then
            print *,""
            print *,"run parameters"
            print *,"^^^^^^^^^^^^^^"
            print *,""
            call report(imb)
            print *,""
        end if


!******************************************************************************    
!*    
!*      DUMP .xyz file pof atom positions if requested
!*    
!*    
!******************************************************************************    


        !---    dump of atom positions and z-slice
            if (isState(OPXYZ) .and. (rank==0)) then   
                print *,"run_SimulatedTEM info - writing atom positions """//trim(opxyz_file_prefix)//".xyz"//""""
                open(unit=500,file=trim(opxyz_file_prefix)//".xyz",action="write")
                    as = getAtomSpace(imb)
                    call getBounds(getImagingSpace(imb),bb)
                    dd = 1/getia(as)
                    a_super = 0
                    ix = getNbuf( imb )     !   buffer used to ensure image is correctly generated
                    !print *,"ix = ",ix
                    a_super(1,1) = (bb(2,1) - 2*ix)*dd
                    a_super(2,2) = (bb(2,2) - 2*ix)*dd
                    a_super(3,3) = (bb(2,3) - 2*ix)*dd
                    !print *,bb(:,1) 
                    !print *,bb(:,2) 
                    iy = 0 ; Nx = bb(2,1) + 1 ; Ny = bb(2,2) + 1 
                    do ii = 1,getnAtoms(imb)
                        rr = getAtomPosition(imb,ii) 
                        if ( (rr(1)>=ix).and.(rr(1)<=Nx-1-ix).and.(rr(2)>=ix).and.(rr(2)<=Ny-1-ix) ) iy = iy + 1
                    end do

                    write(unit=500,fmt='(i8)') iy !getnAtoms(imb)
                    write(unit=500,fmt='(a,9f12.6,a,3f12.6,a)') "Lattice=""",a_super,""" Origin=""",(/ix,ix,ix/)*dd,""" Properties=Species:S:1:Pos:R:3"
                    do ii = 1,getnAtoms(imb)
                        rr = getAtomPosition(imb,ii) 
                        if ( (rr(1)>=ix).and.(rr(1)<=Nx-1-ix).and.(rr(2)>=ix).and.(rr(2)<=Ny-1-ix) )    &
                        write(unit=500,fmt='(a,10f16.8)') getElement(imb),rr*dd
                    end do
                close(unit=500)
            end if






        if (.not. tiltonly) then

!******************************************************************************    
!*    
!*      INTEGRATE
!*    
!*    
!******************************************************************************    

        

        if (rank==0) then
            print *,""
            print *,"propagating beams"
            print *,"^^^^^^^^^^^^^^^^^"
        end if
        tt(T_INTEGRATE) = Callipers_ctor()
            
        call integrate( imb )
        call pause(tt(T_INTEGRATE))
        
!******************************************************************************    
!*    
!*      OUTPUT
!*    
!*    
!******************************************************************************    


        
        
        if (isState(OPDIFF)) then
            call start(tt(T_OUTPUT))
            call outputDiffractionPattern( imb,trim(outfile)//"_after",1024,maxIntensity=.false.,yflip=yflip )
            call pause(tt(T_OUTPUT))
        end if

        ! if (opDiffPatt) then
        ! !    call outputDiffractionPattern( imb,trim(outfile),1024 )
        ! end if

        if (isState(OPPNG)) then
            call start(tt(T_OUTPUT))
            if (rank==0) then
                print *,""
                print *,"output .png"
                print *,"^^^^^^^^^^^"
            end if
                        
            if (com_beams) then
                if (rank==0) print *,"run_SimulatedTEM info - generating centre-of-mass diffraction pattern image"
                call getCentreOfMassImage( imb,com_img,df )    !   gathers result on rank 0
            else if (df) then
                if (all_beams) then
                    if (rank==0) print *,"run_SimulatedTEM info - generating image using all beams"
                    call getIntensity( imb,img )            !   gathers result on rank 0
                else
                    if (rank==0) print *,"run_SimulatedTEM info - generating image using selected beam"
                    call getIntensity( imb,hkl,img )        !   gathers result on rank 0
                end if
            else
                hkl = 0
                if (rank==0) print *,"run_SimulatedTEM info - generating image using hkl=[000]"
                call getIntensity( imb,hkl,img )            !   gathers result on rank 0
            end if

            if (rank==0) then

                if (com_beams) then

                    Nx = size(com_img,dim=2)
                    Ny = size(com_img,dim=3)
                    print *,"run_SimulatedTEM info - image size ",Nx,"x",Ny 
                    write(*,fmt='(a,3f10.6)') " run_SimulatedTEM info - writing to """//trim(outfile)//".png"""

                    if (pngBlur) then
                        dd = getSigma( imb )/getA( imb )
                        do ii = 1,3
                            call blurImg( com_img(ii,:,:), dd ) 
                        end do
                    end if

                    if (yflip) then
                        !   flip image in y-axis
                        do iy = 0,int( (Ny-1)/2 )      !   round down
                            iyp = Ny-1-iy
                            do ix = 0,Nx-1
                                do ii = 1,3
                                    dd = com_img(ii,ix,iy)
                                    com_img(ii,ix,iy) = com_img(ii,ix,iyp)
                                    com_img(ii,ix,iyp) = dd
                                end do
                            end do
                        end do
                    end if

                    ix = getNbuf( imb )     !   buffer used to ensure image is correctly generated
                    call write_rgb_png( trim(outfile)//".png",com_img(:,ix:Nx-1-ix,ix:Ny-1-ix) )

                else

                    Nx = size(img,dim=1)
                    Ny = size(img,dim=2)
                    print *,"run_SimulatedTEM info - image size ",Nx,"x",Ny 
                    write(*,fmt='(a,3f10.6)') " run_SimulatedTEM info - minmaxavg img ",minval(img),maxval(img),sum(img)/(Nx*Ny)
                    write(*,fmt='(a,3f10.6)') " run_SimulatedTEM info - writing to """//trim(outfile)//".png"""

                    if (pngBlur) then
                        dd = getSigma( imb )/getA( imb )
                        call blurImg( img, dd ) 
                    end if

                    if (yflip) then
                        !   flip image in y-axis
                        do iy = 0,int( (Ny-1)/2 )      !   round down
                            iyp = Ny-1-iy
                            do ix = 0,Nx-1
                                dd = img(ix,iy)
                                img(ix,iy) = img(ix,iyp)
                                img(ix,iyp) = dd
                            end do
                        end do
                    end if

                    ix = getNbuf( imb )     !   buffer used to ensure image is correctly generated

                    call writePng( trim(outfile)//".png",img(ix:Nx-1-ix,ix:Ny-1-ix),normalise=pngNorm )

                end if


            end if
            call pause(tt(T_OUTPUT))
        end if
        
        else

            if (isState(OPDIFF)) then
                call start(tt(T_OUTPUT))
                call outputDiffractionPattern( imb,trim(outfile)//"_after",1024,maxIntensity=.false.,yflip=yflip )
                call pause(tt(T_OUTPUT))
            end if            
            if (rank==0) print *,"run_SimulatedTEM info - compute tilt then quit"
        end if      !   tiltonly

!******************************************************************************    
!*    
!*      BYE BYE
!*    
!*    
!******************************************************************************    

        if (rank==0) then
            print *,""
            print *,"timing data (s)"
            print *,"^^^^^^^^^^^^^^^"
            write(*,fmt='(a,f16.6)') "   initialisation ",elapsed(tt(T_INIT))
            write(*,fmt='(a,f16.6)') "   tilt foil      ",elapsed(tt(T_TILT))
            if (.not. tiltonly) write(*,fmt='(a,f16.6)') "   integration    ",elapsed(tt(T_INTEGRATE))
            write(*,fmt='(a,f16.6)') "   output         ",elapsed(tt(T_OUTPUT))
            write(*,fmt='(a,f16.6)') "   total          ",elapsed(tt(T_TOTAL))
        end if

        if (rank==0) then
            print *,""
        end if
        call errorExit("done")
         

               
    contains
!---^^^^^^^^
!        
        
        
        subroutine blurImg(img,sigmaona)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      perform a simple 2d blur operation on an image.
    !*      This can be useful in low contrast images, where tiny fluctuations in 
    !*      atom density or grad x lead to (apparent) differences in the output image.
    !*      blur radius is sigma/a pixels.
            real(kind=real64),dimension(0:,0:),intent(inout)        ::      img
            real(kind=real64),intent(in)                            ::      sigmaona      !   blurring
            
            integer             ::      Nx,Ny
            integer             ::      ix,iy,ik,kernel_width
            real(kind=real64),dimension(:),allocatable      ::      img_tmp
            real(kind=real64),dimension(:),allocatable      ::      kernel
            real(kind=real64)   ::      ff,ww,i2s2

            Nx = size(img,dim=1)
            Ny = size(img,dim=2)
            allocate(img_tmp(0:max(Nx,Ny)-1))
            
            write(*,fmt='(a,f8.3,a)') " run_SimulatedTEM::blurImg info - blur ",sigmaona," (px)"
            kernel_width = ceiling(3*sigmaona)            
            allocate(kernel(-kernel_width:kernel_width))
            i2s2 = 1/(2*sigmaona*sigmaona)
            do ik = -kernel_width,kernel_width
                kernel(ik) = exp( -ik*ik*i2s2 )
            end do

        !---    blur in x direction for each row iy
            do iy = 0,Ny-1
                img_tmp(0:Nx-1) = img(0:Nx-1,iy)
                do ix = 0,Nx-1
                    ff = 0 ; ww = 0
                    do ik = max(0,ix-kernel_width),min(Nx-1,ix+kernel_width)
                        ww = ww + kernel(ik-ix)
                        ff = ff + kernel(ik-ix)*img_tmp(ik)
                    end do
                    img(ix,iy) = ff/ww
                end do                
            end do

        !---    blur in y direction for each column ix
            do ix = 0,Nx-1
                img_tmp(0:Ny-1) = img(ix,0:Ny-1)
                do iy = 0,Ny-1
                    ff = 0 ; ww = 0
                    do ik = max(0,iy-kernel_width),min(Ny-1,iy+kernel_width)
                        ww = ww + kernel(ik-iy)
                        ff = ff + kernel(ik-iy)*img_tmp(ik)
                    end do
                    img(ix,iy) = ff/ww
                end do                
            end do
            return
        end subroutine blurImg

        subroutine tiltFoil(imb,tilt_max,hkl,n_g,dark)        
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
    !*      find a good tilt angle 0:tilt_max (radians) so that the reflection
    !*      ng*[hkl] is bright. 
    !*      Optionally also tweak until 1*[hkl] is dark

            type(IntegrateManyBeams),intent(inout)              ::      imb
            real(kind=real64),intent(in)                        ::      tilt_max
            integer,dimension(:),intent(in)                     ::      hkl
            real(kind=real64),intent(in)                        ::      n_g
            logical,intent(in)                                  ::      dark
            real(kind=real64),dimension(3,3)    ::      foil_tilt

            logical                             ::      MillerBravais 
            integer,dimension(size(hkl))        ::      hkl_bright,hkl_dark
            real(kind=real64)                   ::      I_g
            real(kind=real64),dimension(3,3)    ::      foil_tilt2
            character(len=16)                   ::      form
            logical                             ::      maxIntensity = .false.      !   select tilt based on foil exit intensity
            character(len=256)                  ::      rstring
            integer                             ::      ii,jj


            !if ( dark.and.(rank==0) ) print *,"tiltFoil info - ""dark"" option not currently used in tiltFoil"
            if ( dark.and.(rank==0) ) print *,"tiltFoil info - ""dark"" option used"
            hkl_dark = 0
            if (dark) hkl_dark = hkl

            MillerBravais = (size(hkl)==4)
            form = "(a,3i4,a,f16.8)"
            if (MillerBravais) form = "(a,4i4,a,f10.6)"

            if (n_g<1) call errorExit("tiltFoil error - n_g<1" )

            if (rank==0) then
                call perfectLatticeIntensity(imb,hkl,maxIntensity,I_g)
                write(*,fmt=form) " tiltFoil info - perfect lattice intensity ",hkl," before tilt ",I_g
            end if

            if (tilt_max<=0) return
 
            if (specifiedRmat) then

                foil_tilt = reshape( Rmat , (/3,3/) )
                call restoreUnitaryness( foil_tilt )
                call applyFoilTilt( imb,foil_tilt )
                if (rank==0) then
                    ! if (n_g - floor(n_g) < 0.5) then
                    !     hkl_bright = floor( hkl * n_g )
                    ! else
                    !     hkl_bright = ceiling( hkl * n_g )
                    ! end if
                    hkl_bright = nint( n_g ) * hkl
                    call perfectLatticeIntensity(imb,hkl_bright,maxIntensity,I_g)
                    write(*,fmt=form) " tiltFoil info - perfect lattice intensity at ",hkl_bright," after user specified tilt ",I_g
                end if

            else
        
            !    if (.not. ( (n_g==1.0).and.df) ) then
                    if (abs(n_g - nint(n_g)) < 1.0d-6) then
                        if (rank==0) then
                            print *,""
                            if (dark) write(*,fmt=form) " selecting dark reflection g,n_g = ",hkl,",",1.0d0
                            write(*,fmt=form) " selecting bright reflection g,n_g = ",hkl,",",n_g
                            print *,""
                        end if
                        hkl_bright = nint( n_g ) * hkl
                        ! if (dark) then
                        !     call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT, foil_tilt, hkl_dark )
                        ! else
                            call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT, foil_tilt )
                        !end if
                    else
                        if (rank==0) then
                            print *,""
                            if (dark) write(*,fmt=form) " selecting dark reflection g,n_g = ",hkl,",",1.0d0
                            write(*,fmt=form) " selecting bright reflections for g,n_g = ",hkl,",",n_g
                            print *,""
                        end if
                        hkl_bright = floor( n_g ) * hkl
                        ! if (dark) then
                        !     call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT, foil_tilt, hkl_dark )
                        ! else
                            call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT, foil_tilt )
                        !end if
                        hkl_bright = ceiling( n_g ) * hkl
                        ! if (dark) then
                        !     call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT, foil_tilt2, hkl_dark )
                        ! else
                            call selectFoilTilt( imb,hkl_bright,tilt_max,NTILT, foil_tilt2 )
                        !end if
                        
                        foil_tilt = quaternionToRotMat( slerp( Quaternion_ctor(foil_tilt),Quaternion_ctor(foil_tilt2), ceiling(n_g) - n_g ) )       !   ceiling(x)-x = 1 for x = 3.00 and = 0 for 3.99

                        ! if (n_g - floor(n_g) < 0.5) then
                        !     hkl_bright = floor( n_g ) * hkl
                        ! else
                        !     hkl_bright = ceiling( n_g ) * hkl
                        ! end if
                        hkl_bright = nint( n_g ) * hkl

                    end if            
                    if (rank==0) then
                        call perfectLatticeIntensity(imb,hkl_bright,maxIntensity,I_g)
                        write(*,fmt=form) " tiltFoil info - perfect lattice intensity ",hkl_bright," before tilt ",I_g
                    end if
                    if (rank==0) print *,""
                    call applyFoilTilt(imb,foil_tilt)
                    if (rank==0) then
                        call perfectLatticeIntensity(imb,hkl_bright,maxIntensity,I_g)
                        write(*,fmt=form) " tiltFoil info - perfect lattice intensity at ",hkl_bright," after tilt ",I_g
                    end if
            !    end if
                
                 if (dark) then
                ! !   add an additional tweak to make hkl extra dark
                        if (dark) then
                            call selectFoilTilt( imb,hkl_bright,tilt_max*TWEAK,NTILT, foil_tilt2, hkl_dark )
                        else
                            call selectFoilTilt( imb,hkl_bright,tilt_max*TWEAK,NTILT, foil_tilt2 )
                        end if

!                     call selectFoilTilt( imb,hkl,tilt_max*TWEAK,NTILT,.false., foil_tilt )
                     call applyFoilTilt(imb,foil_tilt2)
                     foil_tilt = matmul(foil_tilt,foil_tilt2)
                     if (rank==0) then
                         call perfectLatticeIntensity(imb,hkl_bright,maxIntensity,I_g)
                         write(*,fmt=form)  " tiltFoil info - perfect lattice intensity ",hkl_bright," after tweak  ",I_g
                         if (dark) then
                            call perfectLatticeIntensity(imb,hkl_dark,maxIntensity,I_g)
                            write(*,fmt=form)  " tiltFoil info - perfect lattice intensity ",hkl_dark," after tweak  ",I_g
                         end if
                     end if
                end if

            end if                


            if ( (rank==0).and.(n_g>1) ) then
                call perfectLatticeIntensity(imb,hkl,maxIntensity,I_g)
                write(*,fmt=form)  " tiltFoil info - perfect lattice intensity at ",hkl," after tilt ",I_g
            end if

            if (rank==0) then
            !---    output the selected rotation matrix in a form that can be reentered for next calculation
                write(rstring,fmt='(9(f16.12,a1))') foil_tilt(1,1),",",foil_tilt(2,1),",",foil_tilt(3,1),      &
                                                ",",foil_tilt(1,2),",",foil_tilt(2,2),",",foil_tilt(3,2),      &
                                                ",",foil_tilt(1,3),",",foil_tilt(2,3),",",foil_tilt(3,3)
                jj = 0
                do ii = 1,len(rstring)
                    if (rstring(ii:ii)/=" ") then
                        jj = jj + 1
                        rstring(jj:jj) = rstring(ii:ii)
                    end if
                end do
                

                write(*,fmt='(a)') " tiltFoil info -R "//rstring(1:jj)
            end if

            return
        end subroutine tiltFoil

        subroutine selectImportantBeams( imb, imin, hkl )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      select a beam as "important" if its intensity ever gets above imin
    !*      allocate and return the array of reflections hkl
            type(IntegrateManyBeams),intent(in)                 ::      imb
            real(kind=real64),intent(in)                        ::      imin
            integer,dimension(:,:),allocatable,intent(out)      ::      hkl

            type(Gvectors)                  ::      gv
            real(kind=real64)               ::      i_important
            real(kind=real64),dimension(:),allocatable          ::      I_g
            logical                         ::      MillerBravais 
            integer                         ::      nG
            integer                         ::      kk,nn
            character(len=32)               ::      form
             


        !---    get the set of g-vectors            
            gv = getGvectors(imb)
            nG = getn(gv)
            MillerBravais = isMillerBravais(gv)

            if (rank==0) then
                print *,"select important beams"
                call report(gv)
            end if


        !---    find max intensity of beam inside the (assumed perfect) foil            
            allocate(I_g(0:nG))
            call getIntensity(imb,I_g,maxIntensity = .true.)

        !---    count the important beams
            nn = 0
            do kk = 0,nG
                if (I_g(kk)>=imin) nn = nn + 1
            end do

        !---    is this a sensible number?             
            if (nn<2) then
                kk = maxloc( I_g(1:),dim=1 )            !   given that I_g(0) = I_g(hkl=[000]) is important, what is number 2?
                i_important = I_g(kk) * 0.50d0          !   I'll say that anything within 50% of the number 2 beam is important. 
                if (rank == 0) then
                    print *,"selectImportantBeams info - did not find any beams (except hkl=[000]) over minimum intensity threshold ",imin
                    if (MillerBravais) then
                        print *,"    second most important beam is [hjkl] = ",gethjkl(gv,kk)," with intensity ",I_g(kk)
                    else
                        print *,"    second most important beam is [hkl] = ",gethkl(gv,kk)," with intensity ",I_g(kk)
                    end if
                    print *,"    using intensity threshold ",i_important
                end if
            end if
                    

        !---    allocate memory for important reflections            
            if (MillerBravais) then
                allocate(hkl(4,0:nn-1))
            else
                allocate(hkl(3,0:nn-1))
            end if
            hkl = 0

        !---    store the important reflections. The first is hkl=[000], which is not in gv.
            form = "(a,i4,a,3i4,a,f16.8)"
            if (MillerBravais) form = "(a,i4,a,4i4,a,f10.6)"

            nn = 0
            if (rank == 0) write (*,fmt=form) " selectImportantBeams info - selected beam ",nn," reflection ",hkl(:,nn)," max intensity ",I_g(0)
            do kk = 1,nG
                if (I_g(kk)>=imin) then
                    nn = nn + 1
                    if (MillerBravais) then
                        hkl(:,nn) = gethjkl(gv,kk)
                    else
                        hkl(:,nn) = gethkl(gv,kk)
                    end if
                    if (rank == 0) write (*,fmt=form) " selectImportantBeams info - selected beam ",nn," reflection ",hkl(:,nn)," max intensity ",I_g(kk)
                end if
            end do
            

            return
        end subroutine selectImportantBeams

        subroutine outputDiffractionPattern( imb,prefix,Nx,maxIntensity,yflip )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      output the diffraction pattern as a .png
            type(IntegrateManyBeams),intent(inout)              ::      imb
            character(len=*),intent(in)                         ::      prefix
            integer,intent(in)                                  ::      Nx
            logical,intent(in)                                  ::      maxIntensity        !   if true, then show where the energy is, if false, show diff pat at foil exir
            logical,intent(in)                                  ::      yflip               !   should we flip the diffraction pattern to match the flip in the TEM image?
            real(kind=real64),dimension(:,:),allocatable        ::      img
            type(Gvectors)                  ::      gv
            integer                         ::      ix,iy,ii
            integer                         ::      nG
            real(kind=real64),dimension(:,:),allocatable        ::      gg
            real(kind=real64)               ::      maxg,ming
            real(kind=real64)               ::      spotRadius,ss,I_0,I_x,dd
            real(kind=real64),dimension(2)  ::      dx
            real(kind=real64),dimension(:),allocatable          ::      I_g
            real(kind=real64),parameter     ::      MIN_GREYSCALE = 1.0d0/(256*256)
            logical                         ::      hide_000

        !---    get the g vectors
            gv = getGvectors(imb)
            nG = getn(gv)
            ii = whichg(gv,(/0.0d0,0.0d0,0.0d0/))       !   is [000] one of the set?
            if (ii == LIB_GVECTORS_UNSET) then
                allocate(gg(3,0:nG))
                gg(:,0) = 0
                do ii = 1,ng
                    gg(:,ii) = getG(gv,ii)
                end do
            else
                allocate(gg(3,1:nG))                
                do ii = 1,ng
                    gg(:,ii) = getG(gv,ii)
                end do
            end if


        !---    find the smallest/largest magnitude g-vector
            maxg = -huge(1.0)
            do ii = lbound(gg,dim=2),nG
                maxg = max(maxg,norm2(gg(:,ii)))
            end do
            ming = +huge(1.0)
            do ii = 1,nG
                ming = min(ming,norm2(gg(:,ii)))
            end do


        !---    find the intensity of the beams at the foil exit            
            allocate(I_g(0:nG))
            call getIntensity(imb,I_g,maxIntensity)
            if ((rank==0).and.dbg) then
                do ii = lbound(gg,dim=2),nG
                    print *,"spot ",ii,gethkl(gv,ii),gg(:,ii),I_g(ii)
                end do
            end if

        !---    find a decent scaling of the image
        !    spotRadius = Nx * sqrt(0.1d0/(nG*PI))            
        !    ss = ( Nx/2 - 2*spotRadius )/maxg
        !    spotRadius = min( spotRadius , maxg*ss/4 )

            ss = (Nx/2) / ( maxg + ming/2 )     !   scale such that 1/2 pixel width will contain maxg + 2x radius of spot
            spotRadius = ss* ming/4                 !   want closest spot separation = 2x radius


            if (rank==0 .and. dbg) print *,"run_SimulatedTEM::outputDiffractionPattern info - smallest/largest |g| ",ming,maxg," spotRadius ",spotRadius," scaling ",ss
            if (rank==0) print *,"run_SimulatedTEM::outputDiffractionPattern info - minmax I_g ",minval(I_g),maxval(I_g)

        !---    change the brightness of the central spot to not overwhelm image
            I_0 = I_g(0)
            I_g(0) = 0
            I_x = maxval(I_g)
            if ( I_0 > 2*I_x ) then
                I_g(0) = 2*I_x
                hide_000 = .true.
                if (rank==0) print *,"run_SimulatedTEM::outputDiffractionPattern info - reset I(000)=",I_0," to ",2*I_x," to limit dynamic range "
            else
                I_g(0) = I_0
                hide_000 = .false.
            end if



        !---    construct the image            
            allocate(img(0:Nx-1,0:Nx-1))
            img = 0
            do iy = 0,Nx-1
                do ix = 0,Nx-1

                !---    spot 000 special
                    dx(1:2) = (/ix,iy/) - Nx/2  
                    dd = norm2(dx)
                    if (hide_000) then
                        if ( (dd<=spotRadius).and.(dd>=spotRadius-2) )  img(ix,iy) = I_g(0) 
                    else 
                        if (dd<=spotRadius)   img(ix,iy) = I_g(0)   
                    end if



                    do ii = 1,nG
                        dx(1:2) = (/ix,iy/) - Nx/2 - gg(1:2,ii)*ss 
                        dd = norm2(dx)
                        if (I_g(ii) >= MIN_GREYSCALE) then
                            if (dd<=spotRadius)   img(ix,iy) = I_g(ii) 
                        else if (I_g(ii) > 0) then
                            if (dd<=spotRadius/2) img(ix,iy) = MIN_GREYSCALE                        
                        end if
                        !if ( (norm2(dx)<=spotRadius).and. (ix==Nx/2).and.(iy==Nx/2) ) print *,"spot centre ",ii,ix,iy,dx,img(ix,iy),I_g(ii)
                    end do

                end do
            end do

        !---    write image
            if (rank==0) then
                print *,"run_SimulatedTEM::outputDiffractionPattern info - drawing ",(nG+1)," spots to """//trim(prefix)//".diff.png"//""""
                print *,"run_SimulatedTEM::outputDiffractionPattern info - minmax img ",minval(img),maxval(img)
                if (yflip) then         !   flip the orientation of the diffraction pattern
                    do iy = 0,(Nx/2)-1
                        do ix = 0,Nx-1
                            ss = img(ix,iy)
                            img(ix,iy) = img(ix,Nx-1-iy)
                            img(ix,Nx-1-iy) = ss
                        end  do
                    end do
                end if
                img = max(1.0d0/65536,min(1.0d0-1.0d0/65536,img))
                !img = img ** (1.0d0/3)      !   gamma correction
                !img = log(img)
                call writePng( trim(prefix)//".diff.png",img,normalise=.true. )
            end if

            return
        end subroutine outputDiffractionPattern



        subroutine errorExit(message)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in),optional             ::      message
#ifdef MPI            
            integer             ::      ierror
#endif
            if (rank==0) then
                if (present(message)) print *,"run_SimulatedTEM::"//trim(message)
            end if
#ifdef MPI            
            call MPI_FINALIZE(ierror)
#endif          
            if (rank==0) print *,""
            stop
        end subroutine errorExit
        
    end program run_SimulatedTEM
