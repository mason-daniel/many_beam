
    program test_Lib_CrystalStructureFactor
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        use Lib_ColouredTerminal
        use Lib_CommandLineArguments
        use Lib_Elements        
        use Lib_Lattices
        use Lib_CrystalStructureFactor
        use Lib_RelativisticElectrons
        use Lib_InputTDS
        use Lib_InputReal
        use Lib_FitDoyleTurner
        use iso_fortran_env
        implicit none


        real(kind=real64),parameter     ::      PI =  3.14159265359d0      

        logical                 ::      ok
        character(len=8),dimension(10)          ::      element = "Fe"
        real(kind=real64),dimension(10)         ::      frac = 1.0
        integer                 ::      nElements
        real(kind=real64),dimension(3)          ::     a0
        character(len=256)      ::      outfile = ""
        real(kind=real64)       ::      V = 200000
        real(kind=real64)       ::      T = 300
        integer,dimension(:,:),allocatable      ::      hkl
        integer,dimension(12000) ::      dat
        integer                 ::      nReflections = 0
        character(len=8)        ::      latticename = ""
        integer                 ::      ii,nS
        complex(kind=real64),dimension(:),allocatable   ::      Vg          !   (1:nReflections)
        complex(kind=real64),dimension(:),allocatable   ::      xig         !   (1:nReflections)
        complex(kind=real64),dimension(DT_COEFFS)       ::      aa,bb       !   fitted complex D-T coefficients
        type(CommandLineArguments)      ::      cla
        real(kind=real64)       ::      dwf,ds,ss,f_elastic,f_inelastic
        type(Lattice)           ::      latt
        integer,dimension(3)    ::      hkl_tmp        
        logical                 ::      greater
        character(len=80)       ::      elcat

        cla = CommandLineArguments_ctor(10)
        call setProgramDescription( cla, "test_Lib_CrystalStructureFactor" )
        nElements = 0
        call get( cla,"e",element,nElements  ,LIB_CLA_OPTIONAL,"        element(s)" )  

        if (nElements>1) then
            call get( cla,"f",frac,ii  ,LIB_CLA_OPTIONAL,"        atomic frac element" )           
        end if
        


        call get( cla,"latt",latticename  ,LIB_CLA_OPTIONAL,"     lattice" )  
        ii = 0
        call get( cla,"a0",a0,ii  ,LIB_CLA_OPTIONAL,"       lattice parameter" )  


        call get( cla,"V",V        ,LIB_CLA_OPTIONAL,"        accelerating voltage (V)" )  
        call get( cla,"T",T        ,LIB_CLA_OPTIONAL,"        temperature (K)" )  
        nReflections = 0
        call get( cla,"hkl",dat,nReflections,LIB_CLA_OPTIONAL,"     reflections" )          
        call get( cla,"o",outfile  ,LIB_CLA_OPTIONAL,"        output to crystal structure factor file" )  


    !  check for consistent input
        if (len_trim(latticename)==0) latticename = getLatticeName(element(1))

        nElements = max(1,nElements)

        if (hasArgument(cla,"a0")) then
            if (ii==1) then 
                a0(2:3) = a0(1)
                if (getLatticeType(latticename) == LATTICE_HCP) a0(3) = sqrt(8.0d0/3)*a0(1)
            end if
        else
            a0 = getLatticeConstant(element(1:nElements),frac(1:nElements))      
        end if


        if (getLatticeType(latticename) == LATTICE_CUSTOM) then
            print *,"test_Lib_CrystalStructureFactor error - did not recognise lattice """//trim(latticename)//""""
            print *,"   only coded"
            call listAvailableLattices()
            stop
        end if      
        Latt = Lattice_ctor(latticename)
        if (getLatticeType(latt) == LATTICE_HCP) call setCoverA(latt,a0(3)/a0(1))

        if (hasArgument(cla,"hkl")) then

            if (latticename == "hcp") then
                if (mod(nReflections,4)/=0) stop "test_Lib_CrystalStructureFactor error - expect sets of hikl Miller-Bravais indices"
                nReflections = nReflections/4
                allocate(hkl(4,nReflections))
                hkl = reshape( dat(1:4*nReflections),(/4,nReflections/) )

            !---    check for [0000] - always needed
                ok = .false.
                do ii = 1,nReflections
                    if (all(hkl(:,ii)==0)) ok = .true.
                end do
                if (.not. ok) then
                    deallocate(hkl)
                    allocate(hkl(4,nReflections+1))
                    hkl(:,1) = 0
                    hkl(:,2:nReflections+1) = reshape( dat(1:4*nReflections),(/4,nReflections/) )
                    nReflections = nReflections+1
                end if
            else
                if (mod(nReflections,3)/=0) stop "test_Lib_CrystalStructureFactor error - expect sets of hkl Miller indices"
                nReflections = nReflections/3
                allocate(hkl(3,nReflections))
                hkl = reshape( dat(1:3*nReflections),(/3,nReflections/) )
                
            !---    check for [000] - always needed
                ok = .false.
                do ii = 1,nReflections
                    if (all(hkl(:,ii)==0)) ok = .true.
                end do
                if (.not. ok) then
                    deallocate(hkl)
                    allocate(hkl(3,nReflections+1))
                    hkl(:,1) = 0
                    hkl(:,2:nReflections+1) = reshape( dat(1:3*nReflections),(/3,nReflections/) )
                    nReflections = nReflections+1
                end if
            end if

        else
            
            call permittedReflections( Latt,hkl )
            nReflections = ubound(hkl,dim=2)
        !---    bubble sort 
            do
                ok = .true.
                do ii = 1,nReflections-1
                    hkl_tmp = hkl(:,ii)
                    greater = norm2(real(hkl_tmp))>norm2(real(hkl(:,ii+1)))                                            &
                        .or. ( norm2(real(hkl_tmp))>norm2(real(hkl(:,ii+1))) .and. abs(hkl_tmp(2))>abs(hkl(2,ii+1)) )             &
                        .or. ( norm2(real(hkl_tmp))>norm2(real(hkl(:,ii+1))) .and. abs(hkl_tmp(2))==abs(hkl(2,ii+1)) .and. abs(hkl_tmp(1))>abs(hkl(1,ii+1)) )
                    if (greater) then
                        hkl(:,ii) = hkl(:,ii+1)
                        hkl(:,ii+1) = hkl_tmp
                        ok = .false.
                    end if
                end do
                if (ok) exit
            end do


        end if


        call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop
        call delete(cla)

        
    !---    compute crystal structure factors
        allocate(Vg(1:nReflections)) ; Vg = 0
        ! print *,""""//element//""""
        ! print *,hkl
        ! print *,T
        ! print *,V
         
        !call computeCrystalStructureFactors( element(1),hkl(:,1:nReflections), T,V, Vg , aa,bb)
        call computeCrystalStructureFactors( element(1:nElements),frac(1:nElements), latt,a0, hkl(:,1:nReflections), T,V, Vg , aa,bb)
        allocate(xig(1:nReflections))
        do ii = 1,nReflections
            xig(ii) = PI * HBAR * velocity( V ) / Vg(ii)
        end do

    !---    output result
        if (size(hkl,dim=1)==4) then
            write (*,fmt='(a16,4a16)') "hkil","Re(V_g)","Im(V_g)","Re(xi_g)","Im(xi_g)"
        else
             write (*,fmt='(a12,4a16)') "hkl","Re(V_g)","Im(V_g)","Re(xi_g)","Im(xi_g)"
        end if
        do ii = 1,nReflections
            if (size(hkl,dim=1)==4) then
                write (*,fmt='(4i4,4f16.6)') hkl(:,ii) ,  real(Vg(ii)),aimag(Vg(ii)),  real(xig(ii)),aimag(xig(ii))
            else
                write (*,fmt='(3i4,4f16.6)') hkl(:,ii) ,  real(Vg(ii)),aimag(Vg(ii)),  real(xig(ii)),aimag(xig(ii))
            end if
        end do

        if (nElements==1) then
            print *,""
            ns = 100
            ds = max(LIB_INPUTREAL_SMAX,LIB_INPUTTDS_SMAX)/ns
            dwf = DebyeWallerFactor(element(1),T,latticename) 
            write(*,fmt='(100a16)') "s=q/(4pi)","f_el(s)","f_inel(s)","f_DT^real(s)","f_DT^imag(s)"
            do ii = 0,ns
                ss = ii*ds
                f_elastic = scatteringAmplitude(element(1),ss)
                call getImagScatteringAmplitude(element(1),V,dwf,ss,f_inelastic)
                write(*,fmt='(100f16.8)') ss,f_elastic,f_inelastic,DoyleTurnerSum(ss,real(aa),real(bb)),DoyleTurnerSum(ss,aimag(aa),aimag(bb))
            end do
        end if
        print *,""


    !---    write crystal structure factor file?
        if (len_trim(outfile)>0) then
            elcat = ""
            do ii = 1,nElements
                elcat = trim(elcat) // trim( element(ii) )
            end do
            
            open(unit=600,file=trim(outfile),action="write")
                write(unit=600,fmt='(a)') "# test_Lib_CrystalStructureFactor v.0.0.1"
                write(unit=600,fmt='(a8,3a12,2a16,a20)') "element","g","T (K)","V (V)","Re(U_g)","Im(U_g)","xi_g"
                if (size(hkl,dim=1)==4) then
                    do ii = 1,nReflections
                        write(unit=600,fmt='(a8,3i4,2f12.3,2f16.8,f20.8)') elcat,MillerBravaisToMiller_plane(hkl(:,ii)),T,V,real(Vg(ii)),aimag(Vg(ii)),abs(xig(ii))
                    end do
                else
                    do ii = 1,nReflections
                        write(unit=600,fmt='(a8,3i4,2f12.3,2f16.8,f20.8)') elcat,hkl(:,ii),T,V,real(Vg(ii)),aimag(Vg(ii)),abs(xig(ii))
                    end do
                end if
            close(unit=600)
        end if










        ok = .true.
        
        if (ok) then
            print *,colour(GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if


    end program test_Lib_CrystalStructureFactor