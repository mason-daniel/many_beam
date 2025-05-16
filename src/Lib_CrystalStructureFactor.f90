
    module Lib_CrystalStructureFactor
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
!*      very short helper module describing the crystal structure factor
!*      in terms of Doyle-Turner ( expansion of gaussian terms ) coefficients
!*
!*      L.M. PENG,  G. REN, S. L. DUDAREV AND M. J. WHELAN
!*      Debye-Waller Factors and Absorptive Scattering Factors of Elemental Crystals
!*      Acta Cryst. (1996). A52, 456-470
!*
!*      call with
!*           Vg = crystalStructureFactor( a0, lattice, g, a,b,dwf )
!*      with
!*          a0        lattice parameter  
!*          lattice = "bcc" or "fcc" 
!*          g         g vector (reduced units eg [110])   
!*          a(:,:),b(:,:) being complex Doyle Turner coefficients for each motif point
!*          dwf       Debye Waller factor (A^-2)
  

!*      units eV/fs/A
!*
        use Lib_Elements
        use Lib_Lattices
        use Lib_RelativisticElectrons
        use Lib_inputReal
        use Lib_inputTDS
        use Lib_FitDoyleTurner
        use iso_fortran_env
        implicit none
        private

        character(len=8),public,parameter       ::      LIB_CRYSTALSTRUCTUREFACTOR_VERSION = "0.0.1"        
        real(kind=real64),parameter,private     ::      PI =  3.14159265359d0      
            
        
    !---    
    
        public      ::      computeCrystalStructureFactors
        public      ::      crystalStructureFactor

        interface       computeCrystalStructureFactors
            module procedure        computeCrystalStructureFactors0
            module procedure        computeCrystalStructureFactors1
            module procedure        computeCrystalStructureFactors2
        end interface
  
                                                                                                                                                                                
        interface           crystalStructureFactor
            module procedure            crystalStructureFactor1
        end interface    
        
    contains
!---^^^^^^^^
 
            
        subroutine computeCrystalStructureFactors0( element, hkl, T,V ,Vg ,a0_in, a,b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the crystal structure factor from the element name and reflection at the given temperature and accelerating voltage
                    
            character(len=*),intent(in)                     ::      element
            integer,dimension(:,:),intent(in)               ::      hkl         !   (3,nReflections)
            real(kind=real64),intent(in)                    ::      T       
            real(kind=real64),intent(in)                    ::      V
            complex(kind=real64),dimension(:),intent(out)   ::      Vg          !   (1:nReflections)
            real(kind=real64),dimension(3),intent(in),optional          ::      a0_in       !   input lattice parameter, if known
            complex(kind=real64),dimension(:),intent(out),optional      ::      a,b

            character(len=LATTNAME_LEN)                     ::      latticename
!            real(kind=real64)                               ::      dwf
            real(kind=real64),dimension(3)                  ::      a0
            type(Lattice)                                   ::      latt
!            integer                                         ::      nMotif,nReflections,ii
            

            ! integer                                         ::      nS
            ! real(kind=real64)                               ::      smin,smax,ds        
            ! real(kind=real64),dimension(DT_COEFFS)          ::      ar,br,ai,bi             !   Doyle-Turner coefficients
            ! real(kind=real64),dimension(:),allocatable      ::      ff,ss           !   Scattering amplitudes
            ! real(kind=real64)                               ::      mser,msei
            ! complex(kind=real64),dimension(:,:),allocatable ::      aa,bb

        !---    find the lattice
            latticename =  getLatticeName(element) 
            if (getLatticeType(latticename) == LATTICE_CUSTOM) then
                print *,"Lib_CrystalStructureFactor::computeCrystalStructureFactors0 error - did not recognise lattice """//trim(latticename)//""""
                print *,"   only coded"
                call listAvailableLattices()
                stop
            end if            
            latt = Lattice_ctor(latticename)
            if (present(a0_in)) then
                a0 = a0_in
            else
                a0 = getLatticeConstant(element)      
            end if
            if (getLatticeType(latt) == LATTICE_HCP) call setCoverA(latt,a0(3)/a0(1))
            if (present(a)) then
                call computeCrystalStructureFactors1( element, latt, a0, hkl, T,V ,Vg , a,b)
            else
                call computeCrystalStructureFactors1( element, latt, a0, hkl, T,V ,Vg )
            end if
            return

            ! nMotif = getConventionalnMotif(latt)
            ! allocate(aa(DT_COEFFS,nMotif))
            ! allocate(bb(DT_COEFFS,nMotif))         
            ! print *,"Lib_CrystalStructureFactor::computeCrystalStructureFactors0 info -"
            ! call report(latt)

        ! !---    Debye Waller factor and lattice constant
        !     dwf = DebyeWallerFactor(element,T,latticename)
        !     a0 = getLatticeConstant(element)          
        !     write(*,fmt='(a,f12.6,a)')  " Lib_CrystalStructureFactor::computeCrystalStructureFactors0 info - Debye-Waller factor ",dwf," (A^2)"
        !     write(*,fmt='(a,3f12.6,a)') " Lib_CrystalStructureFactor::computeCrystalStructureFactors0 info - lattice parameters  ",a0," (A)"

        ! !---    find real DT coeffs
        !     smin = LIB_INPUTREAL_SMIN
        !     smax = LIB_INPUTREAL_SMAX
        !     nS = LIB_INPUTREAL_NS
        !     ds = (smax-smin)/nS
        !     allocate(ss(0:ns))
        !     allocate(ff(0:nS))
        !     do ii = 0,nS
        !         ss(ii) = smin + ii*ds
        !     end do
        !     call getScatteringAmplitude(element,ss,ff)
        !     ar = getrealDoyleTurner_a(element)
        !     br = getrealDoyleTurner_b(element)
        !     call fitDoyleTurner( ss,ff,ar,br, mser )
        !     deallocate(ss)
        !     deallocate(ff)
        !    ! print *,"Lib_CrystalStructureFactor::computeCrystalStructureFactors info - Doyle-Turner coeffs for real part of scattering "
        !     ! write(*,fmt='(a6,2a16)') "coeff"," a_r(ii)"," b_r(ii)"
        !     ! do ii = 1,DT_COEFFS
        !     !     write(*,fmt='(i6,2f16.8))') ii,ar(ii),br(ii)
        !     ! end do

            
        ! !---    find imag DT coeffs
        !     smin = LIB_INPUTTDS_SMIN
        !     smax = LIB_INPUTTDS_SMAX
        !     nS = LIB_INPUTTDS_NS
        !     ds = (smax-smin)/nS
        !     allocate(ss(0:ns))
        !     allocate(ff(0:nS))
        !     call getImagScatteringAmplitude(element,V,dwf,ss,ff)
        !     ai = getscatDoyleTurner_a(element,T)
        !     bi = getscatDoyleTurner_b(element,T)
        !     call fitDoyleTurner( ss,ff,ai,bi, msei )
        !     deallocate(ss)
        !     deallocate(ff)
        !     print *,"Lib_CrystalStructureFactor::computeCrystalStructureFactors0 info - Doyle-Turner coeffs "
        !     write(*,fmt='(a6,4a16)') "coeff"," a_r(ii)"," b_r(ii)"," a_i(ii)"," b_i(ii)"
        !     do ii = 1,DT_COEFFS
        !         write(*,fmt='(i6,4f16.8)') ii,ar(ii),br(ii),ai(ii),bi(ii)
        !     end do
        !     write(*,fmt='(3(a,f12.8))')" Lib_CrystalStructureFactor::computeCrystalStructureFactors0 info - mean square error ",mser," (real) ",msei," (imag)"
            
        ! !---    construct complex DT coeffs
            
        !     do ii = 1,DT_COEFFS
        !         aa(ii,:) = cmplx( ar(ii),ai(ii),kind=real64 )
        !         bb(ii,:) = cmplx( br(ii),bi(ii),kind=real64 )
        !     end do

        !     nReflections = size(hkl,dim=2)
        !     if (size(hkl,dim=1)==3) then
        !         !   hkl given as Miller indices
        !         do ii = 1,nReflections
        !             Vg(ii) = crystalStructureFactor1( a0, latticename, hkl(:,ii), aa,bb,dwf )
        !         end do
        !     else
        !         !   hkl given as Miller-Bravais indices
        !         do ii = 1,nReflections
        !             Vg(ii) = crystalStructureFactor1( a0, latticename, MillerBravaisToMiller_plane(hkl(:,ii)), aa,bb,dwf )
        !         end do
        !     end if


        !     if (present(a)) then
        !         do ii = 1,DT_COEFFS
        !             a(ii) = cmplx( ar(ii),ai(ii),kind=real64 )
        !         end do
        !     end if
        !     if (present(b)) then
        !         do ii = 1,DT_COEFFS
        !             b(ii) = cmplx( br(ii),bi(ii),kind=real64 )
        !         end do
        !     end if

            return
        end subroutine computeCrystalStructureFactors0


            
        subroutine computeCrystalStructureFactors1( element, latt, a0, hkl, T,V ,Vg , a,b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the crystal structure factor from the element name and lattice and reflection at the given temperature and accelerating voltage
                    
            character(len=*),intent(in)                     ::      element
            type(Lattice),intent(in)                        ::      latt
            real(kind=real64),dimension(3),intent(in)       ::      a0
            integer,dimension(:,:),intent(in)               ::      hkl         !   (3,nReflections)
            real(kind=real64),intent(in)                    ::      T       
            real(kind=real64),intent(in)                    ::      V
            complex(kind=real64),dimension(:),intent(out)   ::      Vg          !   (1:nReflections)
            complex(kind=real64),dimension(:),intent(out),optional      ::      a,b

            real(kind=real64)                               ::      dwf
            integer                                         ::      nMotif,nReflections,ii
            

            integer                                         ::      nS
            real(kind=real64)                               ::      smin,smax,ds        
            real(kind=real64),dimension(DT_COEFFS)          ::      ar,br,ai,bi             !   Doyle-Turner coefficients
            real(kind=real64),dimension(:),allocatable      ::      ff,ss           !   Scattering amplitudes
            real(kind=real64)                               ::      mser,msei
            complex(kind=real64),dimension(:,:),allocatable ::      aa,bb
 
            nMotif = getConventionalnMotif(latt)
            allocate(aa(DT_COEFFS,nMotif))
            allocate(bb(DT_COEFFS,nMotif))         
            print *,"Lib_CrystalStructureFactor::computeCrystalStructureFactors1 info -"
            call report(latt,verbose=.true.)

        !---    Debye Waller factor and lattice constant
            dwf = DebyeWallerFactor(element,T)
            write(*,fmt='(a,f12.6,a)')  " Lib_CrystalStructureFactor::computeCrystalStructureFactors1 info - Debye-Waller factor ",dwf," (A^2)"
            write(*,fmt='(a,3f12.6,a)') " Lib_CrystalStructureFactor::computeCrystalStructureFactors1 info - lattice parameters  ",a0," (A)"

        !---    find real DT coeffs
            smin = LIB_INPUTREAL_SMIN
            smax = LIB_INPUTREAL_SMAX
            nS = LIB_INPUTREAL_NS
            ds = (smax-smin)/nS
            allocate(ss(0:ns))
            allocate(ff(0:nS))
            do ii = 0,nS
                ss(ii) = smin + ii*ds
            end do
            call getScatteringAmplitude(element,ss,ff)
            ar = getrealDoyleTurner_a(element)
            br = getrealDoyleTurner_b(element)
            call fitDoyleTurner( ss,ff,ar,br, mser )
            deallocate(ss)
            deallocate(ff)
           ! print *,"Lib_CrystalStructureFactor::computeCrystalStructureFactors info - Doyle-Turner coeffs for real part of scattering "
            ! write(*,fmt='(a6,2a16)') "coeff"," a_r(ii)"," b_r(ii)"
            ! do ii = 1,DT_COEFFS
            !     write(*,fmt='(i6,2f16.8))') ii,ar(ii),br(ii)
            ! end do

            
        !---    find imag DT coeffs
            smin = LIB_INPUTTDS_SMIN
            smax = LIB_INPUTTDS_SMAX
            nS = LIB_INPUTTDS_NS
            ds = (smax-smin)/nS
            allocate(ss(0:ns))
            allocate(ff(0:nS))
            call getImagScatteringAmplitude(element,V,dwf,ss,ff)
            ai = getscatDoyleTurner_a(element,T)
            bi = getscatDoyleTurner_b(element,T)
            call fitDoyleTurner( ss,ff,ai,bi, msei )
            deallocate(ss)
            deallocate(ff)
            print *,"Lib_CrystalStructureFactor::computeCrystalStructureFactors1 info - Doyle-Turner coeffs "
            write(*,fmt='(a6,4a16)') "coeff"," a_r(ii)"," b_r(ii)"," a_i(ii)"," b_i(ii)"
            do ii = 1,DT_COEFFS
                write(*,fmt='(i6,4f16.8)') ii,ar(ii),br(ii),ai(ii),bi(ii)
            end do
            write(*,fmt='(3(a,f12.8))')" Lib_CrystalStructureFactor::computeCrystalStructureFactors1 info - mean square error ",mser," (real) ",msei," (imag)"
            
        !---    construct complex DT coeffs
            
            do ii = 1,DT_COEFFS
                aa(ii,:) = cmplx( ar(ii),ai(ii),kind=real64 )
                bb(ii,:) = cmplx( br(ii),bi(ii),kind=real64 )
            end do

            nReflections = size(hkl,dim=2)
            if (size(hkl,dim=1)==3) then
                !   hkl given as Miller indices
                do ii = 1,nReflections
                    Vg(ii) = crystalStructureFactor1( a0, latt, hkl(:,ii), aa,bb,dwf )
                end do
            else
                !   hkl given as Miller-Bravais indices
                do ii = 1,nReflections
                    Vg(ii) = crystalStructureFactor1( a0, latt, MillerBravaisToMiller_plane(hkl(:,ii)), aa,bb,dwf )
                end do
            end if


            if (present(a)) then
                do ii = 1,DT_COEFFS
                    a(ii) = cmplx( ar(ii),ai(ii),kind=real64 )
                end do
            end if
            if (present(b)) then
                do ii = 1,DT_COEFFS
                    b(ii) = cmplx( br(ii),bi(ii),kind=real64 )
                end do
            end if

            return
        end subroutine computeCrystalStructureFactors1        



        subroutine computeCrystalStructureFactors2( element, frac, latt, a0, hkl, T,V ,Vg ,a,b)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      find the crystal structure factor from the element name and lattice and elemental fraction and reflection at the given temperature and accelerating voltage
    !*      DOES NOT INCLUDE THE RANDOM OCCUPATION ADDITIONAL SCATTERING TERM 
    !*      See High Energy Electron Diffraction and Microscopy Peng, Dudarev, Whelan, OUP 2004 Chap 5.5
            character(len=*),dimension(:),intent(in)        ::      element
            real(kind=real64),dimension(:),intent(in)       ::      frac
            type(Lattice),intent(in)                        ::      latt
            real(kind=real64),dimension(3),intent(in)       ::      a0
            integer,dimension(:,:),intent(in)               ::      hkl         !   (3,nReflections)
            real(kind=real64),intent(in)                    ::      T       
            real(kind=real64),intent(in)                    ::      V
            complex(kind=real64),dimension(:),intent(out)   ::      Vg          !   (1:nReflections)
            complex(kind=real64),dimension(:),intent(out),optional      ::      a,b
            complex(kind=real64),dimension(size(hkl,dim=2)) ::      Vg_el
            integer         ::      ii

            if (present(a)) then
                a = 0
                b = 0
            end if
            Vg = 0
            do ii = 1,size(element)
                if ( (size(element) == 1).and.present(a) ) then
                    call computeCrystalStructureFactors1( element(ii), latt, a0, hkl, T,V ,Vg_el,a,b )
                else
                    call computeCrystalStructureFactors1( element(ii), latt, a0, hkl, T,V ,Vg_el )
                end if
                Vg(:) = Vg(:) + frac(ii) * Vg_el(:)
            end do
            Vg(:) = Vg(:) / sum( frac(1:size(element)) )
            
            return
        end subroutine computeCrystalStructureFactors2




        function crystalStructureFactor1( a0, latt, hkl, a,b,dwf ) result (Vg)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the crystal structure factor from the Doyle Turner Coefficients
    !*      on input g is in reduced coordinates eg [110] and converted to A^-1 
                    
            real(kind=real64),dimension(3),intent(in)       ::      a0          !   for cubic crystals only use a0(1). For hcp will also set c/a ratio.
            type(Lattice),intent(in)                        ::      latt
            complex(kind=real64),dimension(:,:),intent(in)  ::      a,b         !   (5,nMotif)
            real(kind=real64),intent(in)                    ::      dwf
            integer,dimension(3),intent(in)                 ::      hkl
            complex(kind=real64)                            ::      Vg
            
            !type(Lattice)                                   ::      latt
            real(kind=real64)                               ::      omega , aa
            real(kind=real64)                               ::      ss , gdotr
            real(kind=real64),dimension(3)                  ::      rr
            real(kind=real64),dimension(3)                  ::      gg
            real(kind=real64),dimension(3,3)                ::      bb
            complex(kind=real64)        ::      cDT,expgr 
            integer         ::      nMotif,ii
             
            
            ! if (getLatticeType(latticename) == LATTICE_CUSTOM) then
            !     print *,"Lib_CrystalStructureFactor::crystalStructureFactor1 error - did not recognise lattice """//trim(latticename)//""""
            !     print *,"   only coded"
            !     call listAvailableLattices()
            !     stop
            ! end if            
            ! latt = Lattice_ctor(latticename)
            aa = a0(1)                          !   characteristic length scale for the lattice
            ! if (getLatticeType(latt) == LATTICE_HCP) call setCoverA(latt,a0(3)/aa)

        !---    find the reciprocal lattice vectors. Note the reciprocal lattice stored in latt is for the crystal type, so need to set a length scale too.
            bb = getReciprocalLatticeVectors(latt) / aa
        
        !---    find the reflection considered
            gg = bb(:,1)*hkl(1) + bb(:,2)*hkl(2) + bb(:,3)*hkl(3) 

        !---    find the s parameter, so that V_g = V_g^0 Exp[ -B s^2 ] where B is the Debye-Waller factor   
        !       Acta Cryst. (1996). A52, 456-470 Eqn 1     
            ss = norm2(gg)/(4*PI)           
                
        !---    find the motif of the lattice, sum the crystal structure factor
        !       Acta Cryst. (1996). A52, 456-470 Eqn 13    
            nMotif = getConventionalnMotif(latt)
            Vg = 0
            do ii = 1,nMotif
                rr(:) = aa * getConventionalMotif( latt,ii,removeOffset=.true. )    !   remove offset to return motif point 1 at [000] not [¼¼¼]
                gdotr = dot_product( gg,rr )                
                expgr = cmplx( cos(gdotr),-sin(gdotr),kind=real64 )                           !   Exp( -i g.r )
                cDT = complexDoyleTurnerSum( a(:,ii),b(:,ii), dwf, ss )          
                Vg = Vg + expgr * cDT
            end do



        !---    find the volume per conventional cell and scale Vg
            omega = getOmega0(latt) * aa**3 

            Vg = Vg * (HBAR*HBAR/(2*ME))*(4*PI/(nMotif*Omega))

           
 
            return
        end function crystalStructureFactor1
        
         
                    
            
        pure function complexDoyleTurnerSum( a,b, dwf, s ) result( f )            
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
    !*      find the Doyle-Turner sum with complex coefficients
    !*          f = sum_i Re(a)_i Exp[ - (Re(a)_i + dwf) s^2 ] 
    !*            + i sum_i Im(a)_i Exp[ - (Im(a)_i + dwf/2 ) s^2 ] 
    !*      for complex input coefficients a_i,b_i 
    !*      and real Debye-Waller factor dwf
    !*      note that s = g/(4 pi)
    !*      see Acta Cryst. (1996). A52, 456-470 equation 13
    !*      ( note that the constant term and exp[ -i g.r ] terms are needed to compute Vg )

            complex(kind=real64),dimension(:),intent(in)    ::      a,b     !   complex Doyle-Turner coefficient
            real(kind=real64),intent(in)                    ::      dwf     !   Debye-Waller factor
            real(kind=real64),intent(in)                    ::      s       !   scattering angle s = g/4pi = sin(theta)/lambda
            
            complex(kind=real64)                            ::      f
            
            real(kind=real64)               ::      rDTc,iDTc
            
            rDTc = DoyleTurnerSum( s, real(a) ,real(b) + dwf    )
            iDTc = DoyleTurnerSum( s, aimag(a),aimag(b) + dwf/2 )
            f = cmplx( rDTc,iDTc,kind=real64 )
        
            return
        end function complexDoyleTurnerSum
                    
         
        
    end module Lib_CrystalStructureFactor     
    
    
