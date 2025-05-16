!   A module which defines radio-button
!   on-off flag states for a simulation
!   and can report them.
!   A simulation may have flags such as "linear approximation" or "data loaded successfully"
!   or anything really. The design decision made here is to outsource
!   the definition of the data to a header file "dashboard.h"
!   which will provide domain-specific names for the flags.
!   This is portable, and reasonably safe- you won't be able to compile the
!   "asteroids_game" code with a set of flags appropriate for "dislocation_loop_energy"
!
!   There is currently a maximum of 32 flags which can be set, so they store as a single 32 bit int.
!   A good program state is shared by all processors. This I can't guarantee from here - so beware!
!   Note that I also can't guarantee from here that just because you claim to have "data loaded successfully"
!   that this is actually true. So some effort is expected on the part of the user to ensure such a flag is
!   only set when it is valid to do so.
!   I am making no attempt here to have any dynamic testing of states
!
!   The program state is a public static variable
!       integer(kind=PROGRAM_STATE),public      ::      state
!
!   dashboard.h MUST contain the following datafields
!       integer,private,parameter   ::      NUMBER_FLAGS
!       character(len=FLAG_NAME_LEN),dimension(NUMBER_FLAGS),private,parameter  ::      FLAG_NAME
!   and SHOULD contain mnemonics for bits defined as the flags eg
!       integer(kind=PROGRAM_STATE),public,parameter    ::      FLAG_ONE = int(b'00001',kind=PROGRAM_STATE)     !   note BOZ constant
!       ...
!       
!
!   Daniel Mason
!   (c) UKAEA Feb 2025
!   
!   version history
!       0.0.1       Feb 2025        First working version
!

    module Lib_Dashboard
!---^^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        implicit none
        private

                
        include "dashboard.h"
        

    !---    this is the current program state - a public static variable
        integer(kind=PROGRAM_STATE),public      ::      state = 0_PROGRAM_STATE            

    !---    accessor functions
        public          ::      setState
        public          ::      unsetState
        public          ::      isState
    !---    simple screen dump of the program state
        public          ::      reportDashboard


        interface setState
            module procedure setState0
            module procedure setState1
        end interface
        interface unsetState
            module procedure unsetState0
        end interface
        interface isState
            module procedure isState0
        end interface
        interface reportDashboard
            module procedure reportState0
        end interface   




    contains
!---^^^^^^^^
    
    
        subroutine setState0(flag)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set the state flag to true
    !*      this call can, in fact set multiple flags at once, but its probably best 
    !*      to do one at a time for clarity
            integer(kind=PROGRAM_STATE),intent(in)      ::      flag
            state = ior( state,flag )
            return
        end subroutine setState0   
    
        subroutine setState1(flag,onoff)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set the state flag to on or off
    !*      this call can, in fact set multiple flags at once, but its probably best 
    !*      to do one at a time for clarity
            integer(kind=PROGRAM_STATE),intent(in)      ::      flag
            logical,intent(in)                          ::      onoff
            if (onoff) then
                call setState(flag)
            else
                call unsetState(flag)
            end if
            return
        end subroutine setState1   

        subroutine unsetState0(flag)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      set the state flag to false
    !*      this call can, in fact set multiple flags at once, but its probably best 
    !*      to do one at a time for clarity
            integer(kind=PROGRAM_STATE),intent(in)      ::      flag
            state = ior( state,flag ) - flag
            return
        end subroutine unsetState0   

        pure logical function isState0(flag)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if the state flag is true
    !*      ( or all the flags are true )
            integer(kind=PROGRAM_STATE),intent(in)      ::      flag
            isState0 = iand( state,flag ) == flag
            return
        end function isState0 

        subroutine reportState0(u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple screen dump to unit u ( default screen )
    !*      with offset o (default 0)
            integer,intent(in),optional     ::      u
            integer,intent(in),optional     ::      o
            integer     ::      uu,oo
            integer     ::      bit
            
            uu = OUTPUT_UNIT ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o 
            do bit = 0,NUMBER_FLAGS-1
                write(unit=uu,fmt='(a)') repeat(" ",oo)//reportFlagbit( bit,screen=(uu==OUTPUT_UNIT) )
            end do 
            return
        end subroutine reportState0

        pure function reportFlagbit( bit,screen ) result(text)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                          ::      bit
            logical,intent(in)                          ::      screen
            character(len=FLAG_NAME_LEN+1+5+9)          ::      text            
            if (screen) then
                if (btest(state,bit)) then
                    text = FLAG_NAME(bit+1)//" "//colour("92","true")
                else
                    text = FLAG_NAME(bit+1)//" "//colour("31","false")
                end if
            else
                if (btest(state,bit)) then
                    text = FLAG_NAME(bit+1)//" true"
                else
                    text = FLAG_NAME(bit+1)//" false"
                end if
            end if
            return
        end function reportFlagbit

        pure function colour(COLOUR_CHARCODE,text) result(ctext)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
    !*      produce a text string "text" decorated with terminal markup so it appears in the requested colour.
            character(len=2),intent(in)         ::      COLOUR_CHARCODE
            character(len=*),intent(in)         ::      text
            character(len=len_trim(text)+9)     ::      ctext
            
            ctext = achar(27)//"["//COLOUR_CHARCODE//"m"//trim(text)//achar(27)//"[0m"
            
            return
        end function colour

    end module Lib_Dashboard