 



add_executable(hbtem-cp ${CMAKE_SOURCE_DIR}/src/run_SimulatedTEM.F90)
target_link_libraries(hbtem-cp ManyBeam_Utils CSF_Utils PCA_Png CommandLine)

#Install the exec to bin folder of install path
install(TARGETS hbtem-cp
    RUNTIME
    DESTINATION bin)


