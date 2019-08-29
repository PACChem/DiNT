####################################################
#                                                  #
# Compiles DiNT (Direct Nonadiabatic Trajectories) #
#   for Linux                                      #
#                                                  #
####################################################
project(DiNT Fortran)
cmake_minimum_required(VERSION 1.1)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_library(ljlibs  #Lennard-Jones
    ${PROJECT_SOURCE_DIR}/src/LJ.f)

#add_library(tbpluslibs  #
#        ${PROJECT_SOURCE_DIR}/pot/tbplusexp6.f
#        ${PROJECT_SOURCE_DIR}/pot/tbplusexp6dt.f
#        ${PROJECT_SOURCE_DIR}/pot/tbplusexp6all.f
#        ${PROJECT_SOURCE_DIR}/pot/tbplusexp6-2.f)

#add_library(ch4libs  #Methane
#        ${PROJECT_SOURCE_DIR}/pot/ch4hALL.f
#        ${PROJECT_SOURCE_DIR}/pot/ch4oALL.f
#        ${PROJECT_SOURCE_DIR}/pot/ch4ohALL.f
#        ${PROJECT_SOURCE_DIR}/pot/ch4o.f
#        ${PROJECT_SOURCE_DIR}/pot/ch4oh.f)

#add_library(ch3flibs  #Fluoromethane
#        ${PROJECT_SOURCE_DIR}/pot/ch3f-rv2.f
#        ${PROJECT_SOURCE_DIR}/pot/ch3f-rv3aj2.f
#        ${PROJECT_SOURCE_DIR}/pot/ch3f-rv3aj2st.f
#        ${PROJECT_SOURCE_DIR}/pot/ch3f-rv3aj.f
#        ${PROJECT_SOURCE_DIR}/pot/ch3f-rv3.f
#        ${PROJECT_SOURCE_DIR}/pot/ch3f-rv4.f
#        ${PROJECT_SOURCE_DIR}/pot/ch3f-rv.f
#        ${PROJECT_SOURCE_DIR}/pot/ch3f-tmp.f)

#add_library(ho2libs  #HO2
#        ${PROJECT_SOURCE_DIR}/pot/ho2m.f
#        ${PROJECT_SOURCE_DIR}/pot/ho2h2oLS2.f
#        ${PROJECT_SOURCE_DIR}/pot/ho2heLS2.f
#        ${PROJECT_SOURCE_DIR}/pot/ho2heLS.f
#        ${PROJECT_SOURCE_DIR}/pot/ho2n2LS.f
#        ${PROJECT_SOURCE_DIR}/pot/ho2n2LS3C.f
#        ${PROJECT_SOURCE_DIR}/pot/ho2n2LSB.f
#        ${PROJECT_SOURCE_DIR}/pot/ho2n2LSC.f
#        ${PROJECT_SOURCE_DIR}/pot/ho2n2LSD.f
#        ${PROJECT_SOURCE_DIR}/pot/ho2n2LSE.f
#        ${PROJECT_SOURCE_DIR}/pot/j.f)

message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

add_executable(dint-${POT}.x
    ${PROJECT_SOURCE_DIR}/src/dint.F)

target_link_libraries(dint-${POT}.x 
        ljlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

install(TARGETS dint-${POT}.x DESTINATION bin)
