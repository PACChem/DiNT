####################################################
#                                                  #
# Compiles DiNT (Direct Nonadiabatic Trajectories) #
#   for Linux                                      #
#                                                  #
####################################################
project(DiNT Fortran)
cmake_minimum_required(VERSION 3.0)

# Define potential to compile
set(POT
#"exp6"
#"h3o"
#"tbplusexp6-2"
#"h2oM"
"ch4hALL"
#"ch4ohALL"

##TINKER PES
#NOTE: obj_tinker contains link to mm3 PES
#NOTE: obj_tinker2 contains link to amoeba09 PES
#NOTE: obj_tinker3 contains link to oplsaa PES
#rohMP
#roohT
#alkT
)

set(CMAKE_LIBDIR "${PROJECT_SOURCE_DIR}/lib")
set(CMAKE_OBJDIR "${PROJECT_SOURCE_DIR}/src/obj_opt")
#set(CMAKE_INSTALL_DIR "${PROJECT_SOURCE_DIR}/exe")
set(SPRNG_LIBDIR "${CMAKE_SOURCE_DIR}/sprng/lib")

option(USE_TINKER "Use Tinker software" OFF)

if(USE_TINKER)
  list(APPEND EXTRA_INCLUDES "/home/ajasper/KTP/dint/tinker/tinker/source")
endif(USE_TINKER)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

add_subdirectory("src")
#add_subdirectory("pot")
#add_subdirectory("sprng")
list(APPEND EXTRA_INCLUDES "${CMAKE_SOURCE_DIR}/src")
list(APPEND EXTRA_INCLUDES "${CMAKE_SOURCE_DIR}/pot")
list(APPEND EXTRA_INCLUDES "${CMAKE_SOURCE_DIR}/sprng/SRC")

# Include libraries
#list(APPEND EXTRA_LIBS "${CMAKE_SOURCE_DIR}/pot/libch4libs.a")
list(APPEND EXTRA_LIBS "${SPRNG_LIBDIR}/libcmrg.a")
list(APPEND EXTRA_LIBS "${SPRNG_LIBDIR}/liblcg64.a")
list(APPEND EXTRA_LIBS "${SPRNG_LIBDIR}/liblcg.a")
list(APPEND EXTRA_LIBS "${SPRNG_LIBDIR}/liblfg.a")
list(APPEND EXTRA_LIBS "${SPRNG_LIBDIR}/libsprngtest.a")
#list(APPEND EXTRA_LIBS "${CMAKE_BINARY_DIR}/src/CMakeFiles/objs.dir/")
#list(APPEND EXTRA_LIBS "${CMAKE_BINARY_DIR}/pot/CMakeFiles/libch4libs.a")

add_executable(dint-${POT}.x $<TARGET_OBJECTS:objs>
    ${CMAKE_SOURCE_DIR}/src/)

target_link_libraries(dint-${POT}.x PUBLIC
    ${EXTRA_LIBS} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

target_include_directories(dint-${POT}.x PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          ${EXTRA_INCLUDES}
                          )

install(TARGETS dint-${POT}.x DESTINATION exe)
