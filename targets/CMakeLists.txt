find_package(MPI REQUIRED)
find_package(GSL REQUIRED)
find_package(ArgTable REQUIRED)
find_package(Gperftools REQUIRED)
set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${MPI_COMPILER_FLAGS})
set(CMAKE_LINK_LIBRARY_FLAG ${CMAKE_LINK_LIBRARY_FLAG} ${MPI_LINK_FLAGS} ${Gperftools_LINK_FLAGS})
include_directories(${GSL_INCLUDE_DIR} ${ARGTABLE_INCLUDE_DIRS} ${MPI_INCLUDE_PATH} ${GPERFTOOLS_INCLUDE_DIR} . ../src ../bstrlib)
set(LIBS ${LIBS} ${MPI_LIBRARIES} ${GSL_LIBRARIES} ${ARGTABLE_LIBRARIES} ${GPERFTOOLS_LIBRARIES} hgt bstrlib)

set(CMAKE_C_COMPILER mpicc)
add_executable(hgt_mpi_moran hgt_mpi_moran.c)
target_link_libraries(hgt_mpi_moran ${LIBS})

set(CMAKE_C_COMPILER mpicc)
add_executable(hgt_mpi_moran_step hgt_mpi_moran_step.c)
target_link_libraries(hgt_mpi_moran_step ${LIBS})

install(TARGETS hgt_mpi_moran_step DESTINATION bin)