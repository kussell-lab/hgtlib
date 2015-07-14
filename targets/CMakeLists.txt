find_package(MPI REQUIRED)
find_package(GSL REQUIRED)
find_package(ArgTable REQUIRED)

include_directories(${GSL_INCLUDE_DIR} ${ARGTABLE_INCLUDE_DIRS} ${MPI_INCLUDE_PATH} . ../src ../bstrlib)
set(LIBS ${LIBS} ${MPI_LIBRARIES} ${GSL_LIBRARIES} ${ARGTABLE_LIBRARIES} hgt bstrlib)
MESSAGE(STATUS ${MPI_INCLUDE_PATH})

set(CMAKE_C_COMPILER cc)
add_executable(hgt_mpi_moran_step hgt_mpi_moran_step.c)
target_link_libraries(hgt_mpi_moran_step ${LIBS})
install(TARGETS hgt_mpi_moran_step DESTINATION bin)

set(CMAKE_C_COMPILER cc)
add_executable(hgt_mpi_moran_const hgt_mpi_moran_const.c)
target_link_libraries(hgt_mpi_moran_const ${LIBS})
install(TARGETS hgt_mpi_moran_const DESTINATION bin)

set(CMAKE_C_COMPILER cc)
add_executable(hgt_mpi_moran_expon hgt_mpi_moran_expon.c)
target_link_libraries(hgt_mpi_moran_expon ${LIBS})
install(TARGETS hgt_mpi_moran_expon DESTINATION bin)

