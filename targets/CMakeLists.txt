find_package(MPI REQUIRED)
find_package(GSL REQUIRED)

include_directories(${GSL_INCLUDE_DIR} ${MPI_INCLUDE_PATH} . ../src ../bstrlib)
set(LIBS ${LIBS} ${MPI_LIBRARIES} ${GSL_LIBRARIES} hgt bstrlib argtable3)
MESSAGE(STATUS ${MPI_INCLUDE_PATH})

set(CMAKE_C_COMPILER cc)
add_executable(hgt_mpi_moran_const hgt_mpi_moran_const.c)
target_link_libraries(hgt_mpi_moran_const ${LIBS})
install(TARGETS hgt_mpi_moran_const DESTINATION bin)

add_test("run_mpi_example1" "mpirun" "-n" "2" "${CMAKE_CURRENT_BINARY_DIR}/hgt_mpi_moran_const" "-C" "${CMAKE_CURRENT_LIST_DIR}/examples/run_mpi_example1.cfg")
add_test("run_mpi_example2" "mpirun" "-n" "2" "${CMAKE_CURRENT_BINARY_DIR}/hgt_mpi_moran_const" "-C" "${CMAKE_CURRENT_LIST_DIR}/examples/run_mpi_example2.cfg")
add_test("run_mpi_example3" "mpirun" "-n" "2" "${CMAKE_CURRENT_BINARY_DIR}/hgt_mpi_moran_const" "-C" "${CMAKE_CURRENT_LIST_DIR}/examples/run_mpi_example3.cfg")