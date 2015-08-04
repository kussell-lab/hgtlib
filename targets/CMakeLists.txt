find_package(MPI REQUIRED)
find_package(GSL REQUIRED)

set(COMM_LIBS ${GSL_LIBRARIES} hgt bstrlib argtable3)
include_directories(${GSL_INCLUDE_DIR} ${MPI_INCLUDE_PATH} . ../src ../bstrlib)

add_executable(hgt_mpi_moran_const hgt_mpi_moran_const.c)
target_link_libraries(hgt_mpi_moran_const ${MPI_LIBRARIES} ${COMM_LIBS})
install(TARGETS hgt_mpi_moran_const DESTINATION bin)

add_executable(hgt_forward_coal_back hgt_forward_coal_back.c hgt_forward_coal_back.h)
target_link_libraries(hgt_forward_coal_back ${COMM_LIBS})
install(TARGETS hgt_forward_coal_back DESTINATION bin)

add_executable(hgt_simu hgt_simu.c)
target_link_libraries(hgt_simu ${COMM_LIBS})
install(TARGETS hgt_simu DESTINATION bin)

add_custom_target(hgt_simu_debug COMMAND "${CMAKE_CURRENT_BINARY_DIR}/Debug/hgt_simu" -C "${CMAKE_CURRENT_LIST_DIR}/examples/run_hgt_simu_example1.cfg" DEPENDS hgt_simu)

if(WIN32)
	message("Using mpiexe for running examples.")
	set(MPIRUN "mpiexec")
else(WIN32)
	message("Using mpirun for running examples.")
	set(MPIRUN "mpirun")
endif(WIN32)

add_test("run_mpi_example1" "${MPIRUN}" "-n" "2" "${CMAKE_CURRENT_BINARY_DIR}/hgt_mpi_moran_const" "-C" "${CMAKE_CURRENT_LIST_DIR}/examples/run_mpi_example1.cfg")
add_test("run_mpi_example2" "${MPIRUN}" "-n" "2" "${CMAKE_CURRENT_BINARY_DIR}/hgt_mpi_moran_const" "-C" "${CMAKE_CURRENT_LIST_DIR}/examples/run_mpi_example2.cfg")
add_test("run_mpi_example3" "${MPIRUN}" "-n" "2" "${CMAKE_CURRENT_BINARY_DIR}/hgt_mpi_moran_const" "-C" "${CMAKE_CURRENT_LIST_DIR}/examples/run_mpi_example3.cfg")