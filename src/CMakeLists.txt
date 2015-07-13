find_package(GSL REQUIRED)
find_package(ArgTable REQUIRED)
set(CMAKE_C_FLAGS "-pthread -pg -g")
include_directories(. ${GSL_INCLUDE_DIR} ${ARGTABLE_INCLUDE_DIRS} ../bstrlib)
file(GLOB hgt_SRC
    "*.h"
    "*.c"
)
add_library(hgt ${hgt_SRC} hgt_linkage.h hgt_linkage.c)

file(GLOB hgt_HEAD "*.h")
install(TARGETS hgt DESTINATION lib)
install(FILES ${hgt_HEAD} DESTINATION include)