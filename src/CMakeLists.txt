find_package(GSL REQUIRED)
find_package(ArgTable REQUIRED)
include_directories(. ${GSL_INCLUDE_DIR} ${ARGTABLE_INCLUDE_DIRS} ../bstrlib)
file(GLOB hgt_SRC
    "*.h"
    "*.c"
)
add_library(hgt ${hgt_SRC})
