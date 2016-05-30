find_package(GSL REQUIRED)
set(ARGTABLE_INCLUDE_DIRS "../argtable3")
set(BSTRLIB_INCLUDE_DIRS "../bstrlib")
IF(NOT CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_C_FLAGS "-O2")
ENDIF(NOT CMAKE_BUILD_TYPE MATCHES Debug)
include_directories(. ${GSL_INCLUDE_DIR} ${ARGTABLE_INCLUDE_DIRS} ${BSTRLIB_INCLUDE_DIRS})
file(GLOB hgt_SRC
    "*.h"
    "*.c"
)
add_library(hgt ${hgt_SRC})
file(GLOB hgt_HEAD "*.h")
install(TARGETS hgt DESTINATION lib)
install(FILES ${hgt_HEAD} DESTINATION include)
