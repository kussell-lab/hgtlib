INCLUDE_DIRECTORIES(.)
file(GLOB bstrlib_SRC
    "*.h"
    "*.c"
)
add_library(bstrlib ${bstrlib_SRC})