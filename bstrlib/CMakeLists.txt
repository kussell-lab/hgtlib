INCLUDE_DIRECTORIES(.)
file(GLOB bstrlib_SRC
    "*.h"
    "*.c"
)
add_library(bstrlib ${bstrlib_SRC})

file(GLOB bstrlib_HEAD "*.h")

install(TARGETS bstrlib DESTINATION lib)
install(FILES ${bstrlib_HEAD} DESTINATION include)