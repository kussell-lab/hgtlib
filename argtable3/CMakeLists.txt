INCLUDE_DIRECTORIES(.)
file(GLOB argtable3_SRC
    "*.h"
    "*.c"
)
add_library(argtable3 ${argtable3_SRC})

file(GLOB argtable3_HEAD "*.h")

install(TARGETS argtable3 DESTINATION lib)
install(FILES ${argtable3_HEAD} DESTINATION include)