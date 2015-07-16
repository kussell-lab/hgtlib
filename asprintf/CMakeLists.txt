INCLUDE_DIRECTORIES(.)
file(GLOB asprintf_SRC
    "*.h"
    "*.c"
)
add_library(asprintf ${asprintf_SRC})

file(GLOB asprintf_HEAD "*.h")

install(TARGETS asprintf DESTINATION lib)
install(FILES ${asprintf_HEAD} DESTINATION include)