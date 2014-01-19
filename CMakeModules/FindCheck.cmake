# -- Try to find check library
# Once done this will define
# 	CHECK_FOUND -- system has check library
# 	CHECK_INCLUDE_DIRS
# 	CHECK_LIBRARIES

find_package(PkgConfig)
pkg_check_modules(PC_CHECK QUIET check)
set(CHECK_DEFINITIONS ${PC_CHECK_CFLAGS_OTHER})

find_path(CHECK_INCLUDE_DIR check.h
          HINTS ${PC_CHECK_INCLUDEDIR} ${PC_CHECK_INCLUDE_DIRS} 
          PATHS ENV CHECK_INC
          PATH_SUFFIXES check )

find_library(CHECK_LIBRARY NAMES check libcheck
             HINTS ${PC_CHECK_LIBDIR} ${PC_CHECK_LIBRARY_DIRS} 
             PATHS ENV CHECK_BINARY )

set(CHECK_LIBRARIES ${CHECK_LIBRARY} )
set(CHECK_INCLUDE_DIRS ${CHECK_INCLUDE_DIR} )
set(CHECK_CFLAGS_OTHER ${PC_CHECK_CFLAGS_OTHER})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CHECK_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Check  DEFAULT_MSG
                                  CHECK_LIBRARY CHECK_INCLUDE_DIR)

mark_as_advanced(CHECK_INCLUDE_DIR CHECK_LIBRARY )

IF(CHECK_LIBRARIES)
  IF(CHECK_INCLUDE_DIR)

    SET(CHECK_FOUND 1)
    
    MESSAGE(STATUS "Using Check from ${CHECK_LIBRARY}")

  ENDIF(CHECK_INCLUDE_DIR)
ENDIF(CHECK_LIBRARIES)