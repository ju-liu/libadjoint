include_directories(${libadjoint_SOURCE_DIR}/include)

# Here comes a hack: The C headers are adapted for inclusion in Fortran files.
# TODO remove this hack
execute_process(COMMAND sh -c "sed -e 's@/\\*@!@' -e 's@\\*/@@' -e 's@ADJ_CONSTANTS_H@ADJ_CONSTANTS_F_H@' -e '/adj_scalar /d' ${libadjoint_SOURCE_DIR}/include/libadjoint/adj_constants.h > ${libadjoint_SOURCE_DIR}/include/libadjoint/adj_constants_f.h")
execute_process(COMMAND sh -c "grep '^#' ${libadjoint_SOURCE_DIR}/include/libadjoint/adj_error_handling.h | grep -v '^#include' | grep -v 'CHKMALLOC' > ${libadjoint_SOURCE_DIR}/include/libadjoint/adj_error_handling_f.h")

if (CMAKE_Fortran_COMPILER_WORKS)
	file(GLOB libadjoint_SRCS
	     "${CMAKE_CURRENT_SOURCE_DIR}/*.c"
	     "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
	     "${CMAKE_CURRENT_SOURCE_DIR}/*.F90"
	     )
	file(GLOB libadjoint_HEADERS
	     "${libadjoint_SOURCE_DIR}/include/libadjoint/*.h"
	     "${libadjoint_SOURCE_DIR}/include/libadjoint/*.mod"
	     )
else()
	file(GLOB libadjoint_SRCS
	     "${CMAKE_CURRENT_SOURCE_DIR}/*.c"
	     "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
	     )
	file(GLOB libadjoint_HEADERS
	     "${libadjoint_SOURCE_DIR}/include/libadjoint/*.h"
	     )
endif()

# Build both shared and static libraries here. Do not share the object files
# (which would be possible using CMake's OBJECT keyword) since objects intended
# for linking into a shared library need to be compiled with -fPIC.
add_library(adjoint SHARED ${libadjoint_SRCS})

set_target_properties(adjoint PROPERTIES
  PUBLIC_HEADER "${libadjoint_HEADERS}"
  VERSION ${libadjoint_VERSION}
  SOVERSION ${libadjoint_MAJOR_VERSION}
  )

add_library(adjoint-static STATIC ${libadjoint_SRCS})
set_target_properties(adjoint-static PROPERTIES
  OUTPUT_NAME adjoint
  )

# Installation of the program
install(TARGETS adjoint adjoint-static
  # We don't actually need the RUNTIME DESTINATION here since we're only
  # installing a library.
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib # .so
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib # .a
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" COMPONENT dev
  #PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/libadjoint" COMPONENT dev
  )

# -----------------------------------------------------------------------------
# Now some hardcore hacks for the Python bindings.
find_package(PythonInterp REQUIRED)

get_filename_component(LIB_REALPATH ${CMAKE_BINARY_DIR}/src/libadjoint${CMAKE_SHARED_LIBRARY_SUFFIX}
                      REALPATH
                      )
add_custom_command(TARGET adjoint
                  POST_BUILD
                  COMMAND cp ${libadjoint_SOURCE_DIR}/python/libadjoint/clibadjoint.py ${CMAKE_BINARY_DIR}/clibadjoint.py
                  COMMAND cp ${libadjoint_SOURCE_DIR}/python/libadjoint/clibadjoint_constants.py ${CMAKE_BINARY_DIR}/clibadjoint_constants.py
                  COMMENT "Generating clibadjoint.py..."
                  )
# Some extra adaptations for OS X.
IF(APPLE)
   add_custom_command(TARGET adjoint
                      POST_BUILD
                      COMMAND sed -i.bak \"s/._pack_ = 4/._pack_ = 8/\" ${CMAKE_BINARY_DIR}/clibadjoint.py
                      COMMENT "Adaptations for OS X..."
                      )
ENDIF()
add_custom_command(TARGET adjoint
                  POST_BUILD
                  COMMAND chmod a-x ${CMAKE_BINARY_DIR}/clibadjoint.py
                  COMMENT "Changing permissions on clibadjoint.py..."
                  )

# Make sure that the installed clibadjoint.py has the
# location of the shared library hardcoded in it.
add_custom_command(TARGET adjoint
                  POST_BUILD
                  COMMAND sed -i.bak -e 's@LIBADJOINT_LIBRARY_FILE@"${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/libadjoint${CMAKE_SHARED_LIBRARY_SUFFIX}"@' ${CMAKE_BINARY_DIR}/clibadjoint.py
                  COMMENT "Adaptations for clibadjoint.py..."
                  )

# Make sure to install the generated files
install(FILES ${CMAKE_BINARY_DIR}/clibadjoint.py
              ${CMAKE_BINARY_DIR}/clibadjoint_constants.py
              DESTINATION lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/libadjoint/
        )
# -----------------------------------------------------------------------------
