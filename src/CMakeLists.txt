include_directories(${libadjoint_SOURCE_DIR}/include)

# Here comes a hack: The C headers are adapted for inclusion in Fortran files.
# TODO remove this hack
execute_process(COMMAND sh -c "sed -e 's@/\\*@!@' -e 's@\\*/@@' -e 's@ADJ_CONSTANTS_H@ADJ_CONSTANTS_F_H@' -e '/adj_scalar /d' ${libadjoint_SOURCE_DIR}/include/libadjoint/adj_constants.h > ${libadjoint_SOURCE_DIR}/include/libadjoint/adj_constants_f.h")
execute_process(COMMAND sh -c "grep '^#' ${libadjoint_SOURCE_DIR}/include/libadjoint/adj_error_handling.h | grep -v '^#include' | grep -v 'CHKMALLOC' > ${libadjoint_SOURCE_DIR}/include/libadjoint/adj_error_handling_f.h")

file(GLOB libadjoint_SRCS
     "${CMAKE_CURRENT_SOURCE_DIR}/*.c"
     "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
     "${CMAKE_CURRENT_SOURCE_DIR}/*.F90"
     )
file(GLOB libadjoint_HEADERS
     "${libadjoint_SOURCE_DIR}/include/libadjoint/*.h"
     "${libadjoint_SOURCE_DIR}/include/libadjoint/*.mod"
     )

add_library(adjoint ${libadjoint_SRCS})

set_target_properties(adjoint PROPERTIES
  PUBLIC_HEADER "${libadjoint_HEADERS}"
  VERSION ${libadjoint_VERSION}
  SOVERSION ${libadjoint_MAJOR_VERSION}
  )

# Installation of the program
install(TARGETS adjoint
  # We don't actually need the RUNTIME DESTINATION here since we're only
  # installing a library.
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib # .so
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib # .a
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" COMPONENT dev
  #PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/libadjoint" COMPONENT dev
  )

