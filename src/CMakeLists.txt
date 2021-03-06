include_directories(${libadjoint_SOURCE_DIR}/include)

file(GLOB libadjoint_SRCS
     "${CMAKE_CURRENT_SOURCE_DIR}/*.c"
     "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
     )
file(GLOB libadjoint_HEADERS
     "${libadjoint_SOURCE_DIR}/include/libadjoint/*.h"
     )

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

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
find_package(PETSc 3.3)
if (PETSC_FOUND)
  find_package(SLEPc 3.3)
  if (SLEPC_FOUND)
    include_directories(SYSTEM ${PETSC_INCLUDE_DIRS})
    include_directories(SYSTEM ${SLEPC_INCLUDE_DIRS})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_SLEPC")
    target_link_libraries(adjoint ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})
    target_link_libraries(adjoint-static ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})
  endif()
endif()


# Installation of the program
install(TARGETS adjoint adjoint-static
  # We don't actually need the RUNTIME DESTINATION here since we're only
  # installing a library.
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib # .so
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib # .a
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" COMPONENT dev
  )

# -----------------------------------------------------------------------------
find_package(PythonInterp REQUIRED)
# -----------------------------------------------------------------------------
