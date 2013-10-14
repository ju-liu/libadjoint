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
install(TARGETS adjoint
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

# 1.) Generate pylibadjoint.h and mess around with it.
add_custom_command(TARGET adjoint
                   POST_BUILD
                   COMMAND ${CMAKE_C_COMPILER} -E -U__BLOCKS__ -DPYTHON_BINDINGS ${libadjoint_SOURCE_DIR}/include/libadjoint/libadjoint.h > ${CMAKE_BINARY_DIR}/pylibadjoint.h
                   COMMENT "Generating pylibadjoint.h..."
                   )
                 add_custom_command(TARGET adjoint
                                    POST_BUILD
                                    COMMAND sed -i s/__builtin___stpncpy_chk/__builtin___strncpy_chk/ ${CMAKE_BINARY_DIR}/pylibadjoint.h
                                    COMMENT "Adapting pylibadjoint.h..."
                                    )
                 add_custom_command(TARGET adjoint
                                    POST_BUILD
                                    COMMAND ${libadjoint_SOURCE_DIR}/python/ctypeslib/scripts/h2xml.py -q -I. ${CMAKE_BINARY_DIR}/pylibadjoint.h -o ${CMAKE_BINARY_DIR}/libadjoint.xml
                                    COMMENT "Generating libadjoint.xml..."
                                    )
                 add_custom_command(TARGET adjoint
                                    POST_BUILD
                                    COMMAND rm -f ${CMAKE_BINARY_DIR}/include/libadjoint/pylibadjoint.h
                                    COMMENT "Removing pylibadjoint.h..."
                                    )
                 get_filename_component(LIB_REALPATH ${CMAKE_BINARY_DIR}/src/libadjoint.so
                                        REALPATH
                                        )
                 add_custom_command(TARGET adjoint
                                    POST_BUILD
                                    COMMAND ${libadjoint_SOURCE_DIR}/python/ctypeslib/scripts/xml2py.py -r '^adj.*' -l ${LIB_REALPATH} ${CMAKE_BINARY_DIR}/libadjoint.xml -o ${CMAKE_BINARY_DIR}/clibadjoint.py
                                    COMMENT "Generating clibadjoint.py..."
                                    )
                 # Some extra adaptations for OS X.
                 IF(APPLE)
                     add_custom_command(TARGET adjoint
                                        POST_BUILD
                                        COMMAND sed -i \"s/._pack_ = 4/._pack_ = 8/\" ${CMAKE_BINARY_DIR}/clibadjoint.py
                                        COMMENT "Adaptations for OS X..."
                                        )
                 ENDIF()
                 add_custom_command(TARGET adjoint
                                    POST_BUILD
                                    COMMAND rm -f ${CMAKE_BINDARY_DIR}/libadjoint.xml
                                    COMMENT "Removing libadjoint.xml..."
                                    )
                 add_custom_command(TARGET adjoint
                                    POST_BUILD
                                    COMMAND chmod a-x ${CMAKE_BINARY_DIR}/clibadjoint.py
                                    COMMENT "Changing permissions on clibadjoint.py..."
                                    )

                 # 2.) Make sure that the installed clibadjoint.py has the
                 #     location of the shared library hardcoded in it.
                 add_custom_command(TARGET adjoint
                                    POST_BUILD
                                    COMMAND sed -i -e "s@CDLL\\('${libadjoint_SOURCE_DIR}/lib/libadjoint.so'\\)@CDLL\\('${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/libadjoint.so'\\)@" -e "s@CDLL\\('${LIB_REALPATH}'\\)@CDLL\\('${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/libadjoint.so'\\)@" ${CMAKE_BINARY_DIR}/clibadjoint.py
                                    COMMENT "Adaptations for clibadjoint.py..."
                                    )

                 # 3.) Generate python/libadjoint/clibadjoint_constants.py via
                 #     a Python script.
                 add_custom_command(TARGET adjoint
                                    POST_BUILD
                                    COMMAND ${PYTHON_EXECUTABLE} ${libadjoint_SOURCE_DIR}/tools/create_python_constants.py ${libadjoint_SOURCE_DIR} ${CMAKE_BINARY_DIR}
                                    COMMENT "Generating clibadjoint_constants.py..."
                                    )

# Make sure to install the generated files
install(FILES ${CMAKE_BINARY_DIR}/clibadjoint.py
              ${CMAKE_BINARY_DIR}/clibadjoint_constants.py
              DESTINATION lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/libadjoint/
        )
# -----------------------------------------------------------------------------
