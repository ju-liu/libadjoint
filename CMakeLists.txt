cmake_minimum_required(VERSION 2.8)

project(libadjoint C CXX Fortran)

# Set the version. (Used for so-naming of the library.)
set(libadjoint_MAJOR_VERSION 0)
set(libadjoint_MINOR_VERSION 9)
set(libadjoint_PATCH_VERSION 2)
set(libadjoint_VERSION
    ${libadjoint_MAJOR_VERSION}.${libadjoint_MINOR_VERSION}.${libadjoint_PATCH_VERSION})

# Build shared library.
SET(BUILD_SHARED_LIBS ON)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR "lib" CACHE PATH
    "Installation directory for libraries")
set(INSTALL_BIN_DIR "bin" CACHE PATH
    "Installation directory for executables")
set(INSTALL_INCLUDE_DIR "include/libadjoint" CACHE PATH
    "Installation directory for header files")

add_subdirectory(src)
add_subdirectory(python)

# -----------------------------------------------------------------------------
# Now some hardcore hacks for the Python bindings.
find_package(PythonInterp REQUIRED)
# 1.) Generate pylibadjoint.h and mess around with it.
execute_process(COMMAND sh -c "${CMAKE_C_COMPILER} -E -U__BLOCKS__ -DPYTHON_BINDINGS ${libadjoint_SOURCE_DIR}/include/libadjoint/libadjoint.h > ${libadjoint_SOURCE_DIR}/include/libadjoint/pylibadjoint.h")
execute_process(COMMAND sh -c "sed -i s/__builtin___stpncpy_chk/__builtin___strncpy_chk/ ${libadjoint_SOURCE_DIR}/include/libadjoint/pylibadjoint.h")
execute_process(COMMAND sh -c "${libadjoint_SOURCE_DIR}/python/ctypeslib/scripts/h2xml.py -q -I. ${libadjoint_SOURCE_DIR}/include/libadjoint/pylibadjoint.h -o ${libadjoint_SOURCE_DIR}/python/libadjoint/libadjoint.xml")
execute_process(COMMAND sh -c "rm -f ${libadjoint_SOURCE_DIR}/include/libadjoint/pylibadjoint.h")
get_filename_component(LIB_REALPATH ${libadjoint_SOURCE_DIR}/lib/libadjoint.so
                       REALPATH
                       )
execute_process(COMMAND sh -c "${libadjoint_SOURCE_DIR}/python/ctypeslib/scripts/xml2py.py -r '^adj.*' -l ${LIB_REALPATH} ${libadjoint_SOURCE_DIR}/python/libadjoint/libadjoint.xml -o ${libadjoint_SOURCE_DIR}/python/libadjoint/clibadjoint.py")
# Some extra adaptations for OS X.
IF(APPLE)
  execute_process(COMMAND sh -c "sed -i \"s/._pack_ = 4/._pack_ = 8/\" ${libadjoint_SOURCE_DIR}/python/libadjoint/clibadjoint.py")
ENDIF()
execute_process(COMMAND sh -c "rm -f ${libadjoint_SOURCE_DIR}/python/libadjoint/libadjoint.xml")
execute_process(COMMAND sh -c "chmod a-x ${libadjoint_SOURCE_DIR}/python/libadjoint/clibadjoint.py")

# 2.) Make sure that the installed clibadjoint.py has the location of the
#     shared library hardcoded in it.
execute_process(COMMAND sh -c "sed -i -e \"s@CDLL('${libadjoint_SOURCE_DIR}/lib/libadjoint.so')@CDLL('${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/libadjoint.so')@\" -e \"s@CDLL('${LIB_REALPATH}')@CDLL('${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/libadjoint.so')@\" ${libadjoint_SOURCE_DIR}/python/libadjoint/clibadjoint.py")

# 3.) Generate python/libadjoint/clibadjoint_constants.py via a Python script.
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${libadjoint_SOURCE_DIR}/tools/create_python_constants.py ${libadjoint_SOURCE_DIR})
# -----------------------------------------------------------------------------
