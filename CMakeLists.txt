cmake_minimum_required(VERSION 2.8)

project(libadjoint C CXX)

# Set the version. (Used for so-naming of the library.)
set(libadjoint_MAJOR_VERSION 1)
set(libadjoint_MINOR_VERSION 6)
set(libadjoint_PATCH_VERSION 0)
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

# Fix random bug in CMake (fix from Jack Hale)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_subdirectory(src)
add_subdirectory(python)
