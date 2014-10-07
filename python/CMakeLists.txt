# Python bindings.
# <http://bloerg.net/2012/11/10/cmake-and-distutils.html>
find_package(PythonInterp 2.7)
if(PYTHONINTERP_FOUND)
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(CLIBADJOINT_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/libadjoint/clibadjoint.py.in")
    set(CLIBADJOINT_PY "${CMAKE_CURRENT_SOURCE_DIR}/libadjoint/clibadjoint.py")
    set(DEPS   "${CMAKE_CURRENT_SOURCE_DIR}/libadjoint/__init__.py")
    set(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build")

    # Generate setup.py
    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    # Generate clibadjoint.py
    configure_file(${CLIBADJOINT_PY_IN} ${CLIBADJOINT_PY})

    # Run `python setup.py -q build`
    add_custom_command(
      OUTPUT ${OUTPUT}/timestamp
      COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY}
      ARGS  build
      COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
      DEPENDS ${DEPS}
    )

    add_custom_target(target ALL DEPENDS ${OUTPUT}/timestamp)

    # Run `python setup.py install`
    install(CODE
            "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX})"
           )
endif()
