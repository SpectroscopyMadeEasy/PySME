cmake_minimum_required(VERSION 3.18)
project(pysme_astro LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

set(PLATFORM "${CMAKE_SYSTEM_NAME}")
string(TOLOWER "${PLATFORM}" PLATFORM)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/platform.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/platform.h
    @ONLY
)

# Prefer submodule checkout in this repository, fall back to sibling SMElib tree.
set(SMELIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/smelib")
if(NOT EXISTS "${SMELIB_DIR}/src/sme/sme_synth_faster.cpp")
    set(SMELIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../SMElib")
endif()
if(NOT EXISTS "${SMELIB_DIR}/src/sme/sme_synth_faster.cpp")
    message(FATAL_ERROR "Could not find SMElib sources. Run with submodules or provide ../SMElib.")
endif()

set(SMELIB_SOURCES
    ${SMELIB_DIR}/src/sme/sme_synth_faster.cpp
    ${SMELIB_DIR}/src/eos/eos.f
    ${SMELIB_DIR}/src/eos/eos_eqns.f
    ${SMELIB_DIR}/src/eos/eos_math_special.f
    ${SMELIB_DIR}/src/eos/kp_q_spln.f
    ${SMELIB_DIR}/src/sme/hlinop.f
    ${SMELIB_DIR}/src/sme/hlinprof.f
)

add_library(sme SHARED ${SMELIB_SOURCES})
target_include_directories(sme PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    ${SMELIB_DIR}/src/sme
)
target_compile_options(sme PRIVATE -O3)

# Keep existing runtime naming convention expected by pysme.smelib.libtools.
if(WIN32)
    set_target_properties(sme PROPERTIES
        PREFIX ""
        OUTPUT_NAME "libsme-5"
    )
elseif(APPLE)
    set_target_properties(sme PROPERTIES
        OUTPUT_NAME "sme"
        INSTALL_RPATH "@loader_path"
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_NAME_DIR "@rpath"
        MACOSX_RPATH TRUE
    )
else()
    set(SME_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
    if(SME_ARCH STREQUAL "amd64")
        set(SME_ARCH "x86_64")
    endif()
    set_target_properties(sme PROPERTIES
        PREFIX ""
        SUFFIX ""
        OUTPUT_NAME "sme_synth.so.linux.${SME_ARCH}.64"
        INSTALL_RPATH "$ORIGIN"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
endif()

Python_add_library(_smelib MODULE src/pysme/smelib/_smelib.cpp)
target_include_directories(_smelib PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    ${SMELIB_DIR}/src/sme
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pysme/smelib
    ${Python_NumPy_INCLUDE_DIRS}
)
target_link_libraries(_smelib PRIVATE sme)

if(APPLE)
    set_target_properties(_smelib PROPERTIES
        INSTALL_RPATH "@loader_path/../lib"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
elseif(UNIX)
    set_target_properties(_smelib PROPERTIES
        INSTALL_RPATH "$ORIGIN/../lib"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
endif()

install(TARGETS sme
    LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/pysme/lib
    RUNTIME DESTINATION ${SKBUILD_PLATLIB_DIR}/pysme/lib
)

install(TARGETS _smelib
    LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/pysme/smelib
)

install(DIRECTORY ${SMELIB_DIR}/src/data/
    DESTINATION ${SKBUILD_PLATLIB_DIR}/pysme/share/libsme
)

install(FILES
    src/pysme/config_default.json
    src/pysme/datafiles_atmospheres.json
    src/pysme/datafiles_nlte.json
    src/pysme/lineprof.dat.gz
    DESTINATION ${SKBUILD_PLATLIB_DIR}/pysme
)
