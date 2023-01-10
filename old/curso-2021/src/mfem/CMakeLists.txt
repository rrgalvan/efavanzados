# Fichero descargado de https://github.com/mfem/mfem/issues/1163
cmake_minimum_required(VERSION 2.8.12)
project(MFEM_first_test)

set(CMAKE_CXX_STANDARD 11)

# # Prohibit in-source builds
# if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
#   message(FATAL_ERROR "In-source builds are prohibited.")
# endif ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

# Import MFEM. The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#  * mfem_DIR - absolute path to where MFEMConfig.cmake is.
message(STATUS "Looking for mfem...")
set(MFEM_DIR "" CACHE PATH "Path to the MFEM build or install prefix.")
if (MFEM_DIR)
   find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
                "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)
else()
   find_package(mfem REQUIRED NAMES MFEM)
endif()
message(STATUS "Found mfem config in: ${mfem_DIR} (version ${MFEM_VERSION})")
# Use the same C++ compiler as MFEM. This is needed when MFEM was built using
# an MPI wrapper and we do not have explicitly the MPI compile and link flags.
message(STATUS "Mfem compiler is: ${MFEM_CXX_COMPILER}")
if (NOT CMAKE_CXX_COMPILER AND MFEM_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER "${MFEM_CXX_COMPILER}")
endif()

enable_language(C)
enable_language(CXX)

include_directories(${MFEM_INCLUDE_DIRS})

add_executable(MFEM-first-test MFEM-first-test.cpp)
target_link_libraries(MFEM-first-test PRIVATE ${MFEM_LIBRARIES})
add_custom_command(TARGET MFEM-first-test POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy
                    ${CMAKE_SOURCE_DIR}/star.mesh $<TARGET_FILE_DIR:MFEM-first-test>)
