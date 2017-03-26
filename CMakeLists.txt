# vim: expandtab:ts=2:sw=2
cmake_minimum_required(VERSION 3.2.0)
project(mcf)

if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif ()
if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
endif()
if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
endif()
if (NOT PYTHON_OUTPUT_DIRECTORY)
  set(PYTHON_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python_lib")
endif()

set(MCF_EXTERNAL_INSTALL_DIR
  "${CMAKE_BINARY_DIR}/external-install"
  CACHE PATH "External depedency installation directory"
  )

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${MCF_EXTERNAL_INSTALL_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(cmake/DefineCompilerFlags.cmake)


find_package(Boost QUIET)
if (Boost_FOUND)
  option(MCF_USE_BOOST "Use Boost library" ON)
else ()
  option(MCF_USE_BOOST "Use Boost library" OFF)
endif ()

find_package(LEMON QUIET)
if (LEMON_FOUND)
  option(MCF_USE_Lemon "Build COIN-OR Lemon solver" ON)
else ()
  option(MCF_USE_Lemon "Build COIN-OR Lemon solver" OFF)
endif ()

find_package(Clp QUIET)
if (Clp_FOUND)
  option(MCF_USE_Clp "Build COIN-OR CLP solver" ON)
else ()
  option(MCF_USE_Clp "Build COIN-OR CLP solver" OFF)
endif ()

find_package(pybind11 QUIET)
if (pybind11_FOUND)
  option(MCF_BUILD_PYTHON "Build python module" ON)
else ()
  option(MCF_BUILD_PYTHON "Build python module" OFF)
endif ()

option(MCF_BUILD_STATIC "Build static library" ON)
option(MCF_BUILD_SHARED "Build shared library" OFF)
option(MCF_BUILD_EXAMPLES "Build examples" ON)

set(MCF_LIB_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcf/graph.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/graph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcf/internal/util.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcf/impl/util_inl.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/util.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcf/internal/k_shortest_path.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/k_shortest_path.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcf/logging.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcf/impl/logging_inl.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcf/solver.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcf/clp_solver.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/clp_solver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcf/k_shortest_path_solver.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/k_shortest_path_solver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcf/lemon_solver.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/lemon_solver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/mcf/batch_processing.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/batch_processing.cpp
  )

set(MCF_DEFINITIONS "")
set(MCF_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(MCF_LIBRARIES "")

if (MCF_USE_BOOST)
  find_package(Boost REQUIRED)
  set(MCF_DEFINITIONS ${MCF_DEFINITIONS} -DMCF_USE_Boost)
  set(MCF_INCLUDE_DIRS ${MCF_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
endif ()
if (MCF_USE_Lemon)
  find_package(LEMON REQUIRED)
  set(MCF_DEFINITIONS ${MCF_DEFINITIONS} -DMCF_USE_Lemon)
  set(MCF_INCLUDE_DIRS ${MCF_INCLUDE_DIRS} ${LEMON_INCLUDE_DIRS})
  set(MCF_LIBRARIES ${MCF_LIBRARIES} ${LEMON_LIBRARIES})
endif ()
if (MCF_USE_Clp)
  find_package(Clp REQUIRED)
  find_package(LAPACK REQUIRED)
  find_package(ZLIB REQUIRED)
  find_package(BZip2 REQUIRED)
  set(MCF_DEFINITIONS ${MCF_DEFINITIONS} -DMCF_USE_Clp)
  set(MCF_INCLUDE_DIRS
    ${MCF_INCLUDE_DIRS}
    ${Clp_INCLUDE_DIRS}
    ${LAPACK_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    ${BZIP2_INCLUDE_DIRS}
    )
  set(MCF_LIBRARIES
    ${MCF_LIBRARIES}
    ${Clp_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${BZIP2_LIBRARIES}
    )
endif ()

include_directories(${MCF_INCLUDE_DIRS})
add_definitions(${MCF_DEFINITIONS})

find_program(CLANG_FORMAT_EXECUTABLE
  NAMES clang-format clang-format-3.5 clang-format-3.6 clang-format-3.7
        clang-format-3.8 clang-format-3.9
  DOC "Path to clang-format executable"
  )
if (CLANG_FORMAT_EXECUTABLE)
  set(CLANG_FORMAT_ARGS "${CLANG_FORMAT_EXECUTABLE};-i;-style=Google")
  set(FILES_TO_FORMAT
    ${MCF_LIB_SRC}
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/min_cost_flow_example.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/batch_processing_example.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/pymcf.cpp
    )
  add_custom_target(${PROJECT_NAME}_format
    COMMAND ${CLANG_FORMAT_ARGS} ${FILES_TO_FORMAT})
else ()
  set(CLANG_FORMAT_ARGS "${CLANG_FORMAT_EXECUTABLE};-i;-style=Google")
  add_custom_target(${PROJECT_NAME}_format COMMAND "")
endif ()

if (MCF_BUILD_STATIC)
  add_library(${PROJECT_NAME}_s STATIC ${MCF_LIB_SRC})
  add_dependencies(${PROJECT_NAME}_s ${PROJECT_NAME}_format)
endif ()

if (MCF_BUILD_SHARED)
  add_library(${PROJECT_NAME} SHARED ${MCF_LIB_SRC})
  target_link_libraries(${PROJECT_NAME} ${MCF_LIBRARIES})
  add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}_format)
endif ()

if (MCF_BUILD_EXAMPLES)
  add_executable(min_cost_flow_example
    ${MCF_LIB_SRC}
    examples/min_cost_flow_example.cpp
  )
  target_link_libraries(min_cost_flow_example ${MCF_LIBRARIES})
  add_dependencies(min_cost_flow_example ${PROJECT_NAME}_format)

  add_executable(batch_processing_example
    ${MCF_LIB_SRC}
    examples/batch_processing_example.cpp
  )
  target_link_libraries(batch_processing_example ${MCF_LIBRARIES})
  add_dependencies(batch_processing_example ${PROJECT_NAME}_format)
endif ()

if (MCF_BUILD_PYTHON)
  find_package(pybind11 REQUIRED)

  pybind11_add_module(pymcf ${MCF_LIB_SRC} python/pymcf.cpp)
  target_link_libraries(pymcf ${MCF_LIBRARIES})
  set_target_properties(pymcf PROPERTIES OUTPUT_NAME mcf)
  set_target_properties(pymcf PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PYTHON_OUTPUT_DIRECTORY}
    )
endif ()
