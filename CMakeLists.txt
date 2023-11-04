cmake_minimum_required(VERSION 3.6)
project(openmp_test_clion)
set(CMAKE_CXX_STANDARD 11)

# Set the default number of threads
set(DEFAULT_NUM_THREADS 12)

# Create a CMake variable to configure the number of threads
set(NUM_THREADS ${DEFAULT_NUM_THREADS} CACHE STRING "Number of threads to use")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Pass the number of threads as a compile definition to the source code
add_executable(openmp_test_clion ompcInsertion.c)
target_compile_definitions(openmp_test_clion PRIVATE NUM_THREADS=${NUM_THREADS})

# Add the "coords" files to the source files
# Copy the .coords files to the build directory
set(COORDS_FILES
        "9_coords.coord"
        "16_coords.coord"
        "4096_coords.coord"
)

foreach(COORD_FILE ${COORDS_FILES})
    add_custom_command(TARGET openmp_test_clion PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_SOURCE_DIR}/${COORD_FILE}"
            "${CMAKE_BINARY_DIR}/${COORD_FILE}"
    )
endforeach()

find_package(OpenMP)
IF(OPENMP_FOUND)
    message("Found OpenMP!")
    # add flags for OpenMP
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ELSE()
    message("Missed OpenMP!")
ENDIF()
