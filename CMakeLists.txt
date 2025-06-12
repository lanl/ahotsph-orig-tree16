cmake_minimum_required(VERSION 3.10)
project(ptw C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Add source files
add_executable(ptw 
    main.c
    ptw.c
    test_ptw.c
)

# Add include directories
target_include_directories(ptw PRIVATE 
    ${CMAKE_SOURCE_DIR}
)

# Link math library since we're using mathematical functions
target_link_libraries(ptw PRIVATE m)