# Creates example executables

set(EXAMPLE_1 exsample1.exe)
add_executable(${EXAMPLE_1} example1.cpp)
target_include_directories(${EXAMPLE_1} PRIVATE ../include)

set(EXAMPLE_2 exsample2.exe)
add_executable(${EXAMPLE_2} example2.cpp)
target_include_directories(${EXAMPLE_2} PRIVATE ../include)
