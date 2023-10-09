# Creates example executables

add_executable(quadratic_roots quadratic_roots.cpp)
target_include_directories(quadratic_roots PRIVATE ../include)

add_executable(cubic_roots cubic_roots.cpp)
target_include_directories(cubic_roots PRIVATE ../include)

add_executable(evaluate evaluate.cpp)
target_include_directories(evaluate PRIVATE ../include)

