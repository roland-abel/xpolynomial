# Creates example executables

add_executable(quadratic_roots quadratic_roots.cpp)
target_include_directories(quadratic_roots PRIVATE ../include)

add_executable(cubic_roots cubic_roots.cpp)
target_include_directories(cubic_roots PRIVATE ../include)

add_executable(evaluate evaluate.cpp)
target_include_directories(evaluate PRIVATE ../include)

add_executable(sturm_sequence sturm_sequence.cpp)
target_include_directories(sturm_sequence PRIVATE ../include)

add_executable(root_finder root_finder.cpp)
target_include_directories(root_finder PRIVATE ../include)

add_executable(chebyshev_1st chebyshev_1st.cpp)
target_include_directories(chebyshev_1st PRIVATE ../include)