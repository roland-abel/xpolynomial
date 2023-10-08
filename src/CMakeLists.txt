
set(SOURCE
    main.cpp
    main.h
)

# Creates a static library
add_library(${LIBRARY_NAME} STATIC ${SOURCE})

target_include_directories(${LIBRARY_NAME} PUBLIC ../include)

