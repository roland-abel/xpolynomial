# Creates example executables

file(GLOB EXAMPLES "**.cpp")

foreach (EXAMPLE ${EXAMPLES})
    get_filename_component(EXAMPLE_NAME ${EXAMPLE} NAME_WE)

    add_executable(${EXAMPLE_NAME} ${EXAMPLE})
    target_link_libraries(${EXAMPLE_NAME} PRIVATE xpolynomial::xpolynomial)
endforeach()
