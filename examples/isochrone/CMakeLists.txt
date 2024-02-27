example_exe(isochrone)

add_custom_command(TARGET isochrone POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/start.dat
    ${CMAKE_CURRENT_BINARY_DIR}/start.dat)
