cmake_minimum_required(VERSION 3.10)
project(samvt C)

set(CMAKE_C_STANDARD 99)

add_library(libsamvt STATIC coverage.c bam.c fa.c)
set_target_properties(libsamvt PROPERTIES OUTPUT_NAME samvt)
target_link_libraries(libsamvt pthread htsm bigWig z curl mt)

add_executable(samvt main.c common.c samvt_coverage.c samvt_mutation.c)
target_link_libraries(samvt libsamvt)
install(TARGETS samvt RUNTIME DESTINATION bin)

add_subdirectory(mt)

include_directories("./htslib-1.15_modified/")
add_custom_target(make_htsm ALL
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib-1.15_modified/
        COMMENT "generate modified modified version of htslib.")
add_library(htsm SHARED IMPORTED)
set_target_properties(htsm PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/htslib-1.15_modified/libhtsm.so)
add_dependencies(htsm make_htsm)
add_dependencies(libsamvt htsm)

include_directories("./libBigWig-0.4.4/")
add_custom_target(make_bigWig ALL
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/libBigWig-0.4.4/
        COMMENT "generate modified libBigWig.")
add_library(bigWig STATIC IMPORTED)
set_target_properties(bigWig PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/libBigWig-0.4.4/libBigWig.a)
add_dependencies(bigWig make_bigWig)
add_dependencies(libsamvt bigWig)
