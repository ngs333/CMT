add_library(libCommon  INTERFACE)

include_directories(${CMAKE_SOURCE_DIR}/../include)

target_include_directories(libCommon INTERFACE ${CMAKE_SOURCE_DIR}/common)
