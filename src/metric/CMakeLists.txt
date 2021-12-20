add_library(libMetric INTERFACE)

include_directories(${CMAKE_SOURCE_DIR}/../include)

target_include_directories(libMetric INTERFACE ${CMAKE_SOURCE_DIR}/metric)
