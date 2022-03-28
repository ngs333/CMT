cmake_minimum_required(VERSION 3.10.2)	
project(CMT)	

set(CMAKE_CXX_STANDARD 17)	
set(CMAKE_CXX_STANDARD_REQUIRED YES)	
#set_property(TARGET tgt PROPERTY CXX_STANDARD 17)	
#[[	
message(STATUS "Setting MSVC flags")	
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc /std:c++latest")	
]]	



set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")	
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")	
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")	

#find_package(Boost REQUIRED)	

add_subdirectory(common)	
add_subdirectory(metric)	
add_subdirectory(search)	
add_subdirectory(apps)	
add_subdirectory(vendor/edlib EXCLUDE_FROM_ALL)	

add_dependencies(libMetric libCommon)	
add_dependencies(libSearch libCommon)	
add_dependencies(EuclidSearchApp  libCommon libMetric libSearch)	
add_dependencies(EditDistSearchApp  libCommon libMetric libSearch)	
add_dependencies(ExtractFastaApp  libCommon libMetric libSearch)	
add_dependencies(HMSearchApp      libCommon libMetric libSearch)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT tapps)
