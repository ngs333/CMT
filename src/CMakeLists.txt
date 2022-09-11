cmake_minimum_required(VERSION 3.10.2)	
project(CMT)	

set(CMAKE_CXX_STANDARD 17)	
set(CMAKE_CXX_STANDARD_REQUIRED YES)	
#set_property(TARGET tgt PROPERTY CXX_STANDARD 17)	
#[[	
message(STATUS "Setting MSVC flags")	
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc /std:c++latest")	
]]	

option(USE_NEXTAFTER "Use nextafter expanded bounding intervals" ON)
option(USE_HALFINTERVALS "Use half intervals instead of full intervals" OFF)
option(USE_PARASAIL_LIB "Use the parasial library to evaluate bioseqeunce distances" OFF)


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

if (USE_NEXTAFTER)
  add_definitions(-DUSE_NEXTAFTER)
  target_compile_definitions(EditDistSearchApp PUBLIC USE_NEXTAFTER=${USE_NEXTAFTER})
  target_compile_definitions(EuclidSearchApp PUBLIC USE_NEXTAFTER=${USE_NEXTAFTER})
endif()

if (USE_HALF_INTERVALS)
  add_definitions(-DUSE_HALF_INTERVALS)
  target_compile_definitions(EditDistSearchApp PUBLIC USE_HALF_INTERVALS=${USE_HALF_INTERVALS})
  target_compile_definitions(EuclidSearchApp PUBLIC USE_HALF_INTERVALS=${USE_HALF_INTERVALS})
endif()


if (USE_PARASAIL_LIB)
  add_definitions(-DUSE_PARASAIL_LIB)
  target_compile_definitions(EditDistSearchApp PUBLIC USE_PARASAIL_LIB=${USE_PARASAIL_LIB})
  target_compile_definitions(ExtractFastaApp PUBLIC USE_PARASAIL_LIB=${USE_PARASAIL_LIB})
endif()


set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT tapps)
