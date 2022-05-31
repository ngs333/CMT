#set(staticlib "-static-libgcc -static-libstdc++")

include_directories(${CMAKE_SOURCE_DIR}/../include)
#parassail lib include location: 
include_directories(/home/mzuniga/usr/local/include)
#include_directories(${CMAKE_SOURCE_DIR}/common)
#include_directories(${CMAKE_SOURCE_DIR}/metric)
#include_directories(${CMAKE_SOURCE_DIR}/search)

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#parsail lib lib location:
link_directories(/home/mzuniga/usr/local/lib)

add_executable(EuclidSearchApp EuclidSearchMain.cpp)
add_executable(EditDistSearchApp EditDistSearchMain.cpp)
add_executable(ExtractFastaApp ExtractFasta.cpp)
add_executable(InDimCalculatorApp InDimCalculator.cpp)
add_executable(HMSearchApp HMSearchMain.cpp)
add_executable(GEDSearchApp GEDSearchMain.cpp)

target_link_libraries(InDimCalculatorApp libCommon)
target_link_libraries(InDimCalculatorApp libMetric)
target_link_libraries(InDimCalculatorApp libSearch)
target_link_libraries(InDimCalculatorApp edlib)

target_link_libraries(ExtractFastaApp libCommon)
target_link_libraries(ExtractFastaApp libMetric)
target_link_libraries(ExtractFastaApp libSearch)
target_link_libraries(ExtractFastaApp edlib)

target_link_libraries(EuclidSearchApp libCommon)
target_link_libraries(EuclidSearchApp libMetric)
target_link_libraries(EuclidSearchApp libSearch)
target_link_libraries(EuclidSearchApp edlib)


target_link_libraries(EditDistSearchApp  libCommon)
target_link_libraries(EditDistSearchApp libMetric)
target_link_libraries(EditDistSearchApp libSearch)
target_link_libraries(EditDistSearchApp edlib)

target_link_libraries(HMSearchApp libCommon)
target_link_libraries(HMSearchApp libMetric)
target_link_libraries(HMSearchApp libSearch)
target_link_libraries(HMSearchApp edlib)

target_link_libraries(GEDSearchApp libCommon)
target_link_libraries(GEDSearchApp libMetric)
target_link_libraries(GEDSearchApp libSearch)
target_link_libraries(GEDSearchApp edlib)
