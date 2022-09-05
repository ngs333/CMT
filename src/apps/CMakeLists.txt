#set(staticlib "-static-libgcc -static-libstdc++")

include_directories(${CMAKE_SOURCE_DIR}/../include)

#parassail lib include location: 
include_directories(/home/mzuniga/usr/local/include)

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

#parsail lib lib location:
link_directories(/home/mzuniga/usr/local/lib)

add_executable(EuclidSearchApp EuclidSearchMain.cpp)
add_executable(EditDistSearchApp EditDistSearchMain.cpp)
add_executable(ExtractFastaApp ExtractFasta.cpp)
add_executable(InDimCalculatorApp InDimCalculator.cpp)

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
