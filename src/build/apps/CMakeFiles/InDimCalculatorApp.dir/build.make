# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/wg25r/new/CMT/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/wg25r/new/CMT/src/build

# Include any dependencies generated for this target.
include apps/CMakeFiles/InDimCalculatorApp.dir/depend.make

# Include the progress variables for this target.
include apps/CMakeFiles/InDimCalculatorApp.dir/progress.make

# Include the compile flags for this target's objects.
include apps/CMakeFiles/InDimCalculatorApp.dir/flags.make

apps/CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.o: apps/CMakeFiles/InDimCalculatorApp.dir/flags.make
apps/CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.o: ../apps/InDimCalculator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wg25r/new/CMT/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object apps/CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.o"
	cd /home/wg25r/new/CMT/src/build/apps && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.o -c /home/wg25r/new/CMT/src/apps/InDimCalculator.cpp

apps/CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.i"
	cd /home/wg25r/new/CMT/src/build/apps && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wg25r/new/CMT/src/apps/InDimCalculator.cpp > CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.i

apps/CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.s"
	cd /home/wg25r/new/CMT/src/build/apps && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wg25r/new/CMT/src/apps/InDimCalculator.cpp -o CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.s

# Object files for target InDimCalculatorApp
InDimCalculatorApp_OBJECTS = \
"CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.o"

# External object files for target InDimCalculatorApp
InDimCalculatorApp_EXTERNAL_OBJECTS =

/home/wg25r/new/CMT/bin/InDimCalculatorApp: apps/CMakeFiles/InDimCalculatorApp.dir/InDimCalculator.cpp.o
/home/wg25r/new/CMT/bin/InDimCalculatorApp: apps/CMakeFiles/InDimCalculatorApp.dir/build.make
/home/wg25r/new/CMT/bin/InDimCalculatorApp: lib/libedlib.a
/home/wg25r/new/CMT/bin/InDimCalculatorApp: apps/CMakeFiles/InDimCalculatorApp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/wg25r/new/CMT/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable /home/wg25r/new/CMT/bin/InDimCalculatorApp"
	cd /home/wg25r/new/CMT/src/build/apps && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/InDimCalculatorApp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/CMakeFiles/InDimCalculatorApp.dir/build: /home/wg25r/new/CMT/bin/InDimCalculatorApp

.PHONY : apps/CMakeFiles/InDimCalculatorApp.dir/build

apps/CMakeFiles/InDimCalculatorApp.dir/clean:
	cd /home/wg25r/new/CMT/src/build/apps && $(CMAKE_COMMAND) -P CMakeFiles/InDimCalculatorApp.dir/cmake_clean.cmake
.PHONY : apps/CMakeFiles/InDimCalculatorApp.dir/clean

apps/CMakeFiles/InDimCalculatorApp.dir/depend:
	cd /home/wg25r/new/CMT/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wg25r/new/CMT/src /home/wg25r/new/CMT/src/apps /home/wg25r/new/CMT/src/build /home/wg25r/new/CMT/src/build/apps /home/wg25r/new/CMT/src/build/apps/CMakeFiles/InDimCalculatorApp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/CMakeFiles/InDimCalculatorApp.dir/depend

