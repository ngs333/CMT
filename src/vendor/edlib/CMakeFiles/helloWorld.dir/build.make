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
CMAKE_BINARY_DIR = /home/wg25r/new/CMT/src

# Include any dependencies generated for this target.
include vendor/edlib/CMakeFiles/helloWorld.dir/depend.make

# Include the progress variables for this target.
include vendor/edlib/CMakeFiles/helloWorld.dir/progress.make

# Include the compile flags for this target's objects.
include vendor/edlib/CMakeFiles/helloWorld.dir/flags.make

vendor/edlib/CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.o: vendor/edlib/CMakeFiles/helloWorld.dir/flags.make
vendor/edlib/CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.o: vendor/edlib/apps/hello-world/helloWorld.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wg25r/new/CMT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object vendor/edlib/CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.o"
	cd /home/wg25r/new/CMT/src/vendor/edlib && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.o -c /home/wg25r/new/CMT/src/vendor/edlib/apps/hello-world/helloWorld.c

vendor/edlib/CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.i"
	cd /home/wg25r/new/CMT/src/vendor/edlib && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/wg25r/new/CMT/src/vendor/edlib/apps/hello-world/helloWorld.c > CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.i

vendor/edlib/CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.s"
	cd /home/wg25r/new/CMT/src/vendor/edlib && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/wg25r/new/CMT/src/vendor/edlib/apps/hello-world/helloWorld.c -o CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.s

# Object files for target helloWorld
helloWorld_OBJECTS = \
"CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.o"

# External object files for target helloWorld
helloWorld_EXTERNAL_OBJECTS =

bin/helloWorld: vendor/edlib/CMakeFiles/helloWorld.dir/apps/hello-world/helloWorld.c.o
bin/helloWorld: vendor/edlib/CMakeFiles/helloWorld.dir/build.make
bin/helloWorld: lib/libedlib.a
bin/helloWorld: vendor/edlib/CMakeFiles/helloWorld.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/wg25r/new/CMT/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/helloWorld"
	cd /home/wg25r/new/CMT/src/vendor/edlib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/helloWorld.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
vendor/edlib/CMakeFiles/helloWorld.dir/build: bin/helloWorld

.PHONY : vendor/edlib/CMakeFiles/helloWorld.dir/build

vendor/edlib/CMakeFiles/helloWorld.dir/clean:
	cd /home/wg25r/new/CMT/src/vendor/edlib && $(CMAKE_COMMAND) -P CMakeFiles/helloWorld.dir/cmake_clean.cmake
.PHONY : vendor/edlib/CMakeFiles/helloWorld.dir/clean

vendor/edlib/CMakeFiles/helloWorld.dir/depend:
	cd /home/wg25r/new/CMT/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wg25r/new/CMT/src /home/wg25r/new/CMT/src/vendor/edlib /home/wg25r/new/CMT/src /home/wg25r/new/CMT/src/vendor/edlib /home/wg25r/new/CMT/src/vendor/edlib/CMakeFiles/helloWorld.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : vendor/edlib/CMakeFiles/helloWorld.dir/depend

