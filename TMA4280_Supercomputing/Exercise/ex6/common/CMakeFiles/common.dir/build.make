# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common

# Include any dependencies generated for this target.
include CMakeFiles/common.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/common.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/common.dir/flags.make

CMakeFiles/common.dir/common.c.o: CMakeFiles/common.dir/flags.make
CMakeFiles/common.dir/common.c.o: common.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/common.dir/common.c.o"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/common.dir/common.c.o   -c /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common/common.c

CMakeFiles/common.dir/common.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/common.dir/common.c.i"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common/common.c > CMakeFiles/common.dir/common.c.i

CMakeFiles/common.dir/common.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/common.dir/common.c.s"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common/common.c -o CMakeFiles/common.dir/common.c.s

CMakeFiles/common.dir/common.c.o.requires:
.PHONY : CMakeFiles/common.dir/common.c.o.requires

CMakeFiles/common.dir/common.c.o.provides: CMakeFiles/common.dir/common.c.o.requires
	$(MAKE) -f CMakeFiles/common.dir/build.make CMakeFiles/common.dir/common.c.o.provides.build
.PHONY : CMakeFiles/common.dir/common.c.o.provides

CMakeFiles/common.dir/common.c.o.provides.build: CMakeFiles/common.dir/common.c.o

# Object files for target common
common_OBJECTS = \
"CMakeFiles/common.dir/common.c.o"

# External object files for target common
common_EXTERNAL_OBJECTS =

libcommon.a: CMakeFiles/common.dir/common.c.o
libcommon.a: CMakeFiles/common.dir/build.make
libcommon.a: CMakeFiles/common.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C static library libcommon.a"
	$(CMAKE_COMMAND) -P CMakeFiles/common.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/common.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/common.dir/build: libcommon.a
.PHONY : CMakeFiles/common.dir/build

CMakeFiles/common.dir/requires: CMakeFiles/common.dir/common.c.o.requires
.PHONY : CMakeFiles/common.dir/requires

CMakeFiles/common.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/common.dir/cmake_clean.cmake
.PHONY : CMakeFiles/common.dir/clean

CMakeFiles/common.dir/depend:
	cd /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common /home/ahoermer/NTNU2014/TMA4280_Supercomputing/Exercise/ex4/common/CMakeFiles/common.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/common.dir/depend

