# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

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

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /export/starexec/sandbox2/solver

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /export/starexec/sandbox2/solver/build

# Include any dependencies generated for this target.
include CMakeFiles/flow_cutter_pace17.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/flow_cutter_pace17.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/flow_cutter_pace17.dir/flags.make

CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o: CMakeFiles/flow_cutter_pace17.dir/flags.make
CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o: ../flowcutter/pace.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /export/starexec/sandbox2/solver/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o"
	/opt/rh/devtoolset-7/root/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o -c /export/starexec/sandbox2/solver/flowcutter/pace.cpp

CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.i"
	/opt/rh/devtoolset-7/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /export/starexec/sandbox2/solver/flowcutter/pace.cpp > CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.i

CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.s"
	/opt/rh/devtoolset-7/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /export/starexec/sandbox2/solver/flowcutter/pace.cpp -o CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.s

CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o.requires:
.PHONY : CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o.requires

CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o.provides: CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o.requires
	$(MAKE) -f CMakeFiles/flow_cutter_pace17.dir/build.make CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o.provides.build
.PHONY : CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o.provides

CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o.provides.build: CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o

# Object files for target flow_cutter_pace17
flow_cutter_pace17_OBJECTS = \
"CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o"

# External object files for target flow_cutter_pace17
flow_cutter_pace17_EXTERNAL_OBJECTS =

flow_cutter_pace17: CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o
flow_cutter_pace17: CMakeFiles/flow_cutter_pace17.dir/build.make
flow_cutter_pace17: flowcutter/liblibflowcutter.a
flow_cutter_pace17: CMakeFiles/flow_cutter_pace17.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable flow_cutter_pace17"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/flow_cutter_pace17.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/flow_cutter_pace17.dir/build: flow_cutter_pace17
.PHONY : CMakeFiles/flow_cutter_pace17.dir/build

CMakeFiles/flow_cutter_pace17.dir/requires: CMakeFiles/flow_cutter_pace17.dir/flowcutter/pace.cpp.o.requires
.PHONY : CMakeFiles/flow_cutter_pace17.dir/requires

CMakeFiles/flow_cutter_pace17.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/flow_cutter_pace17.dir/cmake_clean.cmake
.PHONY : CMakeFiles/flow_cutter_pace17.dir/clean

CMakeFiles/flow_cutter_pace17.dir/depend:
	cd /export/starexec/sandbox2/solver/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /export/starexec/sandbox2/solver /export/starexec/sandbox2/solver /export/starexec/sandbox2/solver/build /export/starexec/sandbox2/solver/build /export/starexec/sandbox2/solver/build/CMakeFiles/flow_cutter_pace17.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/flow_cutter_pace17.dir/depend
