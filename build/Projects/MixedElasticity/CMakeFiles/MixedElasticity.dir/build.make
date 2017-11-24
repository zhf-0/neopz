# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/superdan/Downloads/neopz

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/superdan/Downloads/neopz/build

# Include any dependencies generated for this target.
include Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/depend.make

# Include the progress variables for this target.
include Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/progress.make

# Include the compile flags for this target's objects.
include Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/flags.make

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/flags.make
Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o: ../Projects/MixedElasticity/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/superdan/Downloads/neopz/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o"
	cd /home/superdan/Downloads/neopz/build/Projects/MixedElasticity && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MixedElasticity.dir/main.cpp.o -c /home/superdan/Downloads/neopz/Projects/MixedElasticity/main.cpp

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MixedElasticity.dir/main.cpp.i"
	cd /home/superdan/Downloads/neopz/build/Projects/MixedElasticity && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/superdan/Downloads/neopz/Projects/MixedElasticity/main.cpp > CMakeFiles/MixedElasticity.dir/main.cpp.i

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MixedElasticity.dir/main.cpp.s"
	cd /home/superdan/Downloads/neopz/build/Projects/MixedElasticity && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/superdan/Downloads/neopz/Projects/MixedElasticity/main.cpp -o CMakeFiles/MixedElasticity.dir/main.cpp.s

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o.requires:

.PHONY : Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o.requires

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o.provides: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o.requires
	$(MAKE) -f Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/build.make Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o.provides.build
.PHONY : Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o.provides

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o.provides.build: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o


Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/flags.make
Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o: ../Projects/MixedElasticity/meshgen.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/superdan/Downloads/neopz/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o"
	cd /home/superdan/Downloads/neopz/build/Projects/MixedElasticity && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MixedElasticity.dir/meshgen.cpp.o -c /home/superdan/Downloads/neopz/Projects/MixedElasticity/meshgen.cpp

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MixedElasticity.dir/meshgen.cpp.i"
	cd /home/superdan/Downloads/neopz/build/Projects/MixedElasticity && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/superdan/Downloads/neopz/Projects/MixedElasticity/meshgen.cpp > CMakeFiles/MixedElasticity.dir/meshgen.cpp.i

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MixedElasticity.dir/meshgen.cpp.s"
	cd /home/superdan/Downloads/neopz/build/Projects/MixedElasticity && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/superdan/Downloads/neopz/Projects/MixedElasticity/meshgen.cpp -o CMakeFiles/MixedElasticity.dir/meshgen.cpp.s

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o.requires:

.PHONY : Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o.requires

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o.provides: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o.requires
	$(MAKE) -f Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/build.make Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o.provides.build
.PHONY : Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o.provides

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o.provides.build: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o


Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/flags.make
Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o: ../Projects/MixedElasticity/pzmixedelasmat.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/superdan/Downloads/neopz/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o"
	cd /home/superdan/Downloads/neopz/build/Projects/MixedElasticity && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o -c /home/superdan/Downloads/neopz/Projects/MixedElasticity/pzmixedelasmat.cpp

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.i"
	cd /home/superdan/Downloads/neopz/build/Projects/MixedElasticity && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/superdan/Downloads/neopz/Projects/MixedElasticity/pzmixedelasmat.cpp > CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.i

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.s"
	cd /home/superdan/Downloads/neopz/build/Projects/MixedElasticity && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/superdan/Downloads/neopz/Projects/MixedElasticity/pzmixedelasmat.cpp -o CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.s

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o.requires:

.PHONY : Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o.requires

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o.provides: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o.requires
	$(MAKE) -f Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/build.make Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o.provides.build
.PHONY : Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o.provides

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o.provides.build: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o


# Object files for target MixedElasticity
MixedElasticity_OBJECTS = \
"CMakeFiles/MixedElasticity.dir/main.cpp.o" \
"CMakeFiles/MixedElasticity.dir/meshgen.cpp.o" \
"CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o"

# External object files for target MixedElasticity
MixedElasticity_EXTERNAL_OBJECTS =

Projects/MixedElasticity/MixedElasticity: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o
Projects/MixedElasticity/MixedElasticity: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o
Projects/MixedElasticity/MixedElasticity: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o
Projects/MixedElasticity/MixedElasticity: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/build.make
Projects/MixedElasticity/MixedElasticity: lib/libpz.a
Projects/MixedElasticity/MixedElasticity: /usr/lib/x86_64-linux-gnu/libpthread.so
Projects/MixedElasticity/MixedElasticity: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/superdan/Downloads/neopz/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable MixedElasticity"
	cd /home/superdan/Downloads/neopz/build/Projects/MixedElasticity && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MixedElasticity.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/build: Projects/MixedElasticity/MixedElasticity

.PHONY : Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/build

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/requires: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/main.cpp.o.requires
Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/requires: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/meshgen.cpp.o.requires
Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/requires: Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/pzmixedelasmat.cpp.o.requires

.PHONY : Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/requires

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/clean:
	cd /home/superdan/Downloads/neopz/build/Projects/MixedElasticity && $(CMAKE_COMMAND) -P CMakeFiles/MixedElasticity.dir/cmake_clean.cmake
.PHONY : Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/clean

Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/depend:
	cd /home/superdan/Downloads/neopz/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/superdan/Downloads/neopz /home/superdan/Downloads/neopz/Projects/MixedElasticity /home/superdan/Downloads/neopz/build /home/superdan/Downloads/neopz/build/Projects/MixedElasticity /home/superdan/Downloads/neopz/build/Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Projects/MixedElasticity/CMakeFiles/MixedElasticity.dir/depend

