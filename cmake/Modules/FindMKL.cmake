################################################################################
#
# \file      cmake/FindMKL.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the Math Kernel Library from Intel
# \date      Thu 26 Jan 2017 02:05:50 PM MST
#
################################################################################

# Find the Math Kernel Library from Intel
#
#  MKL_FOUND - System has MKL
#  MKL_INCLUDE_DIRS - MKL include files directories
#  MKL_LIBRARIES - The MKL libraries
#  MKL_INTERFACE_LIBRARY - MKL interface library
#  MKL_OMP_LIBRARY - MKL iomp5 library
#  MKL_THREAD_LAYER_LIBRARY - MKL Thread layer library
#  MKL_CORE_LIBRARY - MKL core library
#
#  The environment variables MKLROOT and INTEL are used to find the library.
#  Everything else is ignored.
#
#  Example usage:
#
#  find_package(MKL)
#  if(MKL_FOUND)
#    target_link_libraries(TARGET ${MKL_LIBRARIES})
#  endif()

# If already in cache, be silent
if (MKL_INCLUDE_DIRS AND MKL_LIBRARIES AND MKL_INTERFACE_LIBRARY AND
    MKL_SEQUENTIAL_LAYER_LIBRARY AND MKL_CORE_LIBRARY)
  set (MKL_FIND_QUIETLY TRUE)
endif()

set(INT_LIB "mkl_intel_ilp64")
set(OMP_LIB "iomp5")
set(THR_LIB "mkl_intel_thread")
set(COR_LIB "mkl_core")

find_path(MKL_INCLUDE_DIR NAMES mkl.h HINTS $ENV{MKLROOT}/include)
find_library(MKL_INTERFACE_LIBRARY
             NAMES ${INT_LIB}
             PATHS $ENV{MKLROOT}/lib
                   $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)

find_library(MKL_OMP_LIBRARY
             NAMES ${OMP_LIB}
             PATHS $ENV{MKLROOT}/lib
                   $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)

find_library(MKL_THREAD_LAYER_LIBRARY
             NAMES ${THR_LIB}
             PATHS $ENV{MKLROOT}/lib
                   $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)

find_library(MKL_CORE_LIBRARY
             NAMES ${COR_LIB}
             PATHS $ENV{MKLROOT}/lib
                   $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)

set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
if(MKL_INCLUDE_DIRS)
  include_directories(${MKL_INCLUDE_DIRS})
endif()
#I guess libiomp5 is not quite necessary
if(NOT MKL_OMP_LIBRARY)
  set(MKL_OMP_LIBRARY "")
endif()
set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_OMP_LIBRARY} ${MKL_THREAD_LAYER_LIBRARY} ${MKL_CORE_LIBRARY})


if (NOT MKL_INCLUDE_DIR OR
    NOT MKL_INTERFACE_LIBRARY OR
    NOT MKL_THREAD_LAYER_LIBRARY OR
    NOT MKL_CORE_LIBRARY)
  set(MKL_INCLUDE_DIRS "")
  set(MKL_LIBRARIES "")
  set(MKL_INTERFACE_LIBRARY "")
  set(MKL_OMP_LIBRARY "")
  set(MKL_THREAD_LAYER_LIBRARY "")
  set(MKL_CORE_LIBRARY "")

endif()

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS MKL_INTERFACE_LIBRARY MKL_THREAD_LAYER_LIBRARY MKL_CORE_LIBRARY)

mark_as_advanced(MKL_INCLUDE_DIRS MKL_LIBRARIES MKL_INTERFACE_LIBRARY MKL_SEQUENTIAL_LAYER_LIBRARY MKL_CORE_LIBRARY)