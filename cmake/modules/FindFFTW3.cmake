# - Find fftw3 library
# Find the native FFTW3 includes and library
# This module defines
#  FFTW3_INCLUDE, where to find fftw3.h, etc.
# ---
#  FFTW3_LIBRARY, library to link against to use FFTW3
#  FFTW3_OMP_LIBRARY, library to link against to use FFTW3_omp
#  FFTW3_THREADS_LIBRARY, library to link against to use FFTW3_threads
#  FFTW3_FOUND, if false, do not try to use FFTW3.
#  FFTW3_OMP_FOUND, if false, do not try to use OpenMP FFTW3.
#  FFTW3_THREADS_FOUND, if false, do not try to use threaded FFTW3.
# ---
#  FFTW3L_LIBRARY, library to link against to use FFTW3l
#  FFTW3L_OMP_LIBRARY, library to link against to use FFTW3l_omp
#  FFTW3L_THREADS_LIBRARY, library to link against to use FFTW3l_threads
#  FFTW3L_FOUND, if false, do not try to use FFTW3l.
#  FFTW3L_OMP_FOUND, if false, do not try to use OpenMP FFTW3l.
#  FFTW3L_THREADS_FOUND, if false, do not try to use threaded FFTW3l.
# ---
#  FFTW3F_LIBRARY, library to link against to use FFTW3f
#  FFTW3F_OMP_LIBRARY, library to link against to use FFTW3f_omp
#  FFTW3F_THREADS_LIBRARY, library to link against to use FFTW3f_threads
#  FFTW3F_FOUND, if false, do not try to use FFTW3f.
#  FFTW3F_OMP_FOUND, if false, do not try to use OpenMP FFTW3f.
#  FFTW3F_THREADS_FOUND, if false, do not try to use threaded FFTW3f.
# ---
#  FFTW3Q_LIBRARY, library to link against to use FFTW3q
#  FFTW3Q_OMP_LIBRARY, library to link against to use FFTW3q_omp
#  FFTW3Q_THREADS_LIBRARY, library to link against to use FFTW3q_threads
#  FFTW3Q_FOUND, if false, do not try to use FFTW3q.
#  FFTW3Q_OMP_FOUND, if false, do not try to use OpenMP FFTW3q.
#  FFTW3Q_THREADS_FOUND, if false, do not try to use threaded FFTW3q.

find_package(PkgConfig)

function(add_imported_library lib_name library headers)
  add_library(FFTW3::${lib_name} UNKNOWN IMPORTED)
  set_target_properties(FFTW3::${lib_name} PROPERTIES
    IMPORTED_LOCATION ${library}
    INTERFACE_INCLUDE_DIRECTORIES ${headers}
  )
  set(${lib_name}_FOUND 1 CACHE INTERNAL "FFTW3 ${lib_name} found" FORCE)
  set(${lib_name}_LIBRARY ${library}
      CACHE STRING "Path to FFTW3::${lib_name} library" FORCE)
  set(FFTW3_INCLUDE ${headers}
      CACHE STRING "Path to FFTW3 headers" FORCE)
  mark_as_advanced(FORCE ${lib_name}_LIBRARY)
  mark_as_advanced(FORCE FFTW3_INCLUDE)
endfunction()

#as35 if (FFTW3_LIBRARY AND FFTW3_INCLUDE)
#as35   add_imported_library(${FFTW3_LIBRARY} ${FFTW3_INCLUDE})
#as35   if (FFTW3_OMP_LIBRARY AND FFTW3_INCLUDE)
#as35     add_imported_library(${FFTW3_OMP_LIBRARY} ${FFTW3_INCLUDE})
#as35   elseif (FFTW3_THREAD_LIBRARY AND FFTW3_INCLUDE)
#as35     add_imported_library(${FFTW3_THRED_LIBRARY} ${FFTW3_INCLUDE})
#as35   elseif (FFTW3L_LIBRARY AND FFTW3_INCLUDE)
#as35     add_imported_library(${FFTW3L_LIBRARY} ${FFTW3_INCLUDE})
#as35   elseif (FFTW3L_OMP_LIBRARY AND FFTW3_INCLUDE)
#as35     add_imported_library(${FFTW3L_OMP_LIBRARY} ${FFTW3_INCLUDE})
#as35   elseif (FFTW3L_THREAD_LIBRARY AND FFTW3_INCLUDE)
#as35     add_imported_library(${FFTW3L_THRED_LIBRARY} ${FFTW3_INCLUDE})
#as35   elseif (FFTW3F_LIBRARY AND FFTW3_INCLUDE)
#as35     add_imported_library(${FFTW3F_LIBRARY} ${FFTW3_INCLUDE})
#as35   elseif (FFTW3F_OMP_LIBRARY AND FFTW3_INCLUDE)
#as35     add_imported_library(${FFTW3F_OMP_LIBRARY} ${FFTW3_INCLUDE})
#as35   elseif (FFTW3F_THREAD_LIBRARY AND FFTW3_INCLUDE)
#as35     add_imported_library(${FFTW3F_THRED_LIBRARY} ${FFTW3_INCLUDE})
#as35   elseif (FFTW3Q_LIBRARY AND FFTW3_INCLUDE)
#as35     add_imported_library(${FFTW3Q_LIBRARY} ${FFTW3_INCLUDE})
#as35   elseif (FFTW3Q_OMP_LIBRARY AND FFTW3_INCLUDE)
#as35     add_imported_library(${FFTW3Q_OMP_LIBRARY} ${FFTW3_INCLUDE})
#as35   elseif (FFTW3Q_THREAD_LIBRARY AND FFTW3_INCLUDE)
#as35     add_imported_library(${FFTW3Q_THRED_LIBRARY} ${FFTW3_INCLUDE})
#as35   endif()
#as35   return()
#as35 endif()

find_path(FFTW3_INCLUDE NAMES fftw3.h
  HINTS "/usr/include" "/opt/local/include" "/app/include"
)

find_library(FFTW3_LIBRARY fftw3)
find_library(FFTW3_OMP_LIBRARY fftw3_omp)
find_library(FFTW3_THREAD_LIBRARY fftw3_threads)
find_library(FFTW3L_LIBRARY fftw3l)
find_library(FFTW3L_OMP_LIBRARY fftw3l_omp)
find_library(FFTW3L_THREAD_LIBRARY fftw3l_threads)
find_library(FFTW3F_LIBRARY fftw3f)
find_library(FFTW3F_OMP_LIBRARY fftw3f_omp)
find_library(FFTW3F_THREAD_LIBRARY fftw3f_threads)
find_library(FFTW3Q_LIBRARY fftw3q)
find_library(FFTW3Q_OMP_LIBRARY fftw3q_omp)
find_library(FFTW3Q_THREAD_LIBRARY fftw3q_threads)

# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

# FFTW3
find_package_handle_standard_args(FFTW3
                                  DEFAULT_MSG FFTW3_LIBRARY FFTW3_INCLUDE
                                 )
if (FFTW3_FOUND)
  add_imported_library("FFTW3" "${FFTW3_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

# FFTW3_OMP
find_package_handle_standard_args(FFTW3_OMP
                                  REQUIRED_VARS FFTW3_OMP_LIBRARY FFTW3_INCLUDE
                                  HANDLE_COMPONENTS
                                  NAME_MISMATCHED
                                 )
if (FFTW3_OMP_FOUND)
  add_imported_library("FFTW3_OMP" "${FFTW3_OMP_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

# FFTW3_THREAD
find_package_handle_standard_args(FFTW3_THREAD
                                  REQUIRED_VARS FFTW3_THREAD_LIBRARY FFTW3_INCLUDE
                                  HANDLE_COMPONENTS
                                  NAME_MISMATCHED
                                 )
if (FFTW3_THREAD_FOUND)
  add_imported_library("FFTW3_THREAD" "${FFTW3_THREAD_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

# FFTW3L
find_package_handle_standard_args(FFTW3L
                                  REQUIRED_VARS FFTW3L_LIBRARY FFTW3_INCLUDE
                                  HANDLE_COMPONENTS
                                  NAME_MISMATCHED
                                 )
if (FFTW3L_FOUND)
  add_imported_library("FFTW3L" "${FFTW3L_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

# FFTW3L_OMP
find_package_handle_standard_args(FFTW3L_OMP
                                  REQUIRED_VARS FFTW3L_OMP_LIBRARY FFTW3_INCLUDE
                                  HANDLE_COMPONENTS
                                  NAME_MISMATCHED
                                 )
if (FFTW3L_OMP_FOUND)
  add_imported_library("FFTW3L_OMP" "${FFTW3L_OMP_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

# FFTW3L_THREAD
find_package_handle_standard_args(FFTW3L_THREAD
                                  REQUIRED_VARS FFTW3L_THREAD_LIBRARY FFTW3_INCLUDE
                                  HANDLE_COMPONENTS
                                  NAME_MISMATCHED
                                 )
if (FFTW3L_THREAD_FOUND)
  add_imported_library("FFTW3L_THREAD" "${FFTW3L_THREAD_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

# FFTW3F
find_package_handle_standard_args(FFTW3F
                                  REQUIRED_VARS FFTW3F_LIBRARY FFTW3_INCLUDE
                                  HANDLE_COMPONENTS
                                  NAME_MISMATCHED
                                 )
if (FFTW3F_FOUND)
  add_imported_library("FFTW3F" "${FFTW3F_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

# FFTW3F_OMP
find_package_handle_standard_args(FFTW3F_OMP
                                  REQUIRED_VARS FFTW3F_OMP_LIBRARY FFTW3_INCLUDE
                                  HANDLE_COMPONENTS
                                  NAME_MISMATCHED
                                 )
if (FFTW3F_OMP_FOUND)
  add_imported_library("FFTW3F_OMP" "${FFTW3F_OMP_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

# FFTW3F_THREAD
find_package_handle_standard_args(FFTW3F_THREAD
                                  REQUIRED_VARS FFTW3F_THREAD_LIBRARY FFTW3_INCLUDE
                                  HANDLE_COMPONENTS
                                  NAME_MISMATCHED
                                 )
if (FFTW3F_THREAD_FOUND)
  add_imported_library("FFTW3F_THREAD" "${FFTW3F_THREAD_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

# FFTW3Q
find_package_handle_standard_args(FFTW3Q
                                  REQUIRED_VARS FFTW3Q_LIBRARY FFTW3_INCLUDE
                                  HANDLE_COMPONENTS
                                  NAME_MISMATCHED
                                 )
if (FFTW3Q_FOUND)
  add_imported_library("FFTW3Q" "${FFTW3Q_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

# FFTW3Q_OMP
find_package_handle_standard_args(FFTW3Q_OMP
                                  REQUIRED_VARS FFTW3Q_OMP_LIBRARY FFTW3_INCLUDE
                                  HANDLE_COMPONENTS
                                  NAME_MISMATCHED
                                 )
if (FFTW3Q_OMP_FOUND)
  add_imported_library("FFTW3Q_OMP" "${FFTW3Q_OMP_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

# FFTW3Q_THREAD
find_package_handle_standard_args(FFTW3Q_THREAD
                                  REQUIRED_VARS FFTW3Q_THREAD_LIBRARY FFTW3_INCLUDE
                                  HANDLE_COMPONENTS
                                  NAME_MISMATCHED
                                 )
if (FFTW3Q_THREAD_FOUND)
  add_imported_library("FFTW3Q_THREAD" "${FFTW3Q_THREAD_LIBRARY}" "${FFTW3_INCLUDE}")
endif()

pkg_check_modules(FFTW fftw3 QUIET)

if (FFTW3_FOUND)
  message("-- Found FFTW3: ${FFTW3_INCLUDE}, Version: ${FFTW_VERSION}")
endif (FFTW3_FOUND)
