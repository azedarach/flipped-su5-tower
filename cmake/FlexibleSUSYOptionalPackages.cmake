function(fs_package_status_message msg)
  message(STATUS "${Blue}${msg}${ColourReset}")
endfunction()

macro(find_optional_package name enable_name)
  if(${enable_name})
    find_package(${name})
    if(${name}_FOUND)
      fs_package_status_message("Enabling use of ${name}")
    else()
      unset(${enable_name})
      fs_package_status_message("Disabling use of ${name}")
    endif()
  else()
    fs_package_status_message("Disabling use of ${name}")
  endif()
endmacro()

# ==================== Debug mode ====================
function(find_debug)
  if(ENABLE_DEBUG)
    fs_package_status_message("Enabling debug output")
  endif()
endfunction()

# ==================== Threads ====================
function(find_threads)
  set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
  set(THREADS_PREFER_PTHREAD_FLAG TRUE)
  find_optional_package(Threads ENABLE_THREADS)
  if(NOT TARGET Threads::Threads)
    add_library(Threads::Threads INTERFACE IMPORTED)
  endif()
endfunction()

# ==================== Himalaya ====================
function(find_himalaya)
  if(ENABLE_HIMALAYA)
    find_package(Himalaya 1.0 REQUIRED)
    fs_package_status_message("Enabling use of Himalaya.")
  else()
    fs_package_status_message("Disabling use of Himalaya.")
  endif()
  if(NOT TARGET Himalaya::Himalaya)
    add_library(Himalaya::Himalaya INTERFACE IMPORTED)
  endif()
  if(NOT TARGET Himalaya::DSZ)
    add_library(Himalaya::DSZ INTERFACE IMPORTED)
  endif()
endfunction()

# ==================== LAPACK ====================
function(find_lapack)
  find_optional_package(LAPACK ENABLE_LAPACK)
  if(LAPACK_FOUND)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    set_property(TARGET LAPACK::LAPACK PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
  else()
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
  endif()
endfunction()

# ==================== LoopTools ====================
function(find_looptools)
  if(ENABLE_FFLITE AND ENABLE_LOOPTOOLS)
    message(FATAL_ERROR "Cannot use FFLite and LoopTools at the same time!")
  endif()
  if(ENABLE_LOOPTOOLS)
    message(WARNING
      "${Red}Warning: LoopTools is thread-unsafe since it accepts the renormalization scale\
       via a global variable.  This might lead to a race condition if different\
       RGE solvers run in multiple threads.  For such an application, consider the\
       alternative -DENABLE_FFLITE=ON instead.${ColourReset}")
    find_package(LoopTools REQUIRED)
    fs_package_status_message("Enabling use of LoopTools.")
  endif()
  if(NOT TARGET LoopTools::LoopTools)
    add_library(LoopTools::LoopTools INTERFACE IMPORTED)
  endif()
endfunction()

# ==================== Boost's odeint ====================
function(find_odeint)
  if(ENABLE_ODEINT)
    find_file(BOOST_FS_ODEINT_1 boost/numeric/odeint.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_2 boost/numeric/odeint/algebra/algebra_dispatcher.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_3 boost/numeric/odeint/algebra/vector_space_algebra.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_4 boost/numeric/odeint/external/eigen/eigen_algebra_dispatcher.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_5 boost/numeric/odeint/external/eigen/eigen_resize.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_6 boost/numeric/odeint/integrate/integrate_adaptive.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_7 boost/numeric/odeint/stepper/generation.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_8 boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp ${Boost_INCLUDE_DIRS})
    if(BOOST_FS_ODEINT_1 AND
        BOOST_FS_ODEINT_2 AND
        BOOST_FS_ODEINT_3 AND
        BOOST_FS_ODEINT_4 AND
        BOOST_FS_ODEINT_5 AND
        BOOST_FS_ODEINT_6 AND
        BOOST_FS_ODEINT_7 AND
        BOOST_FS_ODEINT_8)
      fs_package_status_message("Enabling use of Boost's odeint")
      mark_as_advanced(BOOST_FS_ODEINT_1 BOOST_FS_ODEINT_2
        BOOST_FS_ODEINT_3 BOOST_FS_ODEINT_4
        BOOST_FS_ODEINT_5 BOOST_FS_ODEINT_6
        BOOST_FS_ODEINT_7 BOOST_FS_ODEINT_8)
    else()
      unset(ENABLE_ODEINT)
      fs_package_status_message("Disabling use of Boost's odeint")
    endif()
  else()
    fs_package_status_message("Disabling use of Boost's odeint")
  endif()
endfunction()

# ==================== SQLite3 ====================
function(find_sqlite)
  if(ENABLE_SQLITE)
    find_package(SQLite3)
    if(SQLITE3_FOUND)
      add_library(SQLite3::SQLite3 INTERFACE IMPORTED)
      set_property(TARGET SQLite3::SQLite3 PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SQLITE3_INCLUDE_DIRS})
      set_property(TARGET SQLite3::SQLite3 PROPERTY INTERFACE_LINK_LIBRARIES ${SQLITE3_LIBRARIES})
      fs_package_status_message("Enabling use of SQLite3")
    endif()
  endif()
  if(NOT TARGET SQLite3::SQLite3)
    add_library(SQLite3::SQLite3 INTERFACE IMPORTED)
    fs_package_status_message("Disabling use of SQLite3")
  endif()
endfunction()
