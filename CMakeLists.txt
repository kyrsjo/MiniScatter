#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(MiniScatter)

#Unfortunately, the Geant4 headers doesn't really work with C++17
# Therefore, for now we stick to to 14 and load the libstd++fs
# manually (GNU only). Using a non-GNU compiler (e.g. Apple-CLANG)
# will loose some functionality.
set(CMAKE_CXX_STANDARD 14)

#CMAKE_POLICY(SET CMP0025 NEW) #Separate AppleClang from Clang

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Find ROOT variables if the variable G4ANALYSIS_USE_ROOT is set
#

set(CMAKE_MODULE_PATH ${Geant4_dir}/Modules)
#find_package(ROOT QUIET)
find_package(ROOT)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
if(ROOT_FOUND)
  message( STATUS "ROOT was found, all is well." )
  add_definitions(-DG4ANALYSIS_USE_ROOT)
  include(${Geant4_USE_FILE})
  include_directories(${ROOT_INCLUDE_DIRS} ${Geant4_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include)
  add_definitions(${ROOT_DEFINITIONS})
else()
  message( FATAL_ERROR "Root was not found, aborting!" )
  include(${Geant4_USE_FILE})
  include_directories(${PROJECT_SOURCE_DIR}/include)
endif()


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
list(FILTER sources EXCLUDE REGEX \\.sys\\.v)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
list(FILTER sources EXCLUDE REGEX \\.sys\\.v)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(MiniScatter MiniScatter.cc ${sources} ${headers})

target_compile_options(MiniScatter PRIVATE "-Wno-overloaded-virtual")

#message(${ROOT_FOUND})
if(ROOT_FOUND)
  target_link_libraries(MiniScatter ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
else()
  target_link_libraries(MiniScatter ${Geant4_LIBRARIES})
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build MiniScatter. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(MiniScatter_SCRIPTS
  vis.mac
  verbose.mac
  gammaFilter.mac
  miniScatterScanner.py
  miniScatterDriver.py
  miniScatterPlots.py
  quicksilver.py
)

foreach(_script ${MiniScatter_SCRIPTS})
#  configure_file(
#    ${PROJECT_SOURCE_DIR}/scripts/${_script}
#    ${PROJECT_BINARY_DIR}/${_script}
#    COPYONLY
#    )
  execute_process( COMMAND
    ${CMAKE_COMMAND} -E create_symlink
    ${PROJECT_SOURCE_DIR}/scripts/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(MiniScatter DEPENDS MiniScatter)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS MiniScatter DESTINATION bin)

## ROOT classes
add_subdirectory(ROOTclasses)
