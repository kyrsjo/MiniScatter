cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# Setup for the ROOT C++ classes used by the Python analysis
# Included from the top-level CMakeLists.txt

find_package(ROOT)
if(ROOT_FOUND)
    include_directories(${ROOT_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "ROOT was not found; aborting!")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

ROOT_GENERATE_DICTIONARY(MiniScatter_ROOTclasses_dict MiniScatter_ROOTclasses.h LINKDEF MiniScatter_ROOTclasses_linkdef.h)
set(MiniScatter_ROOTclasses_src MiniScatter_ROOTclasses.cpp MiniScatter_ROOTclasses_dict.cxx)
add_library(MiniScatter_ROOTclasses SHARED ${MiniScatter_ROOTclasses_src})
target_link_libraries(MiniScatter_ROOTclasses ${ROOT_LIBRARIES})

install(TARGETS MiniScatter_ROOTclasses DESTINATION bin)
