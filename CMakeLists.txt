# Ian Wark July 6 2017

cmake_minimum_required(VERSION 3.1)
project(HFold)

# only continue if we know we already have boost
# create simfold library
add_subdirectory(simfold)

set(SOURCE
    W_final.h W_final.cpp 
    pseudo_loop.h pseudo_loop.cpp 
    V_final.h V_final.cpp 
    VM_final.h VM_final.cpp
    h_common.h h_common.cpp 
    hfold_validation.h hfold_validation.cpp 
    hfold.h HFold.cpp
)

# pass SIMFOLD_HOME to binary
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# set path for executable to be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executable
add_executable(HFold ${SOURCE})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "CParty")

# enable C++11
target_compile_features(HFold PRIVATE cxx_nullptr)

# link to simfold
target_link_libraries(HFold LINK_PUBLIC simfold)
