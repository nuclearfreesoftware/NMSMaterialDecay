#NMSMaterialDecay Library

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(NMSMaterialDecay)

include_directories(/geant/fission)

find_package(Geant4)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_library(NMSMaterialDecay ${sources} ${headers})

install (FILES ${headers} DESTINATION include)
install (TARGETS NMSMaterialDecay DESTINATION lib)

