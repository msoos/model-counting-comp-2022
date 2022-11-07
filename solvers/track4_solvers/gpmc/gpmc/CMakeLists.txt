cmake_minimum_required(VERSION 2.8)
project(mc2021 CXX)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/gmp-lib/lib/ ${CMAKE_CURRENT_SOURCE_DIR}/mpfr-lib/lib/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/flow-cutter-pace17/src /usr/local/include)

add_library(gpmc STATIC 
	core/ComponentCache.cc core/ComponentManager.cc core/Counter.cc
	core/Solver.cc 
	utils/Options.cc utils/System.cc
	preprocessor/Instance.cc preprocessor/Preprocessor.cc preprocessor/TestSolver.cc 
	preprocessor/lib_sharpsat_td/subsumer.cpp
	preprocessor/TreeDecomposition.cc preprocessor/IFlowCutter.cc
	)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -D NDEBUG -Wall -Wno-parentheses -Wno-format -Wno-sign-compare -Wno-literal-suffix -Wno-unused-variable -Wno-unused-result")
add_definitions(-D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS)

add_executable(gpmc-bin core/Main.cc)
set_target_properties(gpmc-bin PROPERTIES OUTPUT_NAME gpmc)

add_subdirectory(flow-cutter-pace17)
target_link_libraries(gpmc-bin gpmc flowcutter z mpfr gmp gmpxx)
set_target_properties(gpmc-bin PROPERTIES LINK_FLAGS "-static")
