SET(flowcutter_SOURCES 
	src/cell.cpp 
	src/greedy_order.cpp
	src/list_graph.cpp
	src/tree_decomposition.cpp)

add_library(flowcutter ${flowcutter_SOURCES})
