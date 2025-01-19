using DataStructures
using Graphs
using MetaGraphsNext
using Test

include("problem_kaplansky.jl")

# Test converting functions
graph1_str = empty_starting_point()
graph1 = convert_string_to_graph(graph1_str)
println(@test graph1_str == convert_graph_to_string(graph1))

graph2_str = "1,0,0,1,1,0"
graph2 = convert_string_to_graph(graph2_str)
println(@test graph2_str == convert_graph_to_string(graph2))

graph3_str = greedy_search_from_startpoint(0, empty_starting_point())[1]
graph3 = convert_string_to_graph(graph3_str)

graph4_str = greedy_search_from_startpoint(0, graph2_str)[1]
graph4 = convert_string_to_graph(graph4_str)
