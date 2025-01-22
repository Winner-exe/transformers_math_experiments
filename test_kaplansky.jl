using DataStructures
using Graphs
using MetaGraphsNext
using Test

include("problem_kaplansky_num_two_cells.jl")

# Test converting functions
graph1_str = empty_starting_point()
graph1 = convert_string_to_graph(graph1_str)
println(@test graph1_str == convert_graph_to_string(graph1))

#graph2_str = greedy_search_from_startpoint(0, empty_starting_point())[1]
#graph2 = convert_string_to_graph(graph2_str)

function generate_random_taiko()
   graph_str = greedy_search_from_startpoint(0, empty_starting_point())[1]
   graph = convert_string_to_graph(graph_str)
   println(graph)
   println(graph_str)

   num_two_cells = count(x -> x != 0, parse.(Int, split(graph_str, ",")))
   println(num_two_cells)
   println(reward_calc(graph_str))
end
