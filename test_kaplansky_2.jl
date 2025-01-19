using DataStructures
using Graphs
using MetaGraphsNext
using Test
using Profile
using ProfileView

include("problem_kaplansky.jl")

function main()
    two_cells_encoded = ((1, 1), (2, 2)), ((1, 2), (2, 3)), ((2, 1), (3, 3)), ((4, 1), (1, 3)), ((3, 1), (5, 4)), ((3, 2), (4, 4)), ((1, 4), (6, 5)), ((2, 4), (7, 3)), ((3, 4), (6, 3)), ((3, 5), (4, 2)), ((4, 3), (7, 5)), ((5, 1), (1, 5)), ((5, 2), (7, 4)), ((5, 3), (8, 5)), ((8, 1), (2, 5)), ((8, 4), (4, 5)), ((5, 5), (7, 2)), ((6, 1), (8, 2)), ((6, 2), (8, 3)), ((7, 1), (6, 4))

    graph = convert_string_to_graph(empty_starting_point())
    for (edge1, edge2) in two_cells_encoded
       if edge1[1] < edge2[1]
           i1 = edge1[1]
           i2 = edge2[1]
           j1 = edge1[2]
           j2 = edge2[2]

           if can_add_two_cell(graph, i1, i2, j1, j2, 1)
               add_two_cell!(graph, i1, i2, j1, j2, 1)
           else
               print("Cannot add two-cell ($i1, $i2, $j1, $j2, 1)")
           end
       else
           i1 = edge2[1]
           i2 = edge1[1]
           j1 = edge2[2]
           j2 = edge1[2]

           if can_add_two_cell(graph, i1, i2, j1, j2, -1)
               add_two_cell!(graph, i1, i2, j1, j2, -1)
           else
               print("Cannot add two-cell ($i1, $i2, $j1, $j2, -1)")
           end
       end
    end

    empty_graph_str = empty_starting_point()
    # graph2_str = greedy_search_from_startpoint(0, empty_graph_str)[1]
    # println(graph2_str)
    # num_two_cells = count(x -> x != 0, parse.(Int, split(convert_graph_to_string(graph), ",")))
    # println(num_two_cells)
    # graph2 = convert_string_to_graph(graph2_str)
end

main()