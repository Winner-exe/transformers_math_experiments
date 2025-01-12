using DataStructures
using Graphs
using MetaGraphsNext
using StaticArrays

include("constants.jl")

const M = 2
const N = 3

function graph_error(error_msg::String, problem_cell=0::Int)::MetaGraph
    """
    Returns a MetaGraph indicating an error occurred with generating a colored oriented taiko.

    The error message stores the index of the problem 2-cell, if any.
    """
    return MetaGraph(
        Graph();
        label_type=Int,
        vertex_data_type=Int,
        edge_data_type=Int8,
        graph_data=error_msg
    )
end

function convert_string_to_graph(graph_str::String)::MetaGraph
    """
    Helper function to convert partitions (in string format) into digraphs.

    Graph is implemented as an undirected MetaGraph. Vertex labels are the same as in the paper (e.g. strings "a1" and "b1")
    and vertex data is given as a string ("a" or "b"). Edges are labelled as 2-tuples of strings with data values given as
    2-tuples of signed integers: the first gives the color and orientation, the second gives the number of 2-cells using this edge.

    If the given string results in a non-orientable or non-colorable taiko, returns an error graph so the local search will reject it later.
    Otherwise, returns the full colored and oriented taiko.
    """
    # Parse string and verify input length
    two_cells = parse.(Int, split(graph_str, ","))
    if size(two_cells, 1) != Int((M*N)*(M-1)*(N-1)/2)
        return graph_error("input error")
    elseif count(x -> x != 0, two_cells) > floor((M*N+1)/2)
        return graph_error("input error")
    end

    # Intialize vertex label-data pairs
    vertices_description = Array{Pair{String, String}, 1}(undef, M+N)
    for i in 1:M
        vertices_description[i] = "a$i" => "a"
    end

    for j in 1:N
        vertices_description[j+M] = "b$(j)" => "b"
    end

    # Initialize edge label-data pairs
    edges_description = Array{Pair{Tuple{String,String}, MVector{2, Int}}}(undef, M*N)
    for i in 1:M
        for j in 1:N
            edges_description[(i-1)*N + j] = ("a$i", "b$j") => MVector(0, 1)
        end
    end

    graph = MetaGraph(complete_bipartite_graph(M,N), vertices_description, edges_description, "coloring+orientation")

    # Add edges to form horizontal graphs L_A and L_B
    color = 1
    assign_color = 1
    index = 1

    # Check valid 2-cells
    for i1 in 1:M
        for j1 in 1:N
            for i2 in i1+1:M
                for j2 in 1:N
                    if j1 != j2
                        # 2-cell not in taiko
                        if two_cells[index] == 0
                            index += 1
                        # 2-cell is in taiko but violates partition condition
                        elseif graph["a$(i1)","b$(j1)"][1] != 0 || graph["a$(i2)","b$(j2)"][1] != 0
                            return graph_error("partition error at $(index)")
                        # 2-cell is in taiko, b graph orientation same
                        elseif j1 < j2
                            # Both horizontal edges filled in, check if orientation match and throw error if otherwise
                            if has_edge(graph, i1, i2) && has_edge(graph, j1 + M, j2 + M)
                                if sign(graph["a$(i1)","a$(i2)"][1]) != sign(graph["b$(j1)","b$(j2)"][1])
                                    return graph_error("orientation error at $(index)")
                                # Colors don't match, update rest of structure to use smaller color
                                elseif graph["a$(i1)","a$(i2)"][1] != graph["b$(j1)","b$(j2)"][1]
                                    new_color = min(abs(graph["a$(i1)","a$(i2)"][1]), abs(graph["b$(j1)","b$(j2)"][1]))
                                    old_color = max(abs(graph["a$(i1)","a$(i2)"][1]), abs(graph["b$(j1)","b$(j2)"][1]))

                                    for edge in edge_labels(graph)
                                        if abs(graph[edge[1], edge[2]][1]) == old_color
                                            sign_mod = sign(graph[edge[1], edge[2]][1])
                                            graph[edge[1], edge[2]][1] = new_color * sign_mod
                                        end
                                    end
                                    graph["a$(i1)","a$(i2)"][2] += 1
                                    graph["b$(j1)","b$(j2)"][2] += 1
                                    graph["a$(i1)","b$(j1)"][1] = new_color
                                    graph["a$(i2)","b$(j2)"][1] = new_color
                                else
                                    assign_color = abs(graph["a$(i1)","a$(i2)"][1])
                                    graph["a$(i1)","a$(i2)"][2] += 1
                                    graph["b$(j1)","b$(j2)"][2] += 1
                                    graph["a$(i1)","b$(j1)"][1] = assign_color
                                    graph["a$(i2)","b$(j2)"][1] = assign_color
                                    index += 1
                                end
                            # Exactly one horizontal edge filled in, make other horizontal edge match color.
                            # If orientation does not match with the given string, throw error.
                            elseif has_edge(graph, i1, i2)
                                if two_cells[index] != sign(graph["a$(i1)","a$(i2)"][1])
                                    return graph_error("orientation error at $(index)")
                                else
                                    assign_color = abs(graph["a$(i1)","a$(i2)"][1])
                                    graph["a$(i1)","a$(i2)"][2] += 1
                                    graph["b$(j1)","b$(j2)"] = MVector(graph["a$(i1)","a$(i2)"][1], 1)
                                    graph["a$(i1)","b$(j1)"][1] = assign_color
                                    graph["a$(i2)","b$(j2)"][1] = assign_color
                                    index += 1
                                end
                            elseif has_edge(graph, j1 + M, j2 + M)
                                if two_cells[index] != sign(graph["b$(j1)","b$(j2)"][1])
                                    return graph_error("orientation error at $(index)")
                                else
                                    assign_color = abs(graph["b$(j1)","b$(j2)"][1])
                                    graph["a$(i1)","a$(i2)"] = MVector(graph["b$(j1)","b$(j2)"][1], 1)
                                    graph["b$(j1)","b$(j2)"][2] += 1
                                    graph["a$(i1)","b$(j1)"][1] = assign_color
                                    graph["a$(i2)","b$(j2)"][1] = assign_color
                                    index += 1
                                end
                            # Neither edge filled in, update both with a new color.
                            else
                                assign_color = color
                                graph["a$(i1)","a$(i2)"] = MVector(assign_color * two_cells[index], 1)
                                graph["b$(j1)","b$(j2)"] = MVector(assign_color * two_cells[index], 1)
                                graph["a$(i1)","b$(j1)"][1] = assign_color
                                graph["a$(i2)","b$(j2)"][1] = assign_color
                                color += 1
                                index += 1
                            end
                        # 2-cell is in taiko, b graph orientation reversed
                        else
                            # Both horizontal edges filled in, check if color/orientation match and throw error if otherwise
                            if has_edge(graph, i1, i2) && has_edge(graph, j1 + M, j2 + M)
                                if sign(graph["a$(i1)","a$(i2)"][1]) != -1 * sign(graph["b$(j1)","b$(j2)"][1])
                                    return graph_error("orientation error at $(index)")
                                # Colors don't match, update rest of structure to use smaller color
                                elseif graph["a$(i1)","a$(i2)"][1] != -1 * graph["b$(j1)","b$(j2)"][1]
                                    new_color = min(abs(graph["a$(i1)","a$(i2)"][1]), abs(graph["b$(j1)","b$(j2)"][1]))
                                    old_color = max(abs(graph["a$(i1)","a$(i2)"][1]), abs(graph["b$(j1)","b$(j2)"][1]))

                                    for edge in edge_labels(graph)
                                        if abs(graph[edge[1], edge[2]][1]) == old_color
                                            sign_mod = sign(graph[edge[1], edge[2]][1])
                                            graph[edge[1], edge[2]][1] = sign_mod * new_color
                                        end
                                    end

                                    graph["a$(i1)","a$(i2)"][2] += 1
                                    graph["b$(j1)","b$(j2)"][2] += 1
                                    graph["a$(i1)","b$(j1)"][1] = new_color
                                    graph["a$(i2)","b$(j2)"][1] = new_color
                                else
                                    assign_color = abs(graph["a$(i1)","a$(i2)"][1])
                                    graph["a$(i1)","a$(i2)"][2] += 1
                                    graph["b$(j1)","b$(j2)"][2] += 1
                                    graph["a$(i1)","b$(j1)"][1] = assign_color
                                    graph["a$(i2)","b$(j2)"][1] = assign_color
                                    index += 1
                                end
                            # Exactly one horizontal edge filled in, make other horizontal edge match color.
                            # If orientation does not match with the given string, throw error.
                            elseif has_edge(graph, i1, i2)
                                if two_cells[index] != sign(graph["a$(i1)","a$(i2)"][1])
                                    return graph_error("orientation error at $(index)")
                                else
                                    assign_color = abs(graph["a$(i1)","a$(i2)"][1])
                                    graph["a$(i1)","a$(i2)"][2] += 1
                                    graph["b$(j1)","b$(j2)"] = MVector(-1 * graph["a$(i1)","a$(i2)"][1], 1)
                                    graph["a$(i1)","b$(j1)"][1] = assign_color
                                    graph["a$(i2)","b$(j2)"][1] = assign_color
                                    index += 1
                                end
                            elseif has_edge(graph, j1 + M, j2 + M)
                                if two_cells[index] != -1 * sign(graph["b$(j1)","b$(j2)"][1])
                                    return graph_error("orientation error at $(index)")
                                else
                                    assign_color = abs(graph["b$(j1)","b$(j2)"][1])
                                    graph["a$(i1)","a$(i2)"] = MVector(-1 * graph["b$(j1)","b$(j2)"][1], 1)
                                    graph["b$(j1)","b$(j2)"][2] += 1
                                    graph["a$(i1)","b$(j1)"][1] = assign_color
                                    graph["a$(i2)","b$(j2)"][1] = assign_color
                                    index += 1
                                end
                            # Neither edge filled in, update both with a new color.
                            else
                                assign_color = color
                                graph["a$(i1)","a$(i2)"] = MVector(assign_color * two_cells[index], 1)
                                graph["b$(j1)","b$(j2)"] = MVector(-1 * assign_color * two_cells[index], 1)
                                graph["a$(i1)","b$(j1)"][1] = assign_color
                                graph["a$(i2)","b$(j2)"][1] = assign_color
                                color += 1
                                index += 1
                            end
                        end
                    end
                end
            end
        end
    end
    return graph
end

function convert_graph_to_string(graph::MetaGraph)::String
    """
    Helper function to convert graphs into partitions (in string format).

    Partitions are represented as strings of (M*N)(M-1)*(N-1)/2 integers separated by commas. The
    string gives the upper triangular part of the adjacency matrix for the 2-skeleton of the
    taiko. The sign of the integer corresponds to the orientation of the A-part of the 2-cell (1 for positive,
    -1 for negative).

    Note that we can take a graph at any stage of the greedy search in here (even with unrefined partitions).

    However, any graphs that violate the partition or orientation conditions will be rejected here.
    """

    if occursin(graph[], "error")
        return "error"
    end

    entries = []
    for i1 in 1:M
        for j1 in 1:N
            for i2 in i1+1:M
                for j2 in 1:N
                    if j1 != j2
                        if has_edge(graph, i1, i2) && has_edge(graph, j1 + M, j2 + M) && abs(graph["a$(i1)","a$(i2)"][1]) == abs(graph["b$(j1)", "b$(j2)"][1])
                            if j1 < j2 && graph["a$(i1)","a$(i2)"][1] == graph["b$(j1)", "b$(j2)"][1]
                                push!(entries, "$(sign(graph["a$(i1)", "a$(i2)"][1])),")
                            elseif j2 < j1 && graph["a$(i1)","a$(i2)"][1] == -1 * graph["b$(j1)", "b$(j2)"][1]
                                push!(entries, "$(sign(graph["a$(i1)", "a$(i2)"][1])),")
                            else
                                push!(entries, "0,")
                            end
                        else
                            push!(entries, "0,")
                        end
                    end
                end
            end
        end
    end

    return chop(join(entries))
end

function is_valid_two_cell(graph::MetaGraph, i1::Int, i2::Int, j1::Int, j2::Int)::Bool
    """
    Returns true if valid oriented 2-cell, false otherwise.
    """

    # Degenerate square
    if i1 == i2 || j1 == j2
        return false
    # Bipartite edges aren't colored
    elseif graph["a$(i1)", "b$(j1)"][1] == 0 || graph["a$(i2)", "b$(j2)"][1] == 0
        return false
    # Bipartite edges don't match color
    elseif graph["a$(i1)", "b$(j1)"][1] != graph["a$(i2)", "b$(j2)"][1]
        return false
    # No horizontal edges
    elseif !(has_edge(graph, i1, i2) && has_edge(graph, j1+M, j2+M))
        return false
    end

    # Horizontal colors/orientation don't match
    if j1 < j2 && (graph["a$(i1)", "a$(i2)"][1] != graph["b$(j1)", "b$(j2)"][1])
        return false
    elseif j2 < j1 && (graph["a$(i1)", "a$(i2)"][1] != -1 * graph["b$(j1)", "b$(j2)"][1])
        return false
    end

    # Passed checks
    return true
end

function remove_two_cell!(graph::MetaGraph, i1::Int, i2::Int, j1::Int, j2::Int)
    """
    Helper function to remove a given 2-cell from the given structure.
    """

    # Check if valid 2-cell
    if !is_valid_two_cell(graph, i1, i2, j1, j2)
        throw("Not a valid 2-cell")
    end

    # Remove colors on bipartite edges
    graph["a$(i1)", "b$(j1)"][1] = 0
    graph["a$(i2)", "b$(j2)"][1] = 0

    # Remove the copies of the horizontal edges
    if graph["a$(i1)", "a$(i2)"][2] > 1
        graph["a$(i1)", "a$(i2)"][2] -= 1
    else
        rem_edge!(graph, i1, i2)
    end

    if graph["b$(j1)", "b$(j2)"][2] > 1
        graph["b$(j1)", "b$(j2)"][2] -= 1
    else
        rem_edge!(graph, j1 + M, j2 + M)
    end
end

function no_fold(graph::MetaGraph)::Bool
    """
    Returns true if the no_fold condition holds, false otherwise.
    """
    # To be implemented
    return true
end

function remove_folds!(graph::MetaGraph)
    """
    Forces the given graph to obey the no-fold condition by removing bad 2-cells.

    Algorithm: For each a-vertex, list all a-neighbors, detect folds (and how many), and delete 2-cells
    until folds are removed. Do the same for b-vertices. Deletion is done randomly.

    Note that positive integers correspond to outgoing edges and vice versa.
    """

    # Find and remove folds in L_A
    for i1 in 1:M
        used_colors = Set{Int}()
        fold_colors = Dict{Int, Int}()

        # Count all folds with corresponding colors/orientation
        for i2 in neighbors(graph, i1)
            if i2 > M
                continue
            elseif i1 < i2
                current_color = graph["a$(i1)", "a$(i2)"][1]
                if current_color in used_colors
                    if haskey(fold_colors, current_color)
                        fold_colors[current_color] += 1
                    else
                        fold_colors[current_color] = 1
                    end
                else
                    push!(used_colors, graph["a$(i1)", "a$(i2)"][1])
            else
                current_color = -1 * graph["a$(i1)", "a$(i2)"][1]
                if current_color in used_colors
                    if haskey(fold_colors, current_color)
                        fold_colors[current_color] += 1
                    else
                        fold_colors[current_color] = 1
                    end
                else
                    push!(used_colors, current_color)
            end
        end

        while !isempty(fold_colors)
            # Randomly choose a color to fold.
            color_choice = rand(collect(keys(fold_colors)))

            # Choose a-neighbor with the least copies to delete
            allowed_a_neighbors = Vector{Int}()
            for i2 in neighbors(graph, i1)
                if i2 > M
                    continue
                elseif i1 < i2 && graph["a$(i1)", "a$(i2)"][1] == color_choice
                    push!(allowed_a_neighbors, i2)
                elseif i2 < i1 && graph["a$(i1)", "a$(i2)"][1] == -1 * color_choice
                    push!(allowed_a_neighbors, i2)
                end
            end

            sort!(allowed_a_neighbors, by=x -> graph["a$(i1)", "a$(x)"][2])
            i2_choice = allowed_a_neighbors[1]

            # Randomly choose an associated two-cell
            allowed_two_cells = Vector{SVector{2, Int}}()
            for j1 in 1:N
                for j2 in 1:N
                    if is_valid_two_cell(i1, i2_choice, j1, j2)
                        push!(allowed_two_cells, SVector{Int}(j1, j2))
                    end
                end
            end
            two_cell_choice = rand(allowed_two_cells)
            j1_choice = two_cell_choice[1]
            j2_choice = two_cell_choice[2]

            remove_two_cell!(graph, i1, i2_choice, j1_choice, j2_choice)

            if !has_edge(graph, i1, i2)
                if fold_colors[color_choice] == 1
                    delete!(fold_colors, color_choice)
                else
                    fold_colors[color_choice] -= 1
                end
            end
        end
    end

    # Find and remove folds in L_B
    for j1 in 1:N
        used_colors = Set{Int}()
        fold_colors = Dict{Int, Int}()

        # Count all folds with corresponding colors/orientation
        for j2 in neighbors(graph, j1)
            if j2 <= M
                continue
            elseif j1 < j2 - M
                current_color = graph["b$(j1)", "b$(j2-M)"][1]
                if current_color in used_colors
                    if haskey(fold_colors, current_color)
                        fold_colors[current_color] += 1
                    else
                        fold_colors[current_color] = 1
                    end
                else
                    push!(used_colors, graph["b$(j1)", "b$(j2-M)"][1])
            else
                current_color = -1 * graph["b$(j1)", "b$(j2-M)"][1]
                if current_color in used_colors
                    if haskey(fold_colors, current_color)
                        fold_colors[current_color] += 1
                    else
                        fold_colors[current_color] = 1
                    end
                else
                    push!(used_colors, current_color)
            end
        end

        while !isempty(fold_colors)
            # Randomly choose a color to fold.
            color_choice = rand(collect(keys(fold_colors)))

            # Choose b-neighbor with the least copies to delete
            allowed_b_neighbors = Vector{Int}()
            for j2 in neighbors(graph, j1)
                if j2 <= M
                    continue
                elseif j1 < j2-M && graph["b$(j1)", "b$(j2-M)"][1] == color_choice
                    push!(allowed_b_neighbors, j2-M)
                elseif j2-M < j1 && graph["b$(j1)", "b$(j2-M)"][1] == -1 * color_choice
                    push!(allowed_b_neighbors, j2-M)
                end
            end

            sort!(allowed_b_neighbors, by=x -> graph["b$(j1)", "b$(x)"][2])
            j2_choice = allowed_b_neighbors[1]

            # Randomly choose an associated two-cell
            allowed_two_cells = Vector{SVector{2, Int}}()
            for i1 in 1:M
                for i2 in 1:M
                    if is_valid_two_cell(i1, i2, j1, j2_choice)
                        push!(allowed_two_cells, SVector{Int}(i1, i2))
                    end
                end
            end
            two_cell_choice = rand(allowed_two_cells)
            i1_choice = two_cell_choice[1]
            i2_choice = two_cell_choice[2]

            remove_two_cell!(graph, i1_choice, i2_choice, j1, j2_choice)

            if !has_edge(graph, j1+M, j2_choice+M)
                if fold_colors[color_choice] == 1
                    delete!(fold_colors, color_choice)
                else
                    fold_colors[color_choice] -= 1
                end
            end
        end
    end
end

function no_pattern(graph::MetaGraph)::Bool
    """
    Returns true if the no_pattern condition holds, false otherwise.
    """
    used_patterns = Set{Pair{Int, Int}}()

    # Check L_A
    for i1 in 1:M
        adj_edges = Vector{Pair{Int, Int}}()
        for i2 in neighbors(graph, i1)
            push!(adj_edges, i1 => i2)

        E = size(adj_edges, 1)
        for e1 in 1:E
            for e2 in e1+1:E
                color1 = 0
                color2 = 0

                i2_1 = adj_edges[e1][2]
                i2_2 = adj_edges[e2][2]

                if i1 < i2_1
                    color1 = graph["a$(i1)", "a$(i2_1)"][1]
                else
                    color1 = -1 * graph["a$(i1)", "a$(i2_1)"][1]
                end

                if i1 < i2_2
                    color2 = graph["a$(i1)", "a$(i2_2)"][1]
                else
                    color2 = -1 * graph["a$(i1)", "a$(i2_2)"][1]
                end

                if (color1 => color2) in used_patterns || (color2 => color1) in used_patterns
                    return false
                else
                    push!(used_patterns, color1 => color2)
                    push!(used_patterns, color2 => color1)
            end
        end
    end

    # Check L_B
    for j1 in 1:N
        adj_edges = Vector{Pair{Int, Int}}()
        for j2 in neighbors(graph, j1)
            push!(adj_edges, j1 => j2-M)

        E = size(adj_edges, 1)
        for e1 in 1:E
            for e2 in e1+1:E
                color1 = 0
                color2 = 0

                j2_1 = adj_edges[e1][2]
                j2_2 = adj_edges[e2][2]

                if j1 < j2_1
                    color1 = graph["b$(j1)", "b$(j2_1)"][1]
                else
                    color1 = -1 * graph["b$(j1)", "b$(j2_1)"][1]
                end

                if j1 < j2_2
                    color2 = graph["b$(j1)", "b$(j2_2)"][1]
                else
                    color2 = -1 * graph["b$(j1)", "b$(j2_2)"][1]
                end

                if (color1 => color2) in used_patterns || (color2 => color1) in used_patterns
                    return false
                else
                    push!(used_patterns, color1 => color2)
                    push!(used_patterns, color2 => color1)
            end
        end
    end

    return true
end

function remove_patterns!(graph::MetaGraph)
    """
    Forces the given graph to obey the no-pattern condition by removing bad 2-cells.

    Algorithm: For each repeated pattern, list all horizontal edges realizing the pattern. Delete 2-cells corresponding
    to the edge of least occurrence until pattern no longer repeats.

    Note that positive integers correspond to outgoing edges and vice versa.
    """

    # To be implemented
end

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::Vector{OBJ_TYPE}
    """
    Main greedy search algorithm.
    It starts and ends with some construction

    Algorithm (taiko construction taken from https://mineyev.web.illinois.edu/art/top-geom-uzd-origami.pdf):

    Input: A string P corresponding to a partition of the edge set.

    (1) Check if current partition is orientable, no-fold, and no-pattern. If not, remove 2-cells until it is.

    Kinds of bad things:
    1. Partition error - remove bad 2-cells and try again
    2. Orientation error - try flipping orientation of bad 2-cell, otherwise remove it

    Note that partition and orientation errors will halt further graph construction, so modifications are made at the string level.

    3. Fold - remove bad 2-cell
    4. Pattern - remove bad 2-cell

    No-fold and no-pattern are checked after graph construction, so modifications are made at the graph level.
    Note that removing a 2-cell will never make it violate partition or orientation conditions any more than it already does.

    (2) Add allowed 2-cells randomly until structure is a taiko. Ensure after each step that the structure is still orientable, no-fold, and no-pattern.
    (i.e. If adding a 2-cell violates any of these conditions, undo. If there are no more moves left, backtrack. (use a memory stack for this))

    (3) Find (an instance of) the shortest cycle, remove the least frequent edge by merging two 2-cells and splitting in the "opposite" way.

    (4) Do similar operations to ensure orientable, no-fold, and no-pattern are satisfied.

    Potentially branch out and list many different ways of accomplishing (1) and (2).
    """
    graph = convert_string_to_graph(obj)

    # Reject any obviously wrong input (wrong length, too many 2-cells)
    if graph[] == "input error"
        return []
    end

    # Fix partition
    while occursin(graph[], "partition error")
        two_cells = split(obj, ",")
        index = parse(Int, last(graph[], 1))
        two_cells[index] = "0"

        graph = convert_string_to_graph(join(two_cells, ","))
    end

    # Fix orientation
    while occursin(graph[], "orientation error")
        two_cells = split(obj, ",")
        index = parse(Int, last(graph[], 1))
        two_cells[index] = string(-1 * parse(Int, two_cells[index]))

        graph = convert_string_to_graph(join(two_cells, ","))

        if occursin(graph[], "orientation error")
            two_cells = split(obj, ",")
            index = parse(Int, last(graph[], 1))
            two_cells[index] = "0"

            graph = convert_string_to_graph(join(two_cells, ","))
        end
    end

    # Remove folds and patterns
    remove_folds!(graph)
    remove_patterns!(graph)

    # To be implemented: (2) adding 2-cells in with memory stack

    # To be implemented: (3) and (4) greedily increasing girth of L_AB
    return [convert_graph_to_string(graph)]
end

function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    """
    Function to calculate the reward of a final construction.

    In our case, calculates the min(girth(L_A), girth(L_B)) if input is a PI struct, or 0 if input is not a PI struct.
    Note that we can bound girth calculations to girths of at most 6 since we win if we reach girth 6.
    """
    return 0
end


function empty_starting_point()::OBJ_TYPE
    """
    If there is no input file, the search starts always with this object
    (E.g. empty graph, all zeros matrix, etc)

    In our case, starts with trivial coloring.
    """
    return chop("0," ^ Int((M*N)*(M-1)*(N-1)/2))
end