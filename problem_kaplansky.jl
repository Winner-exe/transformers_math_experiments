include("constants.jl")

const M = 2
const N = 2

function convert_adjmat_to_string(adjmat::Matrix{Int})::String
    """
    Helper function to convert adjacency matrices into partitions (in string format).
    """
    return ""
end

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::Vector{OBJ_TYPE}
    """
    Main greedy search algorithm.
    It starts and ends with some construction

    Algorithm:

    Input: A string P corresponding to a partition of the edge set.

    (1) Refine the partition to get a PI struct.

    (2) Adjust partition to ensure orientable, no-fold, and no-pattern.

    (3) Find (an instance of) the shortest cycle, remove the least frequent edge by merging two 2-cells and splitting in the "opposite" way.

    (4) Do similar operations to ensure orientable, no-fold, and no-pattern are satisfied w/o adjusting the 2-cells from the prev. step (1).

    Potentially branch out and list many different ways of accomplishing (1) and (2).
    """
end

function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    """
    Function to calculate the reward of a final construction.

    In our case, calculates the min(girth(L_A), girth(L_B)) if input is a PI struct, or 0 if input is not a PI struct.
    """
    return 0
end


function empty_starting_point()::OBJ_TYPE
    """
    If there is no input file, the search starts always with this object
    (E.g. empty graph, all zeros matrix, etc)

    In our case, starts with trivial partition.
    """
    return ""
end