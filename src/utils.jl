using Dubins
using StaticArrays

function wrapToPi(x)
    return atan(sin(x), cos(x))
end

function angle_diff(a, b)
    return wrapToPi(a - b)
end

function angle_add(a, b)
    return wrapToPi(a + b)
end

function make_path_from_waypoints(waypoints::Vector{SVector{3,F}}, r::F) where {F}
    # construct the best path
    path = DubinsPath[]
    sizehint!(path, length(waypoints)-1)
    for i=1:length(waypoints)-1
        e, p = dubins_shortest_path(waypoints[i], waypoints[i+1], r, 1e-3)
        @assert e == Dubins.EDUBOK
        push!(path, p)
    end
    return path
end

function sample_combined_dubins_path(path::Vector{DubinsPath}, s::F) where {F}
    # sample the combined path
    s_init = s
    for p in path
        L = dubins_path_length(p)
        if s < L
            e, q = dubins_path_sample(p, s)
            @assert e == Dubins.EDUBOK
            return q
        end
        s -= L
    end
    e,q = dubins_path_endpoint(path[end])
    @assert e == Dubins.EDUBOK
    return q
end

function dubins_sub_trajectory(path::Vector{DubinsPath}, s::F) where {F}
    # sample the combined path
    sub_path = DubinsPath[]
    for p in path
        L = dubins_path_length(p)
        if s < L
            e, q = dubins_extract_subpath(p, s)
            @assert e == Dubins.EDUBOK
            push!(sub_path, q)
            return sub_path
        else
            push!(sub_path, p)
            s -= L
        end
    end

    return sub_path
end

function dubins_traj_length(path::Vector{DubinsPath})
    L = 0.0
    for p in path
        L += dubins_path_length(p)
    end
    return L
end