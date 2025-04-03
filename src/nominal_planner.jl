using .RRTStar
using .World
using StaticArrays
using Dubins
using NearestNeighbors

# get just the distance between two nodes
function dubins_distance(q1, q2, turning_radius=0.1)

    e, p = Dubins.dubins_shortest_path(q1, q2, turning_radius, 1e-3)
    @assert e == Dubins.EDUBOK

    return Dubins.dubins_path_length(p)
end

@kwdef struct NominalPlannerProblem{F, VW} <: RRTStar.AbstractProblem{SVector{3,F}}
    domain::Tuple{SVector{3,F}, SVector{3,F}} # rectangle defined by opposite corners
    turning_radius::F = 10.0
    max_velocity::F = 10.0
    near_radius::F = 30.0
    unsafe_zones::VW
    destination_indices::Vector{Int64} = Int64[]
end

struct NominalPlanner{TP, TN}
    prob::TP # NominalPlannerProblem{Float64}
    nodes::TN # Vector{RRTStar.Node{SVector{3,Float64}}}
end

function initialize_nominal_planner(goal::TV, domain::Tuple{TV, TV}, zones::Vector{W}, max_iter) where {TV, W <: World.AbstractUnsafeZone}
    P = NominalPlannerProblem(domain = domain, unsafe_zones = zones)
    nodes = Vector{RRTStar.Node{TV}}()
    flipped_yaw = angle_add(goal[3], pi) # flip yaw because we are planning in reverse
    flipped_goal = @SVector [goal[1], goal[2], flipped_yaw]
    push!(nodes, RRTStar.Node(flipped_goal))

    RRTStar.rrt_star!(P, nodes, max_iter)

    return NominalPlanner(P, nodes)
end

function query_nominal_planner(planner::NominalPlanner, start)
    start_flipped = @SVector [start[1], start[2], angle_add(start[3], pi)]
    dist, path = RRTStar.get_best_path(planner.prob, planner.nodes, start_flipped, rev=false)
    
    if !isfinite(dist)
        return false, path
    end
    
    for i = 1:length(path)
        path[i] = @SVector [path[i][1], path[i][2], angle_add(path[i][3], pi)]
    end
    return true, path

end

function sample_domain(P::NominalPlannerProblem)
    v = @SVector rand(3)
    return P.domain[1] + (P.domain[2] - P.domain[1]) .* v
end

function RRTStar.sample_free(problem::NominalPlannerProblem)
    while true
        q = sample_domain(problem)
        # check if q is unsafe
        unsafe = World.is_unsafe(problem.unsafe_zones, q)
        if !unsafe
            return q
        end
    end
end

function RRTStar.nearest(P::NominalPlannerProblem, nodes, x_rand)
    best_dist = Inf
    best_ind = -1

    for i=1:length(nodes)
        if (nodes[i].state[1] - x_rand[1])^2 + (nodes[i].state[2] - x_rand[2])^2 > best_dist^2
            continue
        end
        d = dubins_distance(nodes[i].state, x_rand, P.turning_radius)
        if d < best_dist
            best_dist = d
            best_ind = i
        end
    end
    return best_ind
end

function RRTStar.near(P::NominalPlannerProblem, nodes, x_new)
    ind = Int[]
    # sizehint!(ind, length(nodes)/2)
    for i=1:length(nodes)
        if (nodes[i].state[1] - x_new[1])^2 + (nodes[i].state[2] - x_new[2])^2 > P.near_radius^2
            continue
        end
        d = dubins_distance(nodes[i].state, x_new, P.turning_radius)
        if d < P.near_radius
            push!(ind, i)
        end
    end
    return ind
end

function RRTStar.steer(problem::NominalPlannerProblem, x_nearest, x_rand; max_travel_dist = 5.1)

    # first plan the full dubins path
    errcode, path = Dubins.dubins_shortest_path(x_nearest, x_rand, problem.turning_radius, 1e-3)
    @assert errcode == Dubins.EDUBOK

    L = Dubins.dubins_path_length(path) 

    # get the node 20% of the waythrough
    s = min(L, max_travel_dist)

    errcode, x_new = Dubins.dubins_path_sample(path, s)
    @assert errcode == Dubins.EDUBOK
    
    return SVector{3}(x_new)
    
end

function RRTStar.collision_free(problem::NominalPlannerProblem, x_nearest, x_new; step_size=0.1)

    # create the path
    errcode, path = Dubins.dubins_shortest_path(x_nearest, x_new, problem.turning_radius, 1e-3)
    @assert errcode == Dubins.EDUBOK

    L = Dubins.dubins_path_length(path)

    x = 0.0
    while x < L
        errcode, q = Dubins.dubins_path_sample(path, x)

        @assert errcode == Dubins.EDUBOK errcode
    
        # check each point for collision 
        if q[1] < problem.domain[1][1] || q[1] > problem.domain[2][1] || q[2] < problem.domain[1][2] || q[2] > problem.domain[2][2] || World.is_unsafe(problem.unsafe_zones, q)
            return false
        end
        x += step_size
    end

    return true
end

function RRTStar.path_cost(problem::NominalPlannerProblem, x_near, x_new)

    # create the path
    errcode, path = Dubins.dubins_shortest_path(x_near, x_new, problem.turning_radius, 1e-3)
    @assert errcode == Dubins.EDUBOK    

    L = Dubins.dubins_path_length(path)

    return L
end

function invalidate_nodes!(problem::NominalPlannerProblem, nodes::Vector{RRTStar.Node{TN}}) where {TN}
    for i = 1:length(nodes)
        if nodes[i].parent_index > 0 && !RRTStar.collision_free(problem,  nodes[nodes[i].parent_index].state,nodes[i].state)
            nodes[i] = RRTStar.Node{TN}(nodes[i].state, -1, Inf)
        end
    end
end