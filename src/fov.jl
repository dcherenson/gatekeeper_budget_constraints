using NearestNeighbors
using LinearAlgebra
using StaticArrays

function is_in_fov(x::SVector{3,F}, pt::SVector{2,F}, fovRadius::F, fovAngle::F) where {F}
    return norm(pt - x[SOneTo(2)]) < fovRadius && abs(angle_diff(atan(pt[2]-x[2], pt[1]-x[1]), x[3])) < fovAngle/2
end

function features_in_fov_idx(features::BallTree, x::SVector{3,F}, fovRadius::F, fovAngle::F) where {F}
    features_in_radius = Int16[]
    sizehint!(features_in_radius, length(features.data))
    inrange!(features_in_radius, features, x[SOneTo(2)], fovRadius)
    # check if features are within the fov fovAngle
    features_in_fov = Int16[]
    yaw_vec = @SVector[cos(x[3]), sin(x[3])]
    sizehint!(features_in_fov, length(features_in_radius))
    for i in features_in_radius
        pt = features.data[i] - x[SOneTo(2)]
        pt = pt/norm(pt)
        if pt⋅yaw_vec >= cos(fovAngle/2)
            push!(features_in_fov, i)
        end
    end

    return features_in_fov
end

