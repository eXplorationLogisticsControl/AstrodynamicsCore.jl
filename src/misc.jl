"""Misc helper functions"""

"""
    get_sphere_mesh(r=0.5, n=30)

Build a unit-scale sphere mesh for visualization (GeometryBasics `Mesh`).
"""
function get_sphere_mesh(r::Real=0.5, n::Int=30)
    # Create vertices for a Sphere
    θ = LinRange(0, pi, n)
    φ2 = LinRange(0, 2pi, 2 * n)
    x2 = [r * cos(φv) * sin(θv) for θv in θ, φv in φ2]
    y2 = [r * sin(φv) * sin(θv) for θv in θ, φv in φ2]
    z2 = [r * cos(θv) for θv in θ, φv in 2φ2]
    points = vec([Point3f(xv, yv, zv) for (xv, yv, zv) in zip(x2, y2, z2)])
    
    # The coordinates form a matrix, so to connect neighboring vertices with a face
    # we can just use the faces of a rectangle with the same dimension as the matrix:
    _faces = decompose(QuadFace{GLIndex}, Tessellation(Rect(0, 0, 1, 1), size(z2)))
    # Normals of a centered sphere are easy, they're just the vertices normalized.
    _normals = normalize.(points)
    
    # Now we generate UV coordinates, which map the image (texture) to the vertices.
    # (0, 0) means lower left edge of the image, while (1, 1) means upper right corner.
    function gen_uv(shift)
        return vec(map(CartesianIndices(size(z2))) do ci
            tup = ((ci[1], ci[2]) .- 1) ./ ((size(z2) .* shift) .- 1)
            return Vec2f(reverse(tup))
        end)
    end
    return GeometryBasics.Mesh(points, _faces)
end



"""Trailing window mean of length `n`; output has `length(vs) - n + 1` samples."""
function moving_average(vs::AbstractVector{<:Real}, n::Int)
    n >= 1 || throw(ArgumentError("window length n must be at least 1"))
    m = length(vs) - n + 1
    m < 1 && return [sum(vs) / length(vs)]
    return [sum(@view vs[i:(i + n - 1)]) / n for i in 1:m]
end

"""`moving_average` padded at the front so output length matches `vs`."""
function _moving_average_padded(vs::AbstractVector{<:Real}, n::Int)
    avg = moving_average(vs, n)
    return vcat(fill(avg[1], length(vs) - length(avg)), avg)
end

_wrap_to_pi(x) = mod(x + π, 2π) - π