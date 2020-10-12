module Subdivision

import ForwardDiff
using LinearAlgebra

# Biquartic subdivision
# If the original mesh was a grid with edge length x:
# - the result is a grid with edge length x * 2^-k (where k is the # of iterations)
# - the result is shrunk compared to the original,
#   so there should be a border of width x * 1.5 around the important data

# Masks
mask_tl = [25 50 5; 50 100 10; 5 10 1] / 256
mask_tr = [5 50 25; 10 100 50; 1 10 5] / 256
mask_bl = [5 10 1; 50 100 10; 25 50 5] / 256
mask_br = [1 10 5; 10 100 50; 5 50 25] / 256

function subdivide(points, iterations)
    for _ in 1:iterations
        n, m = size(points)
        next = similar(points, 2(n-2), 2(m-2))
        for i in 2:n-1, j in 2:m-1
            next[2i-3,2j-3] = sum(points[i-1:i+1,j-1:j+1] .* mask_tl)
            next[2i-3,2j-2] = sum(points[i-1:i+1,j-1:j+1] .* mask_tr)
            next[2i-2,2j-3] = sum(points[i-1:i+1,j-1:j+1] .* mask_bl)
            next[2i-2,2j-2] = sum(points[i-1:i+1,j-1:j+1] .* mask_br)
        end
        points, next = next, points
    end
    points
end


# Test functions

function write_surface(points, filename)
    open(filename, "w") do f
        n, m = size(points)
        for i in 1:n, j in 1:m
            p = points[i,j]
            println(f, "v $(p[1]) $(p[2]) $(p[3])")
        end
        for i in 2:n, j in 2:m
            index = (i - 1) * m + j
            println(f, "f $(index) $(index-m) $(index-m-1) $(index-1)")
        end
    end
end

function create_data(f, bbox, resolution)
    points = Matrix{Vector{Float64}}(undef, resolution, resolution)
    for i in 1:resolution, j in 1:resolution
        u = (i - 1) / (resolution - 1)
        v = (j - 1) / (resolution - 1)
        x = bbox[1][1] * (1 - u) + bbox[2][1] * u
        y = bbox[1][2] * (1 - v) + bbox[2][2] * v
        points[i,j] = [x, y, f(x, y)]
    end
    points
end

function test(f, bbox, resolution, iterations, ground_truth_resolution = 50)
    @assert resolution > 4 "Resolution should be at least 5"
    l = max.(bbox[2] - bbox[1]) * (resolution - 1) / (resolution - 4)
    c = (bbox[1] + bbox[2]) / 2
    bbox = [c - l, c + l]
    gradient(x, y) = ForwardDiff.gradient(p -> f(p[1], p[2]), [x, y])
    function normalized(x, y)
        ndf = norm(gradient(x, y))
        if ndf > 0
            f(x, y) / ndf
        else
            f(x, y)
        end
    end
    points = create_data(normalized, bbox, resolution)
    grid = deepcopy(points)
    for p in grid
        p[3] = 0
    end
    write_surface(grid, "/tmp/grid.obj")
    points = subdivide(points, iterations)
    write_surface(points, "/tmp/test.obj")
    gtruth = create_data(f, bbox, ground_truth_resolution)
    write_surface(gtruth, "/tmp/gtruth.obj")
end

test1() = test((x,y) -> x^4+8x^3+y^4-16, ([-6,-2.5], [-1,2.5]), 20, 3)
test2() = test((x,y) -> x^2+y^2-1, ([-1,-1], [1,1]), 10, 3)
test3() = test((x,y) -> y^2-x^3+7x-6, ([-1.5,-5], [3,5]), 20, 3)

end
