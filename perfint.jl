module PerfInt

using LinearAlgebra
using SparseArrays

"""
    write_volume(voxels, filename)

Writes the 3D matrix `voxels` into `filename.raw`
with the metaimage header `filename.mhd`.
"""
function write_volume(voxels, filename)
    n = size(voxels)
    open("$(filename).mhd", "w") do f
        println(f, "NDims = 3")
        println(f, "DimSize = $(n[1]) $(n[2]) $(n[3])")
        println(f, "ElementSize = 1.0 1.0 1.0")
        println(f, "ElementSpacing = 1.0 1.0 1.0")
        println(f, "ElementType = MET_DOUBLE")
        println(f, "ElementByteOrderMSB = False") # LittleEndian
        println(f, "ElementDataFile = $(filename).raw")
    end
    open("$(filename).raw", "w") do f
        for val in voxels
            write(f, val)
        end
    end
end

"""
    bisect(f, xl, xh; iterations, ɛ)

`xl` and `xh` bracket the root of f.
"""
function bisect(f, xl, xh; iterations = 100, ɛ = 1.0e-7)
    xm = xl
    for i = 1:iterations
        xm_old = xm
        xm = (xl + xh) / 2
        if xm != 0 && abs((xm - xm_old) / xm) < ɛ
            break
        end
        test = f(xl) * f(xm)
        if test < 0
            xh = xm
        elseif test > 0
            xl = xm
        else
            break
        end
    end
    xm
end

"""
    find_crossing(f, p1, p2)

Finds the exact position between points `p1` and `p2` where `f` vanishes.
"""
function find_crossing(f, p1, p2)
    g(x) = f(p1 * (1 - x) + p2 * x)
    x = bisect(g, 0, 1, ε=1e-15)
    p1 * (1 - x) + p2 * x
end

# Example:
# PerfInt.distances(p -> norm(p) - 1, ([-1.2,-1.2,-1.2], [2,2,2]), (10,10,10));
function distances(f, bounding_box, resolution; debug_output = false)
    delta = (bounding_box[2] - bounding_box[1]) ./ resolution
    dirs = [[delta[1], 0, 0],
            [0, delta[2], 0],
            [0, 0, delta[3]]]
    values = zeros(resolution .+ 1)
    for i in 0:resolution[1], j in 0:resolution[2], k in 0:resolution[3]
        values[i+1,j+1,k+1] = f(bounding_box[1] + delta .* [i, j, k])
    end

    function deviation_stats()    # just for debug output
        edges = [0, 4, 1, 5, 2, 6, 3, 7,
                 0, 2, 1, 3, 4, 6, 5, 7,
                 0, 1, 2, 3, 4, 5, 6, 7]
        max_dev, avg_dev, count = 0, 0, 0
        for i in 1:resolution[1], j in 1:resolution[2], k in 1:resolution[3]
            vertices = []
            for di in 0:1, dj in 0:1, dk in 0:1
                push!(vertices, values[i+di,j+dj,k+dk])
            end
            origin = bounding_box[1] + [delta[1] * (i - 1), delta[2] * (j - 1), delta[3] * (k - 1)]
            for i in 1:12
                v1, v2 = edges[2i-1], edges[2i]
                a, b = vertices[v1+1], vertices[v2+1]
                a * b > 0 && continue
                denom = abs(b - a)
                x = denom < 1e-15 ? 0.5 : abs(a) / denom
                c = (i - 1) ÷ 4 + 1
                q = dirs[1] * (v1 ÷ 4) + dirs[2] * (v1 % 4 ÷ 2) + dirs[3] * (v1 % 2) + dirs[c] * x
                dev = abs(f(origin + q))
                max_dev = max(max_dev, dev)
                avg_dev += dev
                count += 1
            end
        end
        avg_dev /= count
        "max -> $max_dev, avg -> $avg_dev"
    end

    debug_output && println("Deviation (before): $(deviation_stats())")
    debug_output && write_volume(values, "/tmp/before")

    used_vertices = Dict()
    row, col, nz = Int[], Int[], Float64[]
    n_rows = 0

    function add_constraint!(ijk, c)
        dijk = [x == c ? 1 : 0 for x in 1:3]
        v1, v2 = values[ijk...], values[(ijk+dijk)...]
        v1 * v2 > 0 && return
        p1 = bounding_box[1] + delta .* (ijk .- 1)
        p2 = p1 + sum(dirs .* dijk)
        q = find_crossing(f, p1, p2)
        w1, w2 = (p2 - q)[c], (p1 - q)[c]
        index1 = get!(used_vertices, ijk, length(used_vertices) + 1)
        index2 = get!(used_vertices, ijk + dijk, length(used_vertices) + 1)
        n_rows += 1
        push!(row, n_rows); push!(col, index1); push!(nz, w1)
        push!(row, n_rows); push!(col, index2); push!(nz, w2)
    end

    for i in 1:resolution[1]+1, j in 1:resolution[2]+1, k in 1:resolution[3]+1
        i <= resolution[1] && add_constraint!([i, j, k], 1)
        j <= resolution[2] && add_constraint!([i, j, k], 2)
        k <= resolution[3] && add_constraint!([i, j, k], 3)
    end

    # Add normalizing constraint
    n_rows += 1
    n_cols = length(used_vertices)
    for i in 1:n_cols
        push!(row, n_rows); push!(col, i); push!(nz, 1)
    end

    debug_output && println("Sparse system: $n_rows x $n_cols with $(length(nz)) entries")

    # Setup sparse system
    A = sparse(row, col, nz, n_rows, n_cols)
    b = sparsevec(Dict(n_rows => 1), n_rows)
    A = Array(A)
    b = Array(b)
    x = A \ b

    # Change the corresponding values
    for (ijk, index) in used_vertices
        values[ijk...] = copysign(x[index], values[ijk...])
    end

    debug_output && println("Deviation (after): $(deviation_stats())")
    debug_output && write_volume(values, "/tmp/after")

    values
end

end # module
