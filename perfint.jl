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

function find_crossing(f, p1, p2)
    g(x) = f(p1 * (1 - x) + p2 * x)
    x = bisect(g, 0, 1, ε=1e-15)
    p1 * (1 - x) + p2 * x
end

function distances(f, bounding_box, resolution)
    delta = (bounding_box[2] - bounding_box[1]) ./ resolution
    dirs = [[delta[1], 0, 0],
            [0, delta[2], 0],
            [0, 0, delta[3]]]
    values = zeros(resolution .+ 1)
    for i in 0:resolution[1], j in 0:resolution[2], k in 0:resolution[3]
        values[i+1,j+1,k+1] = f(bounding_box[1] + delta .* [i, j, k])
    end

    write_volume(values, "/tmp/before")

    used_vertices = Dict()
    row, col, nz = Int[], Int[], Float64[]
    n_rows = 0

    function add_constraint!(i, j, k, di, dj, dk)
        v1, v2 = values[i,j,k], values[i+di,j+dj,k+dk]
        v1 * v2 > 0 && return
        ijk, dijk = [i, j, k], [di, dj, dk]
        p1 = bounding_box[1] + [delta[1] * (i - 1), delta[2] * (j - 1), delta[3] * (k - 1)]
        p2 = p1 + dirs[1] * di + dirs[2] * dj + dirs[3] * dk
        q = find_crossing(f, p1, p2)
        w1, w2 = norm((q - p2) .* dijk), norm((p1 - q) .* dijk)
        index1 = get!(used_vertices, ijk, length(used_vertices) + 1)
        index2 = get!(used_vertices, ijk .+ dijk, length(used_vertices) + 1)
        n_rows += 1
        push!(row, n_rows); push!(col, index1); push!(nz, w1)
        push!(row, n_rows); push!(col, index2); push!(nz, w2)
    end

    for i in 1:resolution[1]+1, j in 1:resolution[2]+1, k in 1:resolution[3]+1
        i <= resolution[1] && add_constraint!(i, j, k, 1, 0, 0)
        j <= resolution[2] && add_constraint!(i, j, k, 0, 1, 0)
        k <= resolution[3] && add_constraint!(i, j, k, 0, 0, 1)
    end

    # Add normalizing constraint
    n_rows += 1
    n_cols = length(used_vertices)
    for i in 1:n_cols
        push!(row, n_rows); push!(col, i); push!(nz, 1)
    end

    println("Sparse system: $n_rows x $n_cols with $(length(nz)) entries")

    # Setup sparse system
    A = sparse(row, col, nz, n_rows, n_cols)
    b = sparsevec(Dict(n_rows => 1), n_rows)
    A = Array(A)
    b = Array(b)
    x = A \ b

    # Change the corresponding values
    for (ijk, index) in used_vertices
        values[ijk...] = x[index]
    end

    write_volume(values, "/tmp/after")
    values
end

end # module
