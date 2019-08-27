module Marching2D
using LinearAlgebra

# Global settings

filename = "/tmp/test2d.eps"
sampling_res = 100

### SETUP_BEGIN - do not modify this line
corner = [-1.4, -1.55]
bbox_edge = 3.0
cells = 4
ellipse = [1.2, 0.8]
curve_original = p -> (p[1]/ellipse[1])^2 + (p[2]/ellipse[2])^2 - 1
gradient_original = p -> 2 * [p[1]/ellipse[1]^2, p[2]/ellipse[2]^2]
real_curve = [[ellipse[1]*cos(t), ellipse[2]*sin(t)] for t in 0:0.01:2pi]
show_types = [:real :nolinear :liming_cubic]
### SETUP_END - do not modify this line

# Global constants

fullwidth = 2^-0.25 / 4 * 100 / 2.54 * 72 # A4 = A0 / 4, from centimeters to inches to points
fullheight = fullwidth * sqrt(2)
margin = 0.25 * 72
width = fullwidth - 2 * margin
height = fullheight - 2 * margin
colors = [ 1 0 0 1 1 0 0
           0 1 0 1 0 1 0
           0 0 1 0 1 1 0 ]
point_radius = 3

curve = p -> curve_original(p) / norm(gradient_original(p))
gradient = p -> normalize(gradient_original(p))

# Global variables

global accuracy
check_accuracy = false
global file_handle              # for debug


# Main code

function print_header(f)
    println(f, """%!PS-Adobe-2.0
                  %%BoundingBox: 0 0 $fullwidth $fullheight
                  %%EndComments
                  1 setlinecap
                  1 setlinejoin""")
end

function convert(p)
    q = (p - corner) / bbox_edge * width .+ margin
    q[2] = fullheight - q[2]
    q
end

function print_grid(f)
    println(f, """[3 5] 0 setdash
                  0.5 setgray""")
    for i in 0:cells
        println(f, """newpath
                      $(convert(corner)[1] + i / cells * width) $(height + margin) moveto
                      0 $(-width) rlineto
                      stroke""")
        println(f, """newpath
                      $margin $(convert(corner)[2] - i / cells * width) moveto
                      $width 0 rlineto
                      stroke""")
    end
    println(f, """[] 0 setdash
                  0 setgray""")
end

function print_segments(f, segments)
    isempty(segments) && return
    p = convert(segments[1])
    println(f, """newpath
                  $(p[1]) $(p[2]) moveto""")
    for q in segments[2:end]
        p = convert(q)
        println(f, "$(p[1]) $(p[2]) lineto")
    end
    println(f, "stroke")
end

function print_point(f, p)
    p = convert(p)
    println(f, "$(p[1]) $(p[2]) $point_radius 0 360 arc fill")
end

function find_intersection(p1, p2)
    v1 = curve(p1)
    v2 = curve(p2)
    v1 * v2 > 0 && return nothing
    x = abs(v1) / abs(v2 - v1)
    p1 * (1 - x) + p2 * x
end

function solve_quadratic(a, b, c)
    eps = 1.0e-8
    if abs(a) < eps
        if abs(b) < eps
            @assert abs(c) < eps "No solutions"
            return 0
        end
        return -c / b
    end
    D = b^2 - 4*a*c
    if D < -eps
        @warn "D = $D, using 0 discriminant"
    end
    D = D < 0 ? 0 : sqrt(D)
    ((-b + D) / 2a, (-b - D) / 2a)
end

function liming_parabola(f1, g1, f2, g2)
    # Assuming equation S3^2 - lambda S1 S2 = 0
    # Using vectors as polynomials
    # - linear: [x, y, 1]
    # - quadratic: [x^2, xy, y^2, x, y, 1]
    d = f2 - f1
    n = normalize([-d[2], d[1]])
    if dot(g1, g2) < 0
        @warn "Opposing gradients"
    end
    if dot(n, g1) < 0
        n *= -1
    end
    S1 = [g1[1], g1[2], -dot(f1, g1)]
    S2 = [g2[1], g2[2], -dot(f2, g2)]
    S3 = [n[1], n[2], -dot(f1, n)]
    S1S2 = [S1[1]*S2[1], S1[1]*S2[2]+S1[2]*S2[1], S1[2]*S2[2],
            S1[1]*S2[3]+S1[3]*S2[1], S1[2]*S2[3]+S1[3]*S2[2], S1[3]*S2[3]]
    S3S3 = [S3[1]^2, 2*S3[1]*S3[2], S3[2]^2, 2*S3[1]*S3[3], 2*S3[2]*S3[3], S3[3]^2]

    # Compute correct lambda for parabola, i.e., b^2 - 4ac = 0
    # This is a quadratic equation in lambda, but the constant term vanishes
    a = [S1S2[1], S3S3[1]]
    c = [S1S2[3], S3S3[3]]
    b = [S1S2[2], S3S3[2]]
    lambda = (2*b[1]*b[2]-4*(a[1]*c[2]+a[2]*c[1])) / (b[1]^2-4*a[1]*c[1])

    S3S3 - lambda * S1S2
end

function intersect_implicit(S, p1, p2)
    if p1[1] == p2[1]
        x = p1[1]
        y1, y2 = solve_quadratic(S[3], S[2] * x + S[5], S[1] * x^2 + S[4] * x + S[6])
        lo, hi = min(p1[2], p2[2]), max(p1[2], p2[2])
        y = y1 < lo || y1 > hi ? y2 : y1
        if y < lo || y > hi
            @warn "$y is not in [$lo, $hi]"
        end
        [x, y]
    else
        @assert p1[2] == p2[2] "The points should share exactly one coordinate"
        y = p1[2]
        x1, x2 = solve_quadratic(S[1], S[2] * y + S[4], S[3] * y^2 + S[5] * y + S[6])
        lo, hi = min(p1[1], p2[1]), max(p1[1], p2[1])
        x = x1 < lo || x1 > hi ? x2 : x1
        if x < lo || x > hi
            @warn "$x is not in [$lo, $hi]"
        end
        [x, y]
    end
end

function find_intersection2(p1, p2)
    v1 = curve(p1)
    v2 = curve(p2)
    v1 * v2 > 0 && return nothing
    g1 = gradient(p1)
    g2 = gradient(p2)
    f1 = p1 - g1 * v1
    f2 = p2 - g2 * v2

    # DEBUG: print footpoints
    print_point(file_handle, f1)
    print_point(file_handle, f2)
    print_segments(file_handle, [p1, f1])
    print_segments(file_handle, [p2, f2])
    # t1 = [-g1[2], g1[1]]
    # t2 = [-g2[2], g2[1]]
    # d = norm(f1 - f2)
    # print_segments(file_handle, [f1 - t1 * d, f1 + t1 * d])
    # print_segments(file_handle, [f2 - t2 * d, f2 + t2 * d])

    lp = liming_parabola(f1, g1, f2, g2)
    intersect_implicit(lp, p1, p2)
end

function guess_normal(p, p1, p2)
    p === nothing && return nothing
    x = norm(p - p1) / norm(p2 - p1)
    gradient(p1) * (1 - x) + gradient(p2) * x
end

function sample_bezier(cp, resolution)
    result = []
    n = length(cp) - 1
    for u in range(0, stop=1, length=resolution)
        tmp = copy(cp)
        for k in 1:n, i in 1:n-k+1
            tmp[i] = tmp[i] * (1 - u) + tmp[i+1] * u
        end
        push!(result, tmp[1])
        if check_accuracy
            global accuracy
            d = abs(curve(tmp[1]))
            if d > accuracy
                accuracy = d
            end
        end
    end
    result
end

function print_curve(f, approx_type)
    approx_type === :real && return print_segments(f, real_curve)
    for i in 1:cells, j in 1:cells
        p = corner + [i - 1, j - 1] / cells * bbox_edge
        d = bbox_edge / cells
        points = [p, p + [d, 0], p + [0, d], p + [d, d]]
        if approx_type === :linear
            ints = [find_intersection(points[1], points[2]),
                    find_intersection(points[1], points[3]),
                    find_intersection(points[2], points[4]),
                    find_intersection(points[3], points[4])]
            endpoints = filter(x -> x != nothing, ints)
            !isempty(endpoints) && print_segments(f, sample_bezier(endpoints, sampling_res))
        elseif approx_type === :liming_cubic
            ints = [find_intersection2(points[1], points[2]),
                    find_intersection2(points[1], points[3]),
                    find_intersection2(points[2], points[4]),
                    find_intersection2(points[3], points[4])]
            normals = [guess_normal(ints[1], points[1], points[2]),
                       guess_normal(ints[2], points[1], points[3]),
                       guess_normal(ints[3], points[2], points[4]),
                       guess_normal(ints[4], points[3], points[4])]
            ints = filter(x -> x != nothing, ints)
            normals = filter(x -> x != nothing, normals)
            length(ints) != 2 && continue # TODO
            p1, n1 = ints[1], normals[1]
            p2, n2 = ints[2], normals[2]
            t1, t2 = normalize([-n1[2], n1[1]]), normalize([-n2[2], n2[1]])
            if dot(t1, p2 - p1) < 0
                t1 *= -1
            end
            if dot(t2, p1 - p2) < 0
                t2 *= -1
            end
            d = norm(p2 - p1) / 3
            # Control polygon
            # print_segments(f, [p1, p1 + t1 * d, p2 + t2 * d, p2])
            print_segments(f, sample_bezier([p1, p1 + t1 * d, p2 + t2 * d, p2], sampling_res))
        end
    end
end

function print_curves(f)
    for i in 1:length(show_types)
        c = colors[i,:]
        println(f, "$(c[1]) $(c[2]) $(c[3]) setrgbcolor")
        print_curve(f, show_types[i])
    end
    println(f, "0 0 0 setrgbcolor")
end

function print_settings(f)
    line_height = 24
    line_num = 1
    y = fullheight - width
    println(f, "/Times-Roman findfont 18 scalefont setfont")
    found = false
    for line in split(read("test2d.jl", String), "\n")
        if !found && match(r"^### SETUP_BEGIN.*", line) != nothing
            found = true
        elseif found && match(r"^### SETUP_END.*", line) != nothing
            break
        elseif found
            line_num += 1
            println(f, """$margin $(y - line_num * line_height) moveto
                          ($line) show""")
        end
    end
end

print_footer(f) = println(f, "showpage")


# User functions

function generate()
    open(filename, "w") do f
        global file_handle = f
        print_header(f)
        print_grid(f)
        print_curves(f)
        print_settings(f)
        print_footer(f)
    end
end

function approximate(tolerance; max_cells = 32)
    cells_old = cells
    global cells = 1
    global check_accuracy = true
    global accuracy = Inf
    while cells < max_cells && accuracy > tolerance
        cells += 1
        accuracy = 0
        try
            generate()
        catch
            accuracy = Inf      # Silently pass over failing cases
        end
    end
    println("Achieved accuracy: $accuracy, using $cells x $cells cells")
    check_accuracy = false
    cells = cells_old
    nothing
end

end # module
