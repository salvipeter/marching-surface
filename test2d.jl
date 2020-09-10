module Marching2D

using ForwardDiff
using LinearAlgebra

import ..IpatchApprox

# Global settings

filename = "/tmp/test2d.eps"
sampling_res = 100
ground_truth_resolution = 100

### SETUP_BEGIN - do not modify this line
euclidean_distance = true
cells = 5
offsets = [0.9, 0.6, 0.3, -0.3, -0.6, -0.9]
onsets = [0.9, 0.6, 0.3, -0.3, -0.6, -0.9]

corner = [-9.5, -6]
bbox_edge = 12.0
curve_original = p -> p[1]^4 + 8p[1]^3 + p[2]^4 - 16

show_types = [:real :nolinear :liming_cubic :skip :nooffsets :noonsets]
#              red     green      blue      yellow    pink      cyan
### SETUP_END - do not modify this line


# Sample settings

# Blob
# corner = [-9.5, -6]
# bbox_edge = 12.0
# curve_original = p -> p[1]^4 + 8p[1]^3 + p[2]^4 - 16

# X
# corner = [-4.05, -4.05]
# bbox_edge = 8.0
# curve_original = p -> p[1]^2 - p[2]^2 - p[1]*p[2]*sin(p[1]*p[2])

# Two eggs
# corner = [-4.05, -4.05]
# bbox_edge = 8.0
# curve_original = p -> -4p[1] + 10p[1]^2 + p[2]^(-2) + p[2]^2 - 11

# Ellipse
# corner = [-1.4, -1.55]
# bbox_edge = 3.0
# ellipse = [1.2, 0.8]
# curve_original = p -> (p[1]/ellipse[1])^2 + (p[2]/ellipse[2])^2 - 1

# Egg + line
# corner = [-10.1, -10.1]
# bbox_edge = 20.0
# curve_original = p -> p[2]^2 - p[1]^3 + 7p[1] - 6


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

# Global variables

global accuracy
check_accuracy = false
global curve_old
global gradient_old
global file_handle              # for debug


# Wrapper functions

gradient_original = p -> ForwardDiff.gradient(curve_original, p)
curve = p -> curve_original(p) / norm(gradient_original(p))
gradient = p -> normalize(gradient_original(p))


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

function print_individual_segments(f, segments)
    for (p, q) in segments
        p, q = convert(p), convert(q)
        println(f, """newpath
                      $(p[1]) $(p[2]) moveto
                      $(q[1]) $(q[2]) lineto
                      stroke""")
    end
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
    p = p1 * (1 - x) + p2 * x
    n = safe_normalize(gradient(p1) * (1 - x) + gradient(p2) * x)
    (p, n)
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

evalQuadratic(S, p) =
    S[1] * p[1]^2 + S[2] * p[1] * p[2] + S[3] * p[2]^2 + S[4] * p[1] + S[5] * p[2] + S[6]

function liming_parabola(f1, g1, f2, g2)
    # Assuming equation S3^2 - lambda S1 S2 = 0
    # Using vectors as polynomials
    # - linear: [x, y, 1]
    # - quadratic: [x^2, xy, y^2, x, y, 1]

    if dot(g1, f2 - f1) < 0
        g1 *= -1
    end
    S1 = [g1[1], g1[2], -dot(f1, g1)]

    if dot(g2, f1 - f2) < 0
        g2 *= -1
    end
    S2 = [g2[1], g2[2], -dot(f2, g2)]

    f3 = [S1'; S2'][:,1:2] \ [-S1[3], -S2[3]] # intersection
    d = f2 - f1
    n = normalize([-d[2], d[1]])
    if dot(n, f3 - f1) < 0
        n *= -1
    end
    S3 = [n[1], n[2], -dot(f1, n)]

    S1S2 = [S1[1]*S2[1], S1[1]*S2[2]+S1[2]*S2[1], S1[2]*S2[2],
            S1[1]*S2[3]+S1[3]*S2[1], S1[2]*S2[3]+S1[3]*S2[2], S1[3]*S2[3]]
    S3S3 = [S3[1]^2, 2*S3[1]*S3[2], S3[2]^2, 2*S3[1]*S3[3], 2*S3[2]*S3[3], S3[3]^2]

    p = f1 / 4 + f3 / 2 + f2 / 4 # A point on the curve
    lambda = evalQuadratic(S3S3, p) / evalQuadratic(S1S2, p)

    # Original Liming equation: (1 - m) L1 L2 - m L0^2 = 0
    m = 1 / (lambda + 1)
    if m < 0 || m > 1
        @warn "Liming constant outside [0, 1]: $m"
        lambda = 0
    end

    S3S3 - lambda * S1S2
end

function intersect_implicit(S, p1, p2)
    if p1[1] == p2[1]
        x = p1[1]
        y1, y2 = solve_quadratic(S[3], S[2] * x + S[5], S[1] * x^2 + S[4] * x + S[6])
        lo, hi = min(p1[2], p2[2]), max(p1[2], p2[2])
        y = y1 < lo || y1 > hi ? y2 : y1
        if lo <= y1 <= hi && lo <= y2 <= hi
            # Two valid solutions; take the one closer to the linear approximation
            yl = find_intersection(p1, p2)[1][2]
            y = abs(y1 - yl) < abs(y2 - yl) ? y1 : y2
        elseif y < lo || y > hi
            @warn "$y is not in [$lo, $hi]"
        end
        ([x, y], normalize([2 * S[1] * x + S[2] * y + S[4], 2 * S[3] * y + S[2] * x + S[5]]))
    else
        @assert p1[2] == p2[2] "The points should share exactly one coordinate"
        y = p1[2]
        x1, x2 = solve_quadratic(S[1], S[2] * y + S[4], S[3] * y^2 + S[5] * y + S[6])
        lo, hi = min(p1[1], p2[1]), max(p1[1], p2[1])
        x = x1 < lo || x1 > hi ? x2 : x1
        if lo <= x1 <= hi && lo <= x2 <= hi
            # Two valid solutions; take the one closer to the linear approximation
            xl = find_intersection(p1, p2)[1][1]
            x = abs(x1 - xl) < abs(x2 - xl) ? x1 : x2
        elseif x < lo || x > hi
            @warn "$x is not in [$lo, $hi]"
        end
        ([x, y], normalize([2 * S[1] * x + S[2] * y + S[4], 2 * S[3] * y + S[2] * x + S[5]]))
    end
end

is_zero_vector(v) = norm(v) < 1.0e-8
safe_normalize(p) = is_zero_vector(p) ? p : normalize(p)

function find_intersection2(p1, p2)
    v1 = curve(p1)
    v2 = curve(p2)
    v1 * v2 > 0 && return nothing
    g1 = gradient(p1)
    g2 = gradient(p2)
    f1 = p1 - g1 * v1
    f2 = p2 - g2 * v2

    # DEBUG: print footpoint vectors
    # print_segments(file_handle, [p1, f1])
    # print_segments(file_handle, [p2, f2])

    # When the gradients are nearly parallel (< 25 degrees), fall back to linear
    abs(dot(g1, g2)) > 0.9 && return find_intersection(p1, p2)
    # Also when either is zero
    (is_zero_vector(g1) || is_zero_vector(g2)) && return find_intersection(p1, p2)

    lp = liming_parabola(f1, g1, f2, g2)
    intersect_implicit(lp, p1, p2)
end

function sample_bezier(cp, resolution)
    !check_accuracy && length(cp) == 2 && return cp
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

function print_offset_curve(f, offset)
    for i in 1:cells, j in 1:cells
        p = corner + [i - 1, j - 1] / cells * bbox_edge
        d = bbox_edge / cells
        points = [p, p + [d, 0], p + [0, d], p + [d, d]]
        ints = [find_intersection(points[1], points[2]),
                find_intersection(points[1], points[3]),
                find_intersection(points[2], points[4]),
                find_intersection(points[3], points[4])]
        endpoints = filter(x -> x != nothing, ints)
        if !isempty(endpoints)
            print_segments(f, [endpoints[1][1] + endpoints[1][2] * offset,
                               endpoints[2][1] + endpoints[2][2] * offset])
        end
    end
end

function with_high_resolution(f)
    cells_old = cells
    global cells = ground_truth_resolution
    global curve, curve_old = curve_old, curve
    global gradient, gradient_old = gradient_old, gradient
    f()
    cells = cells_old
    curve, curve_old = curve_old, curve
    gradient, gradient_old = gradient_old, gradient
end


# Local implicit evaluator

coeffDegree(coeff) = (isqrt(9 + 8 * (length(coeff) - 1)) - 3) ÷ 2

pointConstraint(p, degree) = [p[1]^i * p[2]^j for i in 0:degree for j in 0:degree-i]

evalSurface(coeffs, degree, p) = dot(pointConstraint(p, degree), coeffs)

function evalInSubCell(curve, degree, topleft, len, depth)
    p = [topleft, topleft + [len, 0], topleft + [len, len], topleft + [0, len]]
    v = [evalSurface(curve, degree, q) for q in p]
    (all(x -> sign(x) == 1, v) || all(x -> sign(x) == -1, v)) && return []

    result = []
    if depth > 0
        len = len / 2
        depth -= 1
        append!(result, evalInSubCell(curve, degree, topleft, len, depth))
        append!(result, evalInSubCell(curve, degree, topleft + [len, 0], len, depth))
        append!(result, evalInSubCell(curve, degree, topleft + [len, len], len, depth))
        append!(result, evalInSubCell(curve, degree, topleft + [0, len], len, depth))
        return result
    end

    intersections = []
    for i in 1:4
        j = mod1(i + 1, 4)
        if v[i] * v[j] < 0
            push!(intersections, i)
        end
    end
    length(intersections) != 2 && return []

    for i in intersections
        j = mod1(i + 1, 4)
        x = abs(v[i]) / abs(v[j] - v[i])
        push!(result, p[i] * (1 - x) + p[j] * x)
    end

    [result]
end

function evalInCell(curve, topleft, edge, initial_resolution, depth)
    result = []
    degree = coeffDegree(curve)
    len = edge / initial_resolution
    for x in range(topleft[1], step=len, length=initial_resolution)
        for y in range(topleft[2], step=len, length=initial_resolution)
            append!(result, evalInSubCell(curve, degree, [x, y], len, depth))
        end
    end
    result
end

function fit_ipatch(p1, n1, p2, n2, approx_points)
    IpatchApprox.coeffsToVector(IpatchApprox.approxIpatch(p1, n1, p2, n2, approx_points))
end


# Main function

function print_curve(f, approx_type)
    if approx_type === :real
        with_high_resolution() do
            print_curve(f, :linear)
        end
    end
    if approx_type === :offsets
        with_high_resolution() do
            for offset in offsets
                print_offset_curve(f, offset)
            end
        end
    end
    if approx_type === :onsets
        with_high_resolution() do
            curve_old = curve
            global cells = ground_truth_resolution
            for onset in onsets
                curve = p -> curve_old(p) - onset
                print_curve(f, :linear)
            end
            global curve = curve_old
        end
    end
    for i in 1:cells, j in 1:cells
        p = corner + [i - 1, j - 1] / cells * bbox_edge
        d = bbox_edge / cells
        points = [p, p + [d, 0], p + [0, d], p + [d, d]]
        if approx_type === :linear
            ints = [find_intersection(points[1], points[2]),
                    find_intersection(points[1], points[3]),
                    find_intersection(points[2], points[4]),
                    find_intersection(points[3], points[4])]
            endpoints = [x[1] for x in ints if x != nothing]
            !isempty(endpoints) && print_segments(f, sample_bezier(endpoints, sampling_res))
        elseif approx_type === :liming_cubic
            ints = [find_intersection2(points[1], points[2]),
                    find_intersection2(points[1], points[3]),
                    find_intersection2(points[2], points[4]),
                    find_intersection2(points[3], points[4])]
            dirs = [ints[1] === nothing ? nothing : [0, 1],
                    ints[2] === nothing ? nothing : [1, 0],
                    ints[3] === nothing ? nothing : [-1, 0],
                    ints[4] === nothing ? nothing : [0, -1]]
            ints = filter(x -> x != nothing, ints)
            dirs = filter(x -> x != nothing, dirs)
            length(ints) != 2 && continue # TODO
            p1, n1, d1 = ints[1][1], ints[1][2], dirs[1]
            p2, n2, d2 = ints[2][1], ints[2][2], dirs[2]
            t1, t2 = [-n1[2], n1[1]], [-n2[2], n2[1]]
            if dot(t1, d1) < 0
                t1 *= -1
            end
            if dot(t2, d2) < 0
                t2 *= -1
            end
            len = norm(p2 - p1) / 3

            # Compute footpoints
            foots = filter([q - gradient(q) * curve(q) for q in points]) do q
                all(q .>= p) && all(q .<= p + [d, d])
            end

            # DEBUG: print intersection points & tangents
            # print_point(f, p1)
            # print_point(f, p2)
            # print_segments(f, [p1, p1 + t1 * len])
            # print_segments(f, [p2, p2 + t2 * len])

            # Parametric curve
            # print_segments(f, [p1, p1 + t1 * len, p2 + t2 * len, p2]) # control polygon
            # print_segments(f, sample_bezier([p1, p1 + t1 * len, p2 + t2 * len, p2], sampling_res))

            # I-segment
            ipatch = fit_ipatch(p1, n1, p2, n2, foots)
            print_individual_segments(f, evalInCell(ipatch, p, d, 4, 4))
        end
    end
end

function print_curves(f)
    for i in 1:length(show_types)
        c = colors[:,i]
        println(f, "$(c[1]) $(c[2]) $(c[3]) setrgbcolor")
        print_curve(f, show_types[i])
    end
    println(f, "0 0 0 setrgbcolor")
end

function print_settings(f)
    line_height = 16
    line_num = 2
    y = fullheight - width
    println(f, "/Courier findfont 12 scalefont setfont")
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


# Euclidean distance helper functions

function sample_curve()
    result = []
    cells = ground_truth_resolution
    for i in 1:cells, j in 1:cells
        p = corner + [i - 1, j - 1] / cells * bbox_edge
        d = bbox_edge / cells
        points = [p, p + [d, 0], p + [0, d], p + [d, d]]
        ints = [find_intersection(points[1], points[2]),
                find_intersection(points[1], points[3]),
                find_intersection(points[2], points[4]),
                find_intersection(points[3], points[4])]
        endpoints = filter(x -> x != nothing, ints)
        if !isempty(endpoints)
            push!(result, endpoints[1][1], endpoints[2][1])
        end
    end
    result
end

function deviation(points, p)
    min = Inf
    local qmin
    for q in points
        d = norm(p - q)
        if d < min
            min = d
            qmin = q
        end
    end
    qmin - p
end

function with_euclidean_distance(f)
    if euclidean_distance
        global curve_old = curve
        global gradient_old = gradient
        points = sample_curve()
        global curve = p -> copysign(norm(deviation(points, p)), curve_old(p))
        global gradient = p -> -safe_normalize(deviation(points, p)) * sign(curve_old(p))
        f()
        curve = curve_old
        gradient = gradient_old
    else
        f()
    end
end


# User functions

function generate()
    with_euclidean_distance() do
        open(filename, "w") do f
            global file_handle = f
            print_header(f)
            print_grid(f)
            print_curves(f)
            print_settings(f)
            print_footer(f)
        end
    end
    nothing
end

function approximate(tolerance; max_cells = 32)
    with_euclidean_distance() do
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
            if accuracy != Inf
                println("Cells: $cells\tAccuracy: $accuracy")
            end
        end
        println("Achieved accuracy: $accuracy, using $cells x $cells cells")
        check_accuracy = false
        cells = cells_old
    end
    nothing
end

end # module
