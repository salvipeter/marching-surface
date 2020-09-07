module Cell2D

using LinearAlgebra

using Gtk
import Graphics


# Parameters

# GUI parameters
const width = 500
const height = 500
const click_precision = 10
const point_size = 5
const tangent_size = 40
const marching_initial_resolution = 4
const marching_depth = 4
const corners = [[100, 100], [100, 400], [400, 400], [400, 100]]

# GUI variables
show_simple = true
show_implicit = true
show_tangents = true
show_intersections = true
simple_intersections = true

# Global variables
global points = [[151.56365966796875, 196.9283447265625],
                 [53.91229248046875, 240.6383056640625],
                 [428.6517333984375, 239.539306640625],
                 [340.4219970703125, 208.86520385742188]]
global signs = [1, -1, -1, 1]
global simple_curve
global implicit_curve
global intersections


# Implicit fitting

pointConstraint(p, degree) = [p[1]^i * p[2]^j for i in 0:degree for j in 0:degree-i]

evalSurface(coeffs, degree, p) = dot(pointConstraint(p, degree), coeffs)

function gradientConstraint(p, degree)
    pow(x,i) = i == 0 ? 0 : i * x^(i-1)
    ([pow(p[1],i) * p[2]^j for i in 0:degree for j in 0:degree-i],
     [p[1]^i * pow(p[2],j) for i in 0:degree for j in 0:degree-i])
end

function normalConstraint(p, n, degree)
    derivatives = gradientConstraint(p, degree)
    n[1] * derivatives[2] - n[2] * derivatives[1]
end

function fitImplicit(interpolation, approximation, degree)
    rows = []
    for (p, n) in interpolation
        push!(rows, pointConstraint(p, degree))
        push!(rows, normalConstraint(p, n, degree))
    end
    for (p, n) in approximation
        push!(rows, pointConstraint(p, degree))
    end
    A = mapreduce(transpose, vcat, rows)
    F = svd(A)
    x = F.V[:,end]
    x
end

function evalInCell(curve, topleft, axis, depth)
    p = [topleft, topleft + [axis[1], 0], topleft + axis, topleft + [0, axis[2]]]
    v = [evalSurface(curve, 3, q) for q in p]
    (all(x -> sign(x) == 1, v) || all(x -> sign(x) == -1, v)) && return []

    result = []
    if depth > 0
        axis = axis / 2
        depth -= 1
        append!(result, evalInCell(curve, topleft, axis, depth))
        append!(result, evalInCell(curve, topleft + [axis[1], 0], axis, depth))
        append!(result, evalInCell(curve, topleft + axis, axis, depth))
        append!(result, evalInCell(curve, topleft + [0, axis[2]], axis, depth))
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

function evalInCell(curve)
    result = []
    len = (corners[3] - corners[1]) / marching_initial_resolution
    for x in range(corners[1][1], step=len[1], length=marching_initial_resolution)
        for y in range(corners[1][2], step=len[2], length=marching_initial_resolution)
            append!(result, evalInCell(curve, [x, y], len, marching_depth))
        end
    end
    result
end


# Bezier evaluation

"""
    bernstein(n, u)

Computes the Bernstein polynomials of degree `n` at the parameter `u`.
"""
function bernstein(n, u)
    coeff = [1.0]
    for j in 1:n
        saved = 0.0
        for k in 1:j
            tmp = coeff[k]
            coeff[k] = saved + tmp * (1.0 - u)
            saved = tmp * u
        end
        push!(coeff, saved)
    end
    coeff
end

"""
    bernstein_all(n, u)

Computes all Bernstein polynomials up to degree `n` at the parameter `u`.
"""
function bernstein_all(n, u)
    result = [[1.0]]
    coeff = [1.0]
    for j in 1:n
        saved = 0.0
        for k in 1:j
            tmp = coeff[k]
            coeff[k] = saved + tmp * (1.0 - u)
            saved = tmp * u
        end
        push!(coeff, saved)
        push!(result, copy(coeff))
    end
    result
end

"""
    bezier_derivative_controls(curve, d)

Computes the control points for the `d`-th derivative computation.
The Bezier curve is given by its control points.
"""
function bezier_derivative_controls(curve, d)
    n = length(curve) - 1
    dcp = [copy(curve)]
    for k in 1:d
        tmp = n - k + 1
        cp = []
        for i in 1:tmp
            push!(cp, (dcp[k][i+1] - dcp[k][i]) * tmp)
        end
        push!(dcp, cp)
    end
    dcp
end

"""
    bezier_eval(curve, u)

Evaluates a Bezier curve, given by its control points, at the parameter `u`.
"""
function bezier_eval(curve, u)
    n = length(curve) - 1
    coeff = bernstein(n, u)
    sum(curve .* coeff)
end

"""
    bezier_eval(curve, u, d)

Evaluates a Bezier curve, given by its control points, at the parameter `u`, with `d` derivatives.
"""
function bezier_eval(curve, u, d)
    result = []
    n = length(curve) - 1
    du = min(d, n)
    coeff = bernstein_all(n, u)
    dcp = bezier_derivative_controls(curve, du)
    for k in 0:du
        push!(result, sum(dcp[k+1] .* coeff[n-k+1]))
    end
    for k in n+1:d
        push!(result, [0, 0])
    end
    result
end


# Main logic

is_zero_vector(v) = norm(v) < 1.0e-8
safe_normalize(p) = is_zero_vector(p) ? p : normalize(p)

gradient(i) = normalize(points[i] - corners[i])
distance(i) = norm(points[i] - corners[i]) * signs[i]

function find_intersection(i, j)
    v1 = distance(i)
    v2 = distance(j)
    v1 * v2 > 0 && return nothing
    x = abs(v1) / abs(v2 - v1)
    p = corners[i] * (1 - x) + corners[j] * x
    n = safe_normalize(gradient(i) * (x - 1) + gradient(j) * x)
    (p, n)
end

"""
    bisect(f, xl, xh, iterations = 100, ϵ = 1.0e-7)

`xl` and `xh` bracket the root of f.
"""
function bisect(f, xl, xh, iterations = 100, ϵ = 1.0e-7)
    xm = xl
    for i = 1:iterations
        xm_old = xm
        xm = (xl + xh) / 2
        if xm != 0 && abs((xm - xm_old) / xm) < ϵ
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

function intersect_curve_line(cpts, a, b)
    # The line is either horizontal or vertical
    c = a[1] == b[1] ? 1 : 2
    x = a[c]
    (cpts[1][c] - x) * (cpts[end][c] - x) > 0 && return nothing
    bisect(0, 1) do u
        bezier_eval(cpts, u)[c] - x
    end
end

function fit_cubic_bezier(p1, n1, p2, n2)
    len = norm(p1 - p2) / 3
    d1 = [n1[2], -n1[1]]
    d2 = [n2[2], -n2[1]]
    if dot(p2 - p1, d1) < 0
        d1 *= -1
    end
    if dot(p1 - p2, d2) < 0
        d2 *= -1
    end
    [p1, p1 + d1 * len, p2 + d2 * len, p2]
end

inside_segment(p, a, b) = dot(b - a, p - a) > 0 && norm(p - a) <= norm(b - a)

function find_intersection_with_curve(i, j)
    cpts = fit_cubic_bezier(points[i], gradient(i), points[j], gradient(j))

    u = intersect_curve_line(cpts, corners[i], corners[j])
    u === nothing && return nothing

    (p, d) = bezier_eval(cpts, u, 1)
    !inside_segment(p, corners[i], corners[j]) && return nothing
    n = normalize([d[2], -d[1]])
    (p, n)
end

function is_inside_cell(i)
    dir = [(1, 1), (1, -1), (-1, -1), (-1, 1)][i]
    all(x -> x >= 0, (points[i] - corners[i]) .* dir)
end

function generate_curve()
    global simple_curve = []
    global implicit_curve = []
    global intersections = []

    # Simple
    ints = [find_intersection(i, mod1(i + 1, 4)) for i in 1:4]
    ints = filter(x -> x != nothing, ints)
    if length(ints) == 2
        simple_curve = [ints[1][1], ints[2][1]]
        if simple_intersections
            intersections = ints
        end
    end

    # Implicit

    # Compute intersections
    ints = [find_intersection_with_curve(i, mod1(i + 1, 4)) for i in 1:4]
    ints = filter(x -> x != nothing, ints)
    if length(ints) != 2
        # Try the linear approximation
        ints = [find_intersection(i, mod1(i + 1, 4)) for i in 1:4]
        ints = filter(x -> x != nothing, ints)
        length(ints) != 2 && return
    end

    constraints = [(points[i], gradient(i)) for i in 1:4 if is_inside_cell(i)]

    # Parametric version
    # curve = fit_cubic_bezier(ints[1]..., ints[2]...)
    # samples = [bezier_eval(curve, u) for u in range(0, stop=1, length=51)]
    # implicit_curve = [[samples[i], samples[i+1]] for i in 1:50]

    # Implicit version
    curve = fitImplicit(ints, constraints, 3)
    implicit_curve = evalInCell(curve)

    if !simple_intersections
        intersections = ints
    end
end


# Graphics

function draw_polygon(ctx, poly, closep = false)
    if isempty(poly)
        return
    end
    Graphics.new_path(ctx)
    Graphics.move_to(ctx, poly[1][1], poly[1][2])
    for p in poly[2:end]
        Graphics.line_to(ctx, p[1], p[2])
    end
    if closep && length(poly) > 2
        Graphics.line_to(ctx, poly[1][1], poly[1][2])
    end
    Graphics.stroke(ctx)
end

function draw_segments(ctx, poly)
    n = length(poly)
    for i in 1:n
        Graphics.new_path(ctx)
        Graphics.move_to(ctx, poly[i][1][1], poly[i][1][2])
        Graphics.line_to(ctx, poly[i][2][1], poly[i][2][2])
        Graphics.stroke(ctx)
    end
end

draw_callback = @guarded (canvas) -> begin
    ctx = Graphics.getgc(canvas)

    # White background
    Graphics.rectangle(ctx, 0, 0, Graphics.width(canvas), Graphics.height(canvas))
    Graphics.set_source_rgb(ctx, 1, 1, 1)
    Graphics.fill(ctx)

    # Cell
    Graphics.set_source_rgb(ctx, 0, 0, 0)
    Graphics.set_line_width(ctx, 1.0)
    draw_polygon(ctx, corners, true)

    # Generated curves
    if show_simple
        Graphics.set_source_rgb(ctx, 0, 0, 1)
        Graphics.set_line_width(ctx, 1.0)
        draw_polygon(ctx, simple_curve)
    end
    if show_implicit
        Graphics.set_source_rgb(ctx, 0.5, 0.5, 1)
        Graphics.set_line_width(ctx, 5.0)
        draw_segments(ctx, implicit_curve)
    end

    # Input
    for p in [corners; points]
    end
    for i in 1:4
        p = corners[i]
        if signs[i] == 1
            Graphics.set_source_rgb(ctx, 0, 0, 0)
        else
            Graphics.set_source_rgb(ctx, 1, 0, 0)
        end
        Graphics.arc(ctx, p[1], p[2], point_size, 0, 2pi)
        Graphics.fill(ctx)
        # Footpoints
        p = points[i]
        Graphics.set_source_rgb(ctx, 0, 0, 0)
        Graphics.arc(ctx, p[1], p[2], point_size, 0, 2pi)
        Graphics.fill(ctx)
        # Footpoint lines
        Graphics.set_source_rgb(ctx, 0, 0, 0)
        Graphics.set_line_width(ctx, 1.0)
        line = [corners[i], points[i]]
        draw_polygon(ctx, line)
    end

    # Tangents
    if show_tangents
        Graphics.set_source_rgb(ctx, 0, 1, 0)
        Graphics.set_line_width(ctx, 2.0)
        for i in 1:4
            p = points[i]
            n = normalize(points[i] - corners[i])
            d = [n[2], -n[1]] * tangent_size
            line = [p + d, p - d]
            draw_polygon(ctx, line)
        end
    end

    if show_intersections
        for (p, n) in intersections
            Graphics.set_source_rgb(ctx, 0, 1, 1)
            Graphics.arc(ctx, p[1], p[2], point_size, 0, 2pi)
            Graphics.fill(ctx)
            if show_tangents
                Graphics.set_source_rgb(ctx, 1, 0, 1)
                Graphics.set_line_width(ctx, 2.0)
                d = [n[2], -n[1]] * tangent_size
                line = [p + d, p - d]
                draw_polygon(ctx, line)
            end
        end
    end
end


# GUI

mousedown_handler = @guarded (canvas, event) -> begin
    p = [event.x, event.y]
    global clicked = findfirst(points) do q
        norm(p - q) < click_precision
    end
    if clicked === nothing
        signchange = findfirst(corners) do q
            norm(p - q) < click_precision
        end
        if signchange != nothing
            signs[signchange] *= -1
            generate_curve()
            draw(canvas)
        end
    end
end

mousemove_handler = @guarded (canvas, event) -> begin
    clicked === nothing && return
    points[clicked] = [event.x, event.y]
    generate_curve()
    draw(canvas)
end

function clear_variables!()
    global simple_curve = []
    global implicit_curve = []
    global intersections = []
end

function setup_gui()
    clear_variables!()

    win = GtkWindow("2D Cell-Based Curve Reconstruction")
    vbox = GtkBox(:v)

    # Canvas widget
    canvas = GtkCanvas(width, height)
    canvas.mouse.button1press = mousedown_handler
    canvas.mouse.button1motion = mousemove_handler
    draw(draw_callback, canvas)
    push!(win, vbox)
    push!(vbox, canvas)

    hbox = GtkBox(:h)
    push!(vbox, hbox)

    simplep = GtkCheckButton("Show simple approximation")
    simplep.active[Bool] = show_simple
    signal_connect(simplep, "toggled") do cb
        global show_simple = cb.active[Bool]
        draw(canvas)
    end
    push!(hbox, simplep)

    implicitp = GtkCheckButton("Show implicit approximation")
    implicitp.active[Bool] = show_implicit
    signal_connect(implicitp, "toggled") do cb
        global show_implicit = cb.active[Bool]
        draw(canvas)
    end
    push!(hbox, implicitp)

    hbox = GtkBox(:h)
    push!(vbox, hbox)

    tangentp = GtkCheckButton("Show tangents")
    tangentp.active[Bool] = show_tangents
    signal_connect(tangentp, "toggled") do cb
        global show_tangents = cb.active[Bool]
        draw(canvas)
    end
    push!(hbox, tangentp)

    intersectionp = GtkCheckButton("Show intersections")
    intersectionp.active[Bool] = show_intersections
    signal_connect(intersectionp, "toggled") do cb
        global show_intersections = cb.active[Bool]
        draw(canvas)
    end
    push!(hbox, intersectionp)

    simple_intersectionp = GtkCheckButton("Show simple intersections")
    simple_intersectionp.active[Bool] = simple_intersections
    signal_connect(simple_intersectionp, "toggled") do cb
        global simple_intersections = cb.active[Bool]
        generate_curve()
        draw(canvas)
    end
    push!(hbox, simple_intersectionp)

    generate_curve()
    showall(win)
end

run() = begin setup_gui(); nothing end

end
