module Cell2D

using LinearAlgebra

using Gtk
import Graphics


# Parameters

# GUI parameters
const width = 500
const height = 500
const point_size = 5
const tangent_size = 40
const corners = [[100, 100], [100, 400], [400, 400], [400, 100]]

# GUI variables
show_simple = true
show_implicit = true
show_tangents = true
show_intersections = true

# Global variables
global points = [[151.56365966796875, 196.9283447265625],
                 [53.91229248046875, 240.6383056640625],
                 [428.6517333984375, 239.539306640625],
                 [340.4219970703125, 208.86520385742188]]
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
    i = findmax(map(abs, n))[2] # index of max. absolute value in n
    j = mod1(i + 1, 2)
    n[i] * derivatives[j] - n[j] * derivatives[i]
end

function fitCurve(points, normals, degree)
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
        Graphics.set_source_rgb(ctx, 0, 0, 1)
        Graphics.set_line_width(ctx, 1.0)
        draw_polygon(ctx, implicit_curve)
    end

    # Footpoints (input)
    for p in [corners; points]
        Graphics.set_source_rgb(ctx, 0, 0, 0)
        Graphics.arc(ctx, p[1], p[2], point_size, 0, 2pi)
        Graphics.fill(ctx)
    end
    for i in 1:4
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
n        end
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

function generate_curve()
    global simple_curve = []
    global implicit_curve = []
    global intersections = []
end

mousedown_handler = @guarded (canvas, event) -> begin
    p = [event.x, event.y]
    global clicked = findfirst(points) do q
        norm(p - q) < 10
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

    showall(win)
end

run() = begin setup_gui(); nothing end

end
