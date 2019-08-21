module Marching2D

using LinearAlgebra

# Global settings

filename = "/tmp/test2d.eps"

### SETUP_BEGIN - do not modify this line
corner = [-1.4, -1.55]
bbox_edge = 3.0
cells = 3
curve = p -> sqrt(p[1]^2 + p[2]^2) - 1
gradient = p -> normalize(p)
real_curve = [[cos(x), sin(x)] for x in 0:0.01:2pi]
show_types = [:real :linear :parabolic]
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

const Point = Vector{Float64}

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

function find_intersection(p1, p2)
    v1 = curve(p1)
    v2 = curve(p2)
    v1 * v2 > 0 && return nothing
    x = abs(v1) / abs(v2 - v1)
    p1 * (1 - x) + p2 * x
end

function find_intersection2(p1, p2)
    v1 = curve(p1)
    v2 = curve(p2)
    v1 * v2 > 0 && return nothing
    g1 = gradient(p1)
    g2 = gradient(p2)
    f1 = p1 - g1 * v1
    f2 = p2 - g2 * v2
    # lp = liming_parabola(f1, g1, f2, g2)
    (p1 + p2) / 2               # TODO
end

function print_curve(f, type)
    type === :real && return print_segments(f, real_curve)
    for i in 1:cells, j in 1:cells
        p = corner + [i - 1, j - 1] / cells * bbox_edge
        d = bbox_edge / cells
        points = [p, p + [d, 0], p + [0, d], p + [d, d]]
        if type === :linear
            ints = [find_intersection(points[1], points[2]),
                    find_intersection(points[1], points[3]),
                    find_intersection(points[2], points[4]),
                    find_intersection(points[3], points[4])]
            print_segments(f, filter(x -> x != nothing, ints))
        elseif type == :parabolic
            ints = [find_intersection2(points[1], points[2]),
                    find_intersection2(points[1], points[3]),
                    find_intersection2(points[2], points[4]),
                    find_intersection2(points[3], points[4])]
            print_segments(f, filter(x -> x != nothing, ints))
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

function generate()
    open(filename, "w") do f
        print_header(f)
        print_grid(f)
        print_curves(f)
        print_settings(f)
        print_footer(f)
    end
end

end # module
