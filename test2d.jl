module Marching2D

using LinearAlgebra

# Global settings

show_curve = true
filename = "/tmp/test2d.eps"

### SETUP_BEGIN - do not modify this line
corner = [-1.4, -1.55]
length = 3
cells = 3
curve = p -> sqrt(p[1]^2 + p[2]^2) - 1
gradient = p -> normalize(p)
real_curve = [[cos(x), sin(x)] for x in 0:0.01:2pi]
### SETUP_END - do not modify this line

# Global constants

fullwidth = 2^-0.25 / 4 * 100 / 2.54 * 72 # A4 = A0 / 4, from centimeters to inches to points
fullheight = fullwidth * sqrt(2)
margin = 0.25 * 72
width = fullwidth - 2 * margin
height = fullheight - 2 * margin

const Point = Vector{Float64}

function print_header(f)
    println(f, """%!PS-Adobe-2.0
                  %%BoundingBox: 0 0 $fullwidth $fullheight
                  %%EndComments
                  1 setlinecap
                  1 setlinejoin""")
end

function convert(p)
    q = (p - corner) / length * width .+ margin
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

function print_curve(f)
    !show_curve && return
    println(f, "1 0 0 setrgbcolor")
    print_segments(f, real_curve)
    println(f, "0 0 0 setrgbcolor")
end

function find_intersection(p1, p2)
    v1 = curve(p1)
    v2 = curve(p2)
    v1 * v2 > 0 && return nothing
    x = abs(v1) / abs(v2 - v1)
    p1 * (1 - x) + p2 * x
end

function cell_curve(i, j)
    p = corner + [i - 1, j - 1] / cells * length
    d = length / cells
    points = [p, p + [d, 0], p + [0, d], p + [d, d]]
    ints = [find_intersection(points[1], points[2]),
            find_intersection(points[1], points[3]),
            find_intersection(points[2], points[4]),
            find_intersection(points[3], points[4])]
    filter(x -> x != nothing, ints)
end

function print_marching(f)
    println(f, "0 0 1 setrgbcolor")
    for i in 1:cells, j in 1:cells
        print_segments(f, cell_curve(i, j))
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
        print_curve(f)
        print_marching(f)
        print_settings(f)
        print_footer(f)
    end
end

end # module
