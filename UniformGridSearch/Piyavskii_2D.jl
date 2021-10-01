using Plots
using LinearAlgebra


#ax + by + cz + d = 0
struct Hyperplane
    a::Float64
    b::Float64
    c::Float64
    d::Float64
end


function GenerateHyperplane(x)
    n = nullspace([x ones(length(x[:,1]))])
    return Hyperplane(n[1], n[2], n[3], n[4])
end

function GeneratePyramid(fx, x, L)
    offset = 1
    base = [[x[1]+offset, x[2]+offset, fx-L],
            [x[1]-offset, x[2]+offset, fx-L],
            [x[1]+offset, x[2]-offset, fx-L],
            [x[1]-offset, x[2]-offset, fx-L]]

    h = []
    for i in 1:4
        n = 1
        if i == 4
            offset = -3
        end
        x = [base[i] base[i+offset] [x[1], x[2], fx]]'
        push!(h, GenerateHyperplane(x))
    end
    return h;
end

function IntersectFromHyperplane(h1, h2, h3)
    return [h1.a h1.b h1.c; h2.a h2.b h2.c; h3.a h3.b h3.c]\[-h1.d; -h2.d; -h3.d]
end

function GetIntersections(p1, p2, p3)
    points = []
    for i = 1:4
        for j = 1:4
            for k = 1:4
                if i != j && i != k && j != k
                    append!(points, IntersectFromHyperplane(p1[i], p2[j], p3[k]))
                end
            end
        end
    end

    return points
end

f(x) = x[1]^2 + x[2]^2
minX = -5
maxX = 5
xBounds = [-2, 2]
yBounds = [-2, 2]

x = [[xBounds[1], yBounds[1]], [xBounds[1], yBounds[2]], [xBounds[2], yBounds[1]], [xBounds[2], yBounds[2]]]
p1 = GeneratePyramid(f(x[1]), x[1], 8)
p2 = GeneratePyramid(f(x[2]), x[2], 8)
p3 = GeneratePyramid(f(x[3]), x[3], 8)
ints = GetIntersections(p1, p2, p3)


plotf(x,y) = f([x, y])
_x = minX:0.005*abs(maxX - minX):maxX
_y = minX:0.005*abs(maxX - minX):maxX
X = repeat(reshape(_x, 1, :), length(_y), 1)
Y = repeat(_y, 1, length(_x))
Z = map(plotf, X, Y)
p1 = Plots.contour(_x,_y, plotf, fill = true, aspect_ratio=:equal)
plot(p1, xrange = (minX, maxX), yrange = (minX, maxX), legendfontsize = 4, dpi = 400, legend = false)
scatter!(xBounds, yBounds, color=:green)
scatter!(ints)
