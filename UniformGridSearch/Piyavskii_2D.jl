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

function GetIntersections(p1, p2, p3, points)
    for i = 1:4
        for j = 1:4
            for k = 1:4
                A = [p1[i].a p1[i].b p1[i].c; p2[j].a p2[j].b p2[j].c; p3[k].a p3[k].b p3[k].c]
                if abs(det(A)) > 1e-8 #find a better way of doing this
                    push!(points, IntersectFromHyperplane(p1[i], p2[j], p3[k]))
                end
            end
        end
    end

    return points
end

function MinZ(lst)
    min_val = lst[1][3]
    min_pos = 1
    for k in 1:length(lst)
        if lst[k][3] < min_val
            min_val = lst[k][3]
            min_pos = k
        end
    end
    return min_val, min_pos
end

f(x) = x[1]^2 + x[2]^2
minX = -5
maxX = 5
xBounds = [-2, 2.5]
yBounds = [-2, 4]
L = 8

x = [[xBounds[1], yBounds[1]], [xBounds[1], yBounds[2]], [xBounds[2], yBounds[1]], [xBounds[2], yBounds[2]]]
pyr_list = []
int_list = []

#construct first 4 pyramids from the box constraints
for i in 1:length(x)
    push!(pyr_list, GeneratePyramid(f(x[i]), x[i], L))
end

#compute first intersections
int_list = GetIntersections(pyr_list[1], pyr_list[2], pyr_list[3], int_list)
int_list = GetIntersections(pyr_list[1], pyr_list[2], pyr_list[4], int_list)
int_list = GetIntersections(pyr_list[1], pyr_list[3], pyr_list[4], int_list)
int_list = GetIntersections(pyr_list[2], pyr_list[3], pyr_list[4], int_list)


for i in 1:10
    min_int = MinZ(int_list)
    min_x = int_list[min_int[2]][1:2]

    new_pyr = GeneratePyramid(f(min_x), min_x, L)

    for n = 1:length(pyr_list)
        for p = 1:length(pyr_list)
            global int_list = GetIntersections(pyr_list[n], pyr_list[p], new_pyr, int_list)
        end
    end

    push!(pyr_list, new_pyr)

    @show i
    @show f(min_x)
end



plotf(x,y) = f([x, y])
_x = minX:0.005*abs(maxX - minX):maxX
_y = minX:0.005*abs(maxX - minX):maxX
X = repeat(reshape(_x, 1, :), length(_y), 1)
Y = repeat(_y, 1, length(_x))
Z = map(plotf, X, Y)
p1 = Plots.contour(_x,_y, plotf, fill = true, aspect_ratio=:equal)
plot(p1, xrange = (minX, maxX), yrange = (minX, maxX), legendfontsize = 4, dpi = 400, legend = false)
#scatter!(xBounds, yBounds, color=:green)

x_plot = []
y_plot = []
for i in 1:length(int_list[1])
    append!(x_plot, int_list[i][1])
    append!(y_plot, int_list[i][2])
end
scatter!(x_plot, y_plot)
