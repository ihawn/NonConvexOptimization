using Plots
using LinearAlgebra
using Combinatorics
using Profile


#ax + by + cz + d = 0
struct Hyperplane
    a::Float64
    b::Float64
    c::Float64
    d::Float64
end

struct Pyramid
    planes::Vector{Hyperplane}
    peak::Vector{Float64}
end

struct Point
    x::Float64
    y::Float64
    z::Float64
end


function HyperEval(p, h::Hyperplane)
    x = p[1:2]
    z = (h.a*x[1] + h.b*x[2] + h.d)/-h.c

    return p[3] - z < 1e-12
end

function GenerateHyperplane(x)
    n = nullspace([x ones(length(x[:,1]))])
    return Hyperplane(n[1], n[2], n[3], n[4])
end

function GeneratePyramid(fx, x, L)
    offset = 1
    base = [[x[1]+offset, x[2]+offset, fx-L],
            [x[1]-offset, x[2]+offset, fx-L],
            [x[1]-offset, x[2]-offset, fx-L],
            [x[1]+offset, x[2]-offset, fx-L]]

    h = []
    for i in 1:4
        n = 1
        if i == 4
            offset = -3
        end
        n = [base[i] base[i+offset] [x[1], x[2], fx]]'
        push!(h, GenerateHyperplane(n))
    end
    return Pyramid(h, [x[1], x[2], fx]);
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
    return lst[min_pos], min_pos
end


function IntersectFromHyperplane(h1, h2, h3)
    return [h1.a h1.b h1.c; h2.a h2.b h2.c; h3.a h3.b h3.c]\[-h1.d; -h2.d; -h3.d]
end

function GetIntersections(p1, p2, p3, points, prev_min)
    for i = 1:4
        for j = 1:4
            for k = 1:4
                A = [p1.planes[i].a p1.planes[i].b p1.planes[i].c; p2.planes[j].a p2.planes[j].b p2.planes[j].c; p3.planes[k].a p3.planes[k].b p3.planes[k].c]
                if abs(det(A)) > 1e-8 #find a better way of doing this
                    p = IntersectFromHyperplane(p1.planes[i], p2.planes[j], p3.planes[k])
 
                    if p[3] >= prev_min
                        push!(points, p)
                    end
                end
            end
        end
    end

    return points
end

function AddIntersections!(int_lst, c1, c2, c3, L)
    hyper_lst = vcat(c1.planes, c2.planes, c3.planes)

    for i in 1:length(hyper_lst)
        for j in i+1:length(hyper_lst)
            for k in j+1:length(hyper_lst)

                p1, p2, p3 = hyper_lst[[i, j, k]]
                A = [p1.a p1.b p1.c; p2.a p2.b p2.c; p3.a p3.b p3.c]

                if abs(det(A)) > 1e-10
                    sect = IntersectFromHyperplane(p1, p2, p3)

                    #check if on generator cones
                    c_list = Float64[]
                    push!(c_list, c1.peak[3] - L*norm(sect[1:2] - c1.peak[1:2], Inf))
                    push!(c_list, c2.peak[3] - L*norm(sect[1:2] - c2.peak[1:2], Inf))
                    push!(c_list, c3.peak[3] - L*norm(sect[1:2] - c3.peak[1:2], Inf))
                    zhat = sect[3]

                    sum = 0
                    for ii in 1:3
                        if zhat - c_list[ii] < 1e-12 
                            sum+=1
                        end
                    end

                    if sum >= 3
                        push!(int_lst, IntersectFromHyperplane(p1, p2, p3))
                    end
                end  
            end
        end
    end
end

#check if under one of the pyramids in the x_list
function FilterIntersections!(int_lst, pyr_lst)
    i = 0 
    s = length(int_lst)
    del = Int64[]
    for i in 1:s
        for c in pyr_lst
            C = c.peak[3] - L*norm(int_lst[i][1:2] - c.peak[1:2], Inf)

            zhat = int_lst[i][3]
            if C - zhat > 1e-12
                push!(del, i)
                break
            end
        end
    end
    deleteat!(int_lst, del)
end


function IntersectPyramids!(int_lst, pyr_lst, c1, c2, c3, L)

    AddIntersections!(int_lst, c1, c2, c3, L)
    FilterIntersections!(int_lst, pyr_lst)

    return int_lst
end

function ConeSearch()

    f(x) = x[1]^2 + x[2]^2
    minX = -5
    maxX = 5
    xBounds = [-1.8, 1.2]
    yBounds = [-2, 1.1]
    L = 4
    
    x = [[xBounds[1], yBounds[1]], [xBounds[1], yBounds[2]], [xBounds[2], yBounds[1]], [xBounds[2], yBounds[2]]]
    pyr_list = []
    int_list = []

    #construct first 4 pyramids from the box constraints
    for i in 1:length(x)
        push!(pyr_list, GeneratePyramid(f(x[i]), x[i], L))
    end

    sect_list = Vector{Float64}[]
    c = combinations(pyr_list, 3)
    for _c in c
        IntersectPyramids!(sect_list, pyr_list, _c[1], _c[2], _c[3], L)
    end

    minf = f([xBounds[1], yBounds[1]])
    for i in 1:30

        m_pt, m_pos = MinZ(sect_list)
        fx = f(m_pt[1:2])
        if fx < minf
            minf = fx
        end
        push!(pyr_list, GeneratePyramid(fx, m_pt[1:2], L))
        deleteat!(sect_list, m_pos)

        c2 = combinations(pyr_list, 3)
        @show length(c2)
        for _c in c2
            IntersectPyramids!(sect_list, pyr_list, _c[1], _c[2], _c[3], L)
        end

        println()
        @show m_pt[3]
        @show minf
        @show minf - m_pt[3]
        @show length(pyr_list)
        @show length(sect_list)
    end
end

ConeSearch()
# @profile ConeSearch()

# plotf(x,y) = f([x, y])
# _x = minX:0.005*abs(maxX - minX):maxX
# _y = minX:0.005*abs(maxX - minX):maxX
# X = repeat(reshape(_x, 1, :), length(_y), 1)
# Y = repeat(_y, 1, length(_x))
# Z = map(plotf, X, Y)
# p1 = Plots.contour(_x,_y, plotf, fill = true, aspect_ratio=:equal)
# plot(p1, xrange = (minX, maxX), yrange = (minX, maxX), legendfontsize = 4, dpi = 400, legend = false)
# #scatter!(xBounds, yBounds, color=:green)

# x_plot = []
# y_plot = []
# for i in 1:length(pl)
#     append!(x_plot, pl[i][1])
#     append!(y_plot, pl[i][2])
# end
# scatter!(x_plot, y_plot)
