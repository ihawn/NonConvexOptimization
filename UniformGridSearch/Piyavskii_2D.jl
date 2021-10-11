using Plots
using LinearAlgebra


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

function IntersectPyramids(pyr_lst, c1, c2, c3, L)
    int_lst = Vector{Float64}[]
    hyper_lst = vcat(c1.planes, c2.planes, c3.planes)
    for i in 1:length(hyper_lst)
        for j in i+1:length(hyper_lst)
            for k in j+1:length(hyper_lst)
                p1, p2, p3 = hyper_lst[[i, j, k]]
                A = [p1.a p1.b p1.c; p2.a p2.b p2.c; p3.a p3.b p3.c]
              
                if abs(det(A)) > 1e-10
                    push!(int_lst, IntersectFromHyperplane(p1, p2, p3))
                end
                    
            end
        end
    end

    #check if on the generator cones
    i = 1
    s = length(int_lst)
    while i <= s
        C1 = c1.peak[3] - L*norm(int_lst[i][1:2] - c1.peak[1:2], Inf)
        C2 = c2.peak[3] - L*norm(int_lst[i][1:2] - c2.peak[1:2], Inf)
        C3 = c3.peak[3] - L*norm(int_lst[i][1:2] - c3.peak[1:2], Inf)

        zhat = int_lst[i][3]
        if zhat - C1 > 1e-12 || zhat - C2 > 1e-12 || zhat - C3 > 1e-12
            deleteat!(int_lst, i)
            s-=1
        else
            i+=1
        end
    end

    #check if under one of the pyramids in the x_list
    # i = 0 
    # s = length(int_lst)
    # del = Int64[]
    # for i in 1:s
    #     for c in pyr_lst
    #         C = c.peak[3] - L*norm(int_lst[i][1:2] - c.peak[1:2], Inf)

    #         zhat = int_lst[i][3]
    #         @show (C, zhat)
    #         if C - zhat > 1e-12
    #             push!(del, i)
    #             println("deleted")
    #             break
    #         end
    #     end
    #     println()
    # end
    # deleteat!(int_lst, del)
    all_planes = Hyperplane[]
    del = Int64[]
    for p in pyr_lst
        for h in p.planes
            for n = 1:length(int_lst)
                if HyperEval(int_lst[n], h)
                    push!(del, n)
                end
            end

            push!(all_planes, h)
        end
    end
    deleteat!(int_lst, unique(del))

    return int_lst
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
            [x[1]+offset, x[2]-offset, fx-L],
            [x[1]-offset, x[2]-offset, fx-L]]

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

function IntersectFromHyperplane(h1, h2, h3)
    return [h1.a h1.b h1.c; h2.a h2.b h2.c; h3.a h3.b h3.c]\[-h1.d; -h2.d; -h3.d]
end

function GetIntersections(p1, p2, p3, points, prev_min)
    for i = 1:4
        for j = 1:4
            for k = 1:4
                A = [p1[i].a p1[i].b p1[i].c; p2[j].a p2[j].b p2[j].c; p3[k].a p3[k].b p3[k].c]
                if abs(det(A)) > 1e-8 #find a better way of doing this
                    p = IntersectFromHyperplane(p1[i], p2[j], p3[k])
                    @show A
                    if p[3] >= prev_min
                        push!(points, p)
                    end
                end
            end
        end
    end

    return points
end

function MinZ(lst, xb, yb)
    min_val = lst[1][3]
    min_pos = 1
    for k in 1:length(lst)
        if lst[k][3] < min_val && lst[k][1] > xb[1] && lst[k][1] < xb[2] &&
            lst[k][2] > yb[1] && lst[k][2] < yb[2]
            min_val = lst[k][3]
            min_pos = k
        end
    end
    return min_val, min_pos
end

f(x) = x[1]^2 + x[2]^2
minX = -5
maxX = 5
xBounds = [-2, 1]
yBounds = [-2, 1]
L = 4

x = [[xBounds[1], yBounds[1]], [xBounds[1], yBounds[2]], [xBounds[2], yBounds[1]], [xBounds[2], yBounds[2]]]
pyr_list = []
int_list = []



#construct first 4 pyramids from the box constraints
for i in 1:length(x)
    push!(pyr_list, GeneratePyramid(f(x[i]), x[i], L))
end


#o = IntersectPyramids(pyr_list, pyr_list[1], pyr_list[2], pyr_list[3], L)
#compute first intersections
int_list = GetIntersections(pyr_list[1], pyr_list[2], pyr_list[3], int_list, -Inf)
int_list = GetIntersections(pyr_list[1], pyr_list[2], pyr_list[4], int_list, -Inf)
int_list = GetIntersections(pyr_list[1], pyr_list[3], pyr_list[4], int_list, -Inf)
int_list = GetIntersections(pyr_list[2], pyr_list[3], pyr_list[4], int_list, -Inf)

# pl = []
# for i in 1:10
#     min_int = MinZ(int_list, xBounds, yBounds)
#     min_x = int_list[min_int[2]][1:2]
#     deleteat!(int_list, min_int[2])

#     push!(pl,min_x) #for plotting

#     new_pyr = GeneratePyramid(f(min_x), min_x, L)
    

#     for n = 1:length(pyr_list)
#         for p = 1:length(pyr_list)
#             global int_list = GetIntersections(pyr_list[n], pyr_list[p], new_pyr, int_list, min_int[1])
#         end
#     end
#     @show i
#     @show f(min_x)

#     push!(pyr_list, new_pyr)

# end



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
