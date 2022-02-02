using Plots
using LinearAlgebra
using Combinatorics
using Profile
using Calculus
include("test_objectives.jl")


#ax + by + cz + d = 0
struct Hyperplane
    a::Float64
    b::Float64
    c::Float64
    d::Float64
end

mutable struct Pyramid
    planes::Vector{Hyperplane}
    peak::Vector{Float64}
    closed::Bool
end

struct Point
    x::Float64
    y::Float64
    z::Float64
end

function Compute_Gradient(_f, _x)
    return Calculus.gradient(g -> _f(g),_x)
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
    return Pyramid(h, [x[1], x[2], fx], false);
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

function AddIntersections!(int_lst, c1, c2, c3, L, bl)
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

                    should_add = true

                    for b in 1:length(bl)
                        if abs(sect[3] - bl[b]) <= 1e-8
                            should_add = false
                            break
                        end
                    end

                    if sum >= 2 && should_add
                        push!(int_lst, sect)
                    end
                end  
            end
        end
    end
end

#check if under one of the pyramids in the x_list
function FilterIntersections!(int_lst, pyr_lst, L)
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

function IntersectPyramids!(int_lst, pyr_lst, c1, c2, c3, L, bl)
    AddIntersections!(int_lst, c1, c2, c3, L, bl)
    #FilterIntersections!(int_lst, pyr_lst, L)
end

function ConeSearch()

    plot_x = Float64[]
    plot_y = Float64[]

    f(x) = Styblinski_Tang(x, 2)
    xBounds = [-3.8, 3.9]
    yBounds = [-3.6, 3.95]
    L = 86
    
    x = [[xBounds[1], yBounds[1]], [xBounds[1], yBounds[2]], [xBounds[2], yBounds[1]], [xBounds[2], yBounds[2]]]
    pyr_list = []
    int_list = []

    #construct first 4 pyramids from the box constraints
    for i in 1:length(x)
        push!(pyr_list, GeneratePyramid(f(x[i]), x[i], L))
    end

    sect_list = Vector{Float64}[]
    c = combinations(pyr_list, 3)

    
    blacklist = Float64[]

    for _c in c
        IntersectPyramids!(sect_list, pyr_list, _c[1], _c[2], _c[3], L, blacklist)
    end


    minf = [Inf, Inf, Inf]
    push!(plot_x, xBounds[1])
    push!(plot_y, yBounds[1])

    ##
    ##LOOP START
    ##
    for i in 1:3
        m_pt, m_pos = MinZ(sect_list)
        fx = f(m_pt[1:2])

        if fx < minf[3]
            minf = [m_pt[1], m_pt[2], fx]
            push!(plot_x, m_pt[1])
            push!(plot_y, m_pt[2])
        end

        r = length(sect_list)
        y = 1
        while y < r
            if norm(sect_list[y] - m_pt, 2) <= 1e-10
                deleteat!(sect_list, y)
                r-=1
            end
            y+=1
        end

        #deleteat!(sect_list, m_pos)
        #append!(blacklist, m_pt[3])


        #determine each cone triplet with closest proximity
        # all_groups = Vector{Pyramid}[]
        # for n in 1:length(pyr_list)
        #     pyr_group = [pyr_list[n], pyr_list[2], pyr_list[3]]
        #     min_dist = Inf

        #     for m in n+1:length(pyr_list)
        #         for v in m+1:length(pyr_list)
        #             d = norm(pyr_list[n].peak[1:2] - pyr_list[m].peak[1:2], 2) + 
        #                 norm(pyr_list[m].peak[1:2] - pyr_list[v].peak[1:2], 2) + 
        #                 norm(pyr_list[n].peak[1:2] - pyr_list[v].peak[1:2], 2)

        #             if d < min_dist
        #                 pyr_group = [pyr_list[n], pyr_list[m], pyr_list[v]]
        #             end
        #         end
        #     end
        #     push!(all_groups, pyr_group)
        # end

        #@show length(all_groups)


        # for g in 1:length(all_groups)
        #     for h in 1:3
        #         println(all_groups[g][h].peak)
        #     end
        #     println()
        # end

        c2 = combinations(pyr_list, 3)
        @show length(c2)

        for _c in c2
            if !(_c[1].closed && _c[2].closed && _c[3].closed) #Check to make sure this combo of cones hasn't been accounted for
                IntersectPyramids!(sect_list, pyr_list, _c[1], _c[2], _c[3], L, blacklist)
                # _c[1].closed = true
                # _c[2].closed = true
                # _c[3].closed = true
            end
        end

        push!(pyr_list, GeneratePyramid(fx, m_pt[1:2], L))

        println()
        @show m_pt
        @show fx
        @show minf
        @show length(pyr_list)
        @show length(sect_list)
       # @show length(all_groups)
        @show blacklist
        @show sect_list
    end

    pks_x = [pyr_list[1].peak[1]]
    pks_y = [pyr_list[1].peak[2]]
    for p in pyr_list
        append!(pks_x, p.peak[1])
        append!(pks_y, p.peak[2])
    end

    for p in pyr_list
        println(p.peak)
    end

    return plot_x, plot_y, sect_list, pks_x, pks_y
end

@time plot_x, plot_y, sect_list, peaks_x, peaks_y = ConeSearch()
# @profile ConeSearch()
minX = -4
maxX = 4
f(x) = Styblinski_Tang(x, 2)

# max_g = 0
# for t in 1:1e7
#     x = [rand()*(3.9+3.8)-3.8, rand()*(3.6+3.95)-3.6]
#     ∇f = Compute_Gradient(f, x)
#     nr = norm(∇f, 2)
#     if nr > max_g
#         global max_g = nr
#     end
# end
# @show max_g

plotf(x,y) = f([x, y])
_x = minX:0.005*abs(maxX - minX):maxX
_y = minX:0.005*abs(maxX - minX):maxX
X = repeat(reshape(_x, 1, :), length(_y), 1)
Y = repeat(_y, 1, length(_x))
Z = map(plotf, X, Y)
p1 = Plots.contour(_x,_y, plotf, fill = true, aspect_ratio=:equal)
plot(p1, xrange = (minX, maxX), yrange = (minX, maxX), legendfontsize = 4, dpi = 400, legend = false)
@show plot_x
@show plot_y
plot!(plot_x, plot_y, color=:green)
scatter!(plot_x, plot_y, color=:green)
scatter!(peaks_x, peaks_y, color=:red)