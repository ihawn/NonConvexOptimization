using Zygote
using LinearAlgebra
using Plots
include("C:/Users/Isaac/Documents/Optimization/NonConvex/NonConvexOptimization/NonConvexOptimiztion/testobjectives.jl")


#Assumptions: f is differentiable and has Lipschitz first derivative

function Grad(_f, _x)
    return _f'(_x)
end


#inner product of gradient with d normalized by dividing by Δ
function uvp(arg, _d)
    _n = length(arg)
    ∇f = Grad(f, arg)
    sum = 0
    for i = 1:_n
        sum += ∇f[i]*_d[i]
    end
    return sum
end

n = 2
r = 2.0 #r > 1
f(x) = GP(x)
α = [-2.0, -2.0] #Lower left corner of hypercube
β = [2.0, 2.0] #Upper right corner of hypercube
D = (α, β) #Hypercube
d = β - α #length of hypercube edges
Δ = norm(d)
v = f(α)
u = f(β)
vp = dot(Grad(f, α), d)/Δ#uvp(α, d)
up = dot(Grad(f, β), d)/Δ#uvp(β, d)
M_list = []
D_list = [D]
R_list = []
δ_list = []
f_list = []
y_list = []
y1_list = []
y2_list = []
DPlot_list = []
Plot_list = []
Accepted_D_list = []
maxIt = 500
ϵ = Δ*0.01

minVal = Inf
for c = 1:maxIt


    #Compute the absolute values of the slopes over the hypercube
    m1 = abs(up - vp)/Δ
    m2 = 2*(-(u - v) + vp*Δ)/Δ^2
    m3 = 2*((u - v) - up*Δ)/Δ^2
    _M = max(m1, m2, m3)
    append!(M_list, _M)


    #Approximate the Lipschitz constant L
    M = maximum(M_list)
    m = 0.0
    if M == 0
        m = 1.0
    else
        m = r*M
    end


    #Calculate the characteristic for each hypercube
    _m = 0.5m*(1 + (_M/m)^2)
    δ = (-(u - v) + up*Δ + 0.5*_m*Δ^2)/(_m*Δ + (up - vp))
    R = v + vp*δ - 0.5*_m*δ^2
    append!(R_list, R)
    append!(δ_list, δ)


    #Find the hypercube D with minimal R
    minI = findmin(R_list)[2]
    minD = D_list[minI]


    #New trial points
    _α, _β = minD
    _d = _β - _α
    Δ = norm(_d)
    #dl, l = findmax(_d)


    # v = f(_α)
    # u = f(_β)
    # vp = dot(Grad(f, _α), _d)/Δ
    # up = dot(Grad(f, _β), _d)/Δ


    if _β[1] - _α[1] > _β[2] - _α[2]
        dl = _β[1] - _α[1]
        y1 = [_α[1] + δ*dl/Δ, _β[2]]
        y2 = [_α[1] + δ*dl/Δ, _α[2]]
    else
        dl = _β[2] - _α[2]
        y1 = [_β[1], _α[2] + δ*dl/Δ]
        y2 = [_α[1], _α[2] + δ*dl/Δ]
    end

    push!(y1_list, y1)
    push!(y2_list, y2)
    # minδ = δ_list[minI]
    # y1 = copy(_β)
    # y1[l] = _α[l] + minδ*dl/Δ
    # y2 = copy(_α)
    # y2[l] = y1[l]


    #Subdivide D into two hypercubes
    D1 = (_α, y1)
    D2 = (y2, _β)

    append!(f_list, f(y1))
    append!(f_list, f(y2))


    push!(Accepted_D_list, D_list[minI])
    deleteat!(D_list, minI)

    push!(D_list, D1)
    push!(D_list, D2)
    push!(y_list, y1)
    push!(y_list, y2)

    fx = f(y1)
    if fx < minVal
        minVal = fx
    end

    push!(Plot_list, y1)
    push!(DPlot_list, D1)
    push!(DPlot_list, D2)


    println()
    println(minVal)


    if Δ <= ϵ
        println("Condition met at i = ", c)
        break
    end
end

sol = findmin(f_list)
println("\nSolution: ", sol[1])
println("at x = ", y_list[sol[2]])


if n == 2
    zoomMult = 1
    plotf(x,y) = f([x, y])
    _x = zoomMult*α[1]:zoomMult*0.05:zoomMult*β[1]
    _y = zoomMult*α[1]:zoomMult*0.05:zoomMult*β[1]
    X = repeat(reshape(_x, 1, :), length(_y), 1)
    Y = repeat(_y, 1, length(_x))
    Z = map(plotf, X, Y)
    p1 = Plots.contour(_x,_y, plotf, fill = true, aspect_ratio=:equal, levels=[10,100,1000,10000,100000,1200000])
    plot(p1, xrange = (zoomMult*α[1], zoomMult*β[1]), yrange = (zoomMult*α[1], zoomMult*β[1]), legend = false, dpi = 400)
    #scatter!(Plot_list, color=:Green)
    #scatter!(y_list)
    scatter!([y[1] for y in y1_list], [y[2] for y in y1_list], color=:white)
    scatter!([y[1] for y in y2_list], [y[2] for y in y2_list], color=:blue)

    # for i in 1:length(Accepted_D_list)
    #     a, b = Accepted_D_list[i]
    #
    #     p = plot!([a[1], b[1]],[a[2],a[2]],color=:green, linewidth = 1)
    #     p = plot!([b[1], b[1]],[a[2],b[2]],color=:green, linewidth = 1)
    #     p = plot!([a[1], b[1]],[b[2],b[2]],color=:green, linewidth = 1)
    #     p = plot!([a[1], a[1]],[a[2],b[2]],color=:green, linewidth = 1)
    #
    #     display(p)
    # end

    scatter!([y[1] for y in y1_list], [y[2] for y in y1_list], color=:white)
    scatter!([y[1] for y in y2_list], [y[2] for y in y2_list], color=:blue)

end
