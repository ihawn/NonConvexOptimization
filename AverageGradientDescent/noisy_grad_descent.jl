using Plots
using Calculus
using LinearAlgebra
include("C:/Users/Isaac/Documents/Optimization/NonConvex/NonConvexOptimization/NonConvexOptimiztion/testobjectives.jl")


function Compute_Gradient(_f, _x)
    return Calculus.gradient(g -> _f(g),_x)
end


function Compute_Hessian(_f, _x)
    return hessian(h -> _f(h),_x)
end


function Back_Line_Search(_x, _f, _fx, _∇f, _Δx, _α, _β, _κ)
    _t = _κ

    while _f((_x + _t*_Δx)) > _fx + _α*_t*transpose(_∇f)*_Δx
        _t *= _β
    end

    return _t
end


function Newton_Step(_∇f, _∇2f, _x)
    return _∇2f\-_∇f
end


function Newton_Decrement(_∇2f, _Δxnt)
    return transpose(_Δxnt)*_∇2f*_Δxnt
end


function Unconstrained_Newton(_f, _x, _α, _β, _κ, _ϵ, maxIt)

    λ2 = 1
    it = 0
    val = _f(_x)

    while λ2 > _ϵ && it <= maxIt
        val = _f(_x)
        ∇f = Compute_Gradient(_f, _x)
        ∇2f = Compute_Hessian(_f, _x)
        Δxnt = Newton_Step(∇f, ∇2f, _x)
        λ2 = Newton_Decrement(∇2f, Δxnt)
        t = Back_Line_Search(_x, _f, val, ∇f, Δxnt, _α, _β, _κ)
        _x += t*Δxnt
        it+=1
    end

    push!(finalSolX, _x[1])
    push!(finalSolY, _x[2])

    return _x, val
end


function Ave_Grad(_f, _x, _ℓ, _ρ)
    sum = Compute_Gradient(_f, _x)

    for i = 1:_ρ
        vec = (rand(length(_x)) .- 0.5) * _ℓ + _x
        sum += Compute_Gradient(_f, vec)

        push!(noiseX, vec[1])
        push!(noiseY, vec[2])
    end

    return sum/(_ρ + 1)
end


function Ave_Grad_Descent(_f, _x, _α, _β, _ϵ, _η, _κ, _ℓ, _ρ, maxIt)
    it = 0;
    t = 1
    ∇f = Compute_Gradient(_f, _x)
    normGrad = norm(∇f)

    val = _f(_x)

    push!(xPlot, _x[1])
    push!(yPlot, _x[2])

    while normGrad > _η && it < maxIt
        ∇f = Compute_Gradient(_f, _x)
        ∇f_ave = Ave_Grad(_f, _x, _ℓ, _ρ)
        normGrad = norm(∇f)

        bls = Back_Line_Search(_x, _f, val, ∇f, -∇f_ave, _α, _β, _κ)
        t = bls[1]
        _x -= t*∇f_ave
        val = _f(_x)

        push!(xPlot, _x[1])
        push!(yPlot, _x[2])

        it += 1
    end

    println("\nIterations: ", it)

    push!(solPlotX, _x[1])
    push!(solPlotY, _x[2])

    sol = Unconstrained_Newton(_f, _x, _α, _β, _κ, _ϵ, maxIt)

    return sol
end

flush(stdout)

x0 = [-4.0,2.0]
ϵ = 1e-16
η = 1e-2
η = 1e-3
α = 0.5
β = 0.8
κ = 1
ℓ = 5
ρ = 50



xPlot = []
yPlot = []
noiseX = []
noiseY = []
solPlotX = []
solPlotY = []
finalSolX = []
finalSolY = []

var = x0
maxIterations = 1e3

#f(x) = x[1]^2 + x[2]^2 + 7*sin(x[1] + x[2]) + 10*sin(5x[1])
#f(x) = (x[2] - 0.129*x[1]^2 + 1.6*x[1] - 6)^2 + 6.07*cos(x[1]) + 10
#f(x) = (3x[1] + 4x[2])^2 + x[1]^2*(1 - x[1])^2 + x[2]^2*(1 - x[2])^2
#f(x) = Rastrigin(x, 2)
#f(x) = Ackley(x)
#f(x) = Rosenbrock(x, 2)
#f(x) = Beale(x)
#f(x) = Bukin(x)
f(x) = Holder_Table(x)

@time minimum = Ave_Grad_Descent(f, x0, α, β, ϵ, η, κ, ℓ, ρ, maxIterations)
println(minimum)


plotf(x,y) = f([x, y])
_x = -10.0:0.03:10.0
_y = -10.0:0.03:10.0
X = repeat(reshape(_x, 1, :), length(_y), 1)
Y = repeat(_y, 1, length(_x))
Z = map(plotf, X, Y)
p1 = Plots.contour(_x,_y, plotf, fill = true)
plot(p1, legend = false, title = "Global Minimization With Noisy Gradient Descent")
plot!(xPlot, yPlot, color = "white")
scatter!(xPlot, yPlot, markersize = 2, color = "red")
#scatter!(noiseX, noiseY, markersize = 1)
scatter!(solPlotX, solPlotY, color = "green")
scatter!(finalSolX, finalSolY, color = "green")
