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


function Grad_Descent(_f, _x, _α, _β, _ϵ, _κ)
    it = 0;
    t = 1
    ∇f = Compute_Gradient(_f, _x)
    normGrad = norm(∇f)

    val = _f(_x)

    push!(xPlot, _x[1])
    push!(yPlot, _x[2])

    while normGrad > _ϵ
        ∇f = Compute_Gradient(_f, _x)
        normGrad = norm(∇f)

        bls = Back_Line_Search(_x, _f, val, ∇f, -∇f, _α, _β, _κ)
        t = bls[1]
        _x -= t*∇f
        val = _f(_x)

        push!(xPlot, _x[1])
        push!(yPlot, _x[2])
    end

    push!(solPlotX, _x[1])
    push!(solPlotY, _x[2])

    return _x, val
end


function Generate_ℓ_Vector(_x, _ℓ)
    return (rand(length(_x)) .- 0.5) * _ℓ
end


function Metropolis_Accept(val_new, val_old, _T)

    if val_new < val_old
        return true
    else
        return rand() >= exp(min(0, (val_old - val_new)/_T))
    end
end


function Monte_Carlo_Step(old_val, _f, _x, _ℓ, _η, _α, _β, _κ, _T)

    vec = Generate_ℓ_Vector(_x, _ℓ)
    local_sol = Grad_Descent(_f, _x + vec, _α, _β, _η, _κ)
    accept = Metropolis_Accept(local_sol[2], old_val, _T)

    push!(searchX, (local_sol[1] + vec)[1])
    push!(searchY, (local_sol[1] + vec)[2])

    return accept, local_sol
end


function Adjust_Step_Length(target_rate, rate, scale_factor, _ℓ)

    if rate > target_rate #stuck in a local min so increase step length
        _ℓ /= scale_factor
    else #random search is all over the place so take smaller steps
        _ℓ *= scale_factor
    end

    return _ℓ
end


function Basin_Hopping(_f, _x, _α, _β, _η, _ϵ, _κ, _ℓ, _γ, _ρ, _T, target_rate, width, maxIt)

    acceptance = 0
    rejection = 0

    for i = 1:maxIt
        println("Iteration ", i)
        sol = Grad_Descent(_f, _x, _α, _β, _η, _κ)
        accept = Monte_Carlo_Step(sol[2], _f, _x, _ℓ, _η, _α, _β, _κ, _T)

        if accept[1]
            _x = accept[2][1]
            acceptance += 1
        else
            rejection += 1
        end

        acceptance_rate = acceptance/(acceptance + rejection)
        _ℓ = Adjust_Step_Length(target_rate, acceptance_rate, γ, _ℓ)
    end

    sol = Unconstrained_Newton(_f, _x, _α, _β, _κ, _ϵ, maxIt)

    println(sol)

    return sol
end

flush(stdout)

x0 = [-10.0,10.0]
ϵ = 1e-8
η = 1e-3
α = 0.5
β = 0.8
κ = 0.1
ℓ = 1
γ = 0.5
ρ = 2
T = 10
searchWidth = 10
target_acc_rate = 0.5


xPlot = []
yPlot = []
solPlotX = []
solPlotY = []
searchX = []
searchY = []
finalSolX = []
finalSolY = []

var = x0
maxIterations = 5e2

#f(x) = x[1]^2 + x[2]^2 + 7*sin(x[1] + x[2]) + 10*sin(5x[1])
#f(x) = (x[2] - 0.129*x[1]^2 + 1.6*x[1] - 6)^2 + 6.07*cos(x[1]) + 10
#f(x) = Rastrigin(x, 2)
f(x) = Ackley(x)
#f(x) = Bukin(x)
#f(x) = Holder_Table(x)
#f(x) = Schaffer_N2(x)
#f(x) = Styblinski_Tang(x,2)

@time minimum = Basin_Hopping(f, x0, α, β, η, ϵ, κ, ℓ, γ, ρ, T, searchWidth, maxIterations)

plotf(x,y) = f([x, y])
_x = -10.0:0.03:10.0
_y = -10.0:0.03:10.0
#_x = 0.0:0.03:15.0
#_y = -5.0:0.03:20.0
X = repeat(reshape(_x, 1, :), length(_y), 1)
Y = repeat(_y, 1, length(_x))
Z = map(plotf, X, Y)
p1 = Plots.contour(_x,_y, plotf, fill = true)
plot(p1, legend = false, title = "Global Minimization With Basin Hopping")
scatter!(searchX, searchY, markersize = 2.5, color = "blue")
plot!(xPlot, yPlot, color = "white")
scatter!(xPlot, yPlot, markersize = 2, color = "red")
scatter!(solPlotX, solPlotY, color = "green")
scatter!(finalSolX, finalSolY, color = "white")
