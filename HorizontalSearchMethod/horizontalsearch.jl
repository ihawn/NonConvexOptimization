using Plots
using Calculus
using LinearAlgebra

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

function Newton(_f, _x, _ϵ)
    ∇f = Compute_Gradient(_f, _x)
    ∇2f = Compute_Hessian(_f, _x)

    while abs(∇f) > _ϵ
        _x -= ∇f/∇2f
        ∇f = Compute_Gradient(_f, _x)
        ∇2f = Compute_Hessian(_f, _x)

    end

    return _x, _f(_x)
end

function Newton_Root(_f, _x, _ϵ, b)
    fx = _f(_x) - b
    ∇f = Compute_Gradient(_f, _x)

    while abs(fx) > _ϵ
        _x -= fx / ∇f
        ∇f = Compute_Gradient(_f, _x)
        fx = _f(_x) - b
    end

    return _x, _f(_x) - b
end

function Grad_Descent(_f, _x, _α, _β, _ϵ, _κ)
    it = 0;
    t = 1
    ∇f = Compute_Gradient(_f, _x)
    normGrad = norm(∇f)

    push!(xPlot, _x)
    push!(yPlot, _f(_x))

    val = _f(_x)

    while normGrad > _ϵ
        ∇f = Compute_Gradient(_f, _x)
        normGrad = norm(∇f)

        bls = Back_Line_Search(_x, _f, val, ∇f, -∇f, _α, _β, _κ)
        t = bls[1]
        _x -= t*∇f
        val = _f(_x)

        push!(xPlot, _x)
        push!(yPlot, _f(_x))


    end
    push!(xSol, _x)
    push!(ySol, _f(_x))

    return _x, val
end


function Adj_F(_f, _x, b)
    return _f(_x) - b
end


function Search(local_sol, _f, _ℓ, _γ, _η)
    while _f(local_sol[1] + _ℓ) >= local_sol[2] + _η && _f(local_sol[1] - _ℓ) >= local_sol[2] + _η && _ℓ > _η
        _ℓ *= _γ
        print("\nl = ", ℓ)
    end

    if _f(local_sol[1] + _ℓ) < local_sol[2] + _η
        return local_sol[1] + _ℓ
    else
        return local_sol[1] - _ℓ
    end
end


function Horizontal_Search(_f, _x, _α, _β, _η, _ϵ, _κ, _ℓ, _γ, width, maxIt)


    sol = Grad_Descent(_f, _x, _α, _β, _η, _κ)
    _x = Search(sol, _f, _ℓ, _γ, _η)

    sol = Grad_Descent(_f, _x, _α, _β, _η, _κ)


    return sol
end



f(x) = x^2 + sin(3x^2)

x0 =  5.0
ϵ = 1e-8
η = 1e-1
α = 0.5
β = 0.8
κ = 0.02
ℓ = 10
γ = 0.8
searchWidth = 10
xPlot = []
yPlot = []
xSol = []
ySol = []
xSolFinal = []
ySolFinal = []

var = x0
maxIterations = 150


minimum = Horizontal_Search(f, x0, α, β, η, ϵ, κ, ℓ, γ, searchWidth, maxIterations)
xlist = range(-3, 5, length = 1000)
plot( xlist, f.(xlist), legend = false)
scatter!(xSol, ySol)
#scatter!(xSolFinal, ySolFinal, color = "Red")
