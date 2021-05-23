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

    push!(xSolFinal, _x)
    push!(ySolFinal, _f(_x))

    return _x, _f(_x)
end

function Newton_Root(_f, _x, _ϵ)
    fx = _f(_x)
    ∇f = Compute_Gradient(_f, _x)

    while abs(fx) > _ϵ
        _x -= fx / ∇f
        ∇f = Compute_Gradient(_f, _x)
        fx = _f(_x)
    end



    return _x, _f(_x)
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

    push!(xSearchRight, local_sol[1])
    push!(ySearchRight, local_sol[2])
    push!(xSearchLeft, local_sol[1])
    push!(ySearchLeft, local_sol[2])

    push!(xSearchRight, local_sol[1] + _ℓ)
    push!(ySearchRight, local_sol[2])
    push!(xSearchLeft, local_sol[1] - _ℓ)
    push!(ySearchLeft, local_sol[2])

    while _f(local_sol[1] + _ℓ) >= local_sol[2] + _η && _f(local_sol[1] - _ℓ) >= local_sol[2] + _η && _ℓ > _η
        _ℓ *= _γ
        push!(xSearchRight, local_sol[1] + _ℓ)
        push!(ySearchRight, local_sol[2])
        push!(xSearchLeft, local_sol[1] - _ℓ)
        push!(ySearchLeft, local_sol[2])
    end

    if _f(local_sol[1] + _ℓ) < local_sol[2] + _η
        push!(xSearchRight, local_sol[1] + _ℓ)
        push!(ySearchRight, _f(local_sol[1] + _ℓ))
        return local_sol[1] + _ℓ, _ℓ
    else
        push!(xSearchLeft, local_sol[1] - _ℓ)
        push!(ySearchLeft, _f(local_sol[1] - _ℓ))
        return local_sol[1] - _ℓ, _ℓ
    end
end


function Horizontal_Search(_f, _x, _α, _β, _η, _ϵ, _κ, _ℓ, _γ, width, maxIt)

    sol = Grad_Descent(_f, _x, _α, _β, _η, _κ)
    s = Search(sol, _f, _ℓ, _γ, _η)
    x_prev = s[1]
    _ℓ = s[2]

    while abs(norm(_x) - norm(x_prev)) > _η

        s = Search(sol, _f, _ℓ, _γ, _η)
        x_prev = s[1]
        _ℓ = s[2]

        sol = Grad_Descent(_f, x_prev, _α, _β, _η, _κ)
        s = Search(sol, _f, _ℓ, _γ, _η)
        _x = s[1]
        _ℓ = s[2]
    end

    sol = Newton(_f, _x, _ϵ)

    return sol
end



#f(x) = x^2 + sin(3x^2)
#f(x) = x^4 + 3x^3 + x^2 + x + sin(3x^4)
#f(x) = x^4 - x^3 - 50x^2 + 100*sin(30x)
#f(x) = x^2/20 + 10sin(x) + 5sin(5x)
f(x) = x^4/20000 - x^3/1000 - x^2/10 + 10sin(x) + 5sin(5x)

x0 = -50.0
ϵ = 1e-8
η = 1e-2
α = 0.5
β = 0.8
κ = 0.01
ℓ = 35
γ = 0.99
searchWidth = 10
xPlot = []
yPlot = []
xSol = []
ySol = []
xSearchLeft = []
ySearchLeft = []
xSearchRight = []
ySearchRight = []
xSolFinal = []
ySolFinal = []

var = x0
maxIterations = 150


minimum = Horizontal_Search(f, x0, α, β, η, ϵ, κ, ℓ, γ, searchWidth, maxIterations)
xlist = range(-50.0, 65.0, length = 1000)
plot( xlist, f.(xlist), #=xrange = (25,50), yrange = (-115,-40),=# legend = false)
scatter!(xPlot, yPlot, markersize = 2)
scatter!(xSol, ySol)
plot!(xSearchRight, ySearchRight)
plot!(xSearchLeft, ySearchLeft)
scatter!(xSolFinal, ySolFinal, color = "Red")
