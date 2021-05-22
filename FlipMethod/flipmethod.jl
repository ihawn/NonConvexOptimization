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

function Grad_Descent(_f, orig_f, _x, _α, _β, _ϵ, _κ)
    it = 0;
    t = 1
    ∇f = Compute_Gradient(_f, _x)
    normGrad = norm(∇f)

    push!(xPlot, _x)
    push!(yPlot, orig_f(_x))

    val = _f(_x)

    while normGrad > _ϵ
        ∇f = Compute_Gradient(_f, _x)
        normGrad = norm(∇f)

        bls = Back_Line_Search(_x, _f, val, ∇f, -∇f, _α, _β, _κ)
        t = bls[1]
        _x -= t*∇f
        val = _f(_x)

        push!(xPlot, _x)
        push!(yPlot, orig_f(_x))



    end
    push!(xSol, _x)
    push!(ySol, orig_f(_x))

    return _x, val
end

function Normal_Sign(a)
    if a < 0
        return  -1.0
    elseif a > 0
        return  1.0
    else
        return  0.0
    end
end

function Flip_Method(_f, _x, _α, _β, _η, _ϵ, _κ, maxIt, minVal, maxVal)
    neg_f(x) = -f(x)

    x_start = _x
    sol = Grad_Descent(_f, _f, _x, _α, _β, _η, _κ)
    min = sol[2]
    minX = sol[1]

    _x = sol[1] + 2*_η*Normal_Sign(sol[1] - x_start)
    sol = Grad_Descent(neg_f, _f, _x, _α, _β, _η, _κ)

    global keepGoing = true
    i = 0

    while keepGoing && i < maxIt


        _x = sol[1] + 2*_η*Normal_Sign(sol[1] - x_start)
        x_start = sol[1]
        sol = Grad_Descent(_f, _f, _x, _α, _β, _η, _κ)
        val = _f(_x)

        if sol[2] < min
            min = sol[2]
            minX = sol[1]
        else
            keepGoing = false
        end

        if _x < minVal || _x > maxVal
            keepGoing = false
        end

        _x = sol[1] + 2*_η*Normal_Sign(sol[1] - x_start)
        sol = Grad_Descent(neg_f, _f, _x, _α, _β, _η, _κ)

        print("\nf(x) = ", sol[2])

        i += 1
    end

    sol = Newton(_f, minX, _ϵ)

    push!(xSolFinal, sol[1])
    push!(ySolFinal, sol[2])

    print("\nSolution found at x = ", sol[1], "\nf(x) = ", sol[2])

    return sol
end

f(x) = x^2 + sin(3x^2)
#f(x) = x^4 + 3x^3 + x^2 + x + sin(3x^4)
#f(x) = x^4 + x^3 - x^2


x0 = 5.0
ϵ = 1e-8
η = 1e-2
α = 0.5
β = 0.8
κ = 0.01
xPlot = []
yPlot = []
xSol = []
ySol = []
xSolFinal = []
ySolFinal = []

var = x0
maxIterations = 150

xlist = range(-3, 5, length = 1000)
minimum = Flip_Method(f, x0, α, β, η, ϵ, κ, maxIterations, -1e10, 1e10)
plot( xlist, f.(xlist), legend = false)
scatter!(xPlot, yPlot, markersize = 2)
scatter!(xSol, ySol)
scatter!(xSolFinal, ySolFinal, color = "Red")
