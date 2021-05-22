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


function Generate_ℓ_Vector(_x, n, _ℓ)
    vec = zeros(length(_x))
    vec[n] = _ℓ
    return  vec
end


function Search(local_sol, _f, _ℓ, _γ, _η)

    ℓ_start = _ℓ

    for i = 1:length(local_sol[1])
        vec = Generate_ℓ_Vector(_x, i, _ℓ)
        while _f(local_sol[1][i] + vec) >= local_sol[2] + _η && _f(local_sol[1][i] + vec) >= local_sol[2] + _η && _ℓ > _η
            _ℓ *= _γ
        end
    end



    if _f(local_sol[1] + _ℓ) < local_sol[2] + _η
        return local_sol[1] + _ℓ, _ℓ
    else
        return local_sol[1] - _ℓ, _ℓ
    end
end


function Horizontal_Search(_f, _x, _α, _β, _η, _ϵ, _κ, _ℓ, _γ, width, maxIt)

    sol = Grad_Descent(_f, _x, _α, _β, _η, _κ)
#    s = Search(sol, _f, _ℓ, _γ, _η)
    #x_prev = s[1]
    #_ℓ = s[2]

#=    while abs(norm(_x) - norm(x_prev)) > _η

        s = Search(sol, _f, _ℓ, _γ, _η)
        x_prev = s[1]
        _ℓ = s[2]

        sol = Grad_Descent(_f, x_prev, _α, _β, _η, _κ)
        s = Search(sol, _f, _ℓ, _γ, _η)
        _x = s[1]
        _ℓ = s[2]
    end

    sol = Newton(_f, _x, _ϵ)=#

    return sol
end



x0 = [2,3]
ϵ = 1e-8
η = 1e-2
α = 0.5
β = 0.8
κ = 0.5
ℓ = 35
γ = 0.99
searchWidth = 10

xPlot = []
yPlot = []
solPlotX = []
solPlotY = []

var = x0
maxIterations = 150

f(x) = x[1]^2 + x[2]^2 + 5*sin(x[1] + x[2])

minimum = Horizontal_Search(f, x0, α, β, η, ϵ, κ, ℓ, γ, searchWidth, maxIterations)

plotf(x,y) = f([x, y])
_x = -4.0:0.03:4.0
_y = -4.0:0.03:4.0
X = repeat(reshape(_x, 1, :), length(_y), 1)
Y = repeat(_y, 1, length(_x))
Z = map(plotf, X, Y)
p1 = Plots.contour(_x,_y, plotf, fill = true)
plot(p1, legend = false)
plot!(xPlot, yPlot, color = "white")
scatter!(xPlot, yPlot, markersize = 2, color = "red")
scatter!(solPlotX, solPlotY, color = "green")
