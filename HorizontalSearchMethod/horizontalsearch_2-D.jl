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


function Search(local_sol, _f, _ℓ, _γ, _η, _κ, _α, _β)

    ℓ_start = _ℓ
    len = length(local_sol)
    low = local_sol[2]

    vec = Generate_ℓ_Vector(local_sol[1], 1, _ℓ)
    push!(searchX, (local_sol[1] + vec)[1])
    push!(searchY, (local_sol[1] + vec)[2])
    push!(searchX, (local_sol[1] - vec)[1])
    push!(searchY, (local_sol[1] - vec)[2])

    vec = Generate_ℓ_Vector(local_sol[1], 2, _ℓ)
    push!(searchX, (local_sol[1] + vec)[1])
    push!(searchY, (local_sol[1] + vec)[2])
    push!(searchX, (local_sol[1] - vec)[1])
    push!(searchY, (local_sol[1] - vec)[2])


    for k in 1:len

        _ℓ = ℓ_start
        vec = Generate_ℓ_Vector(local_sol[1], k, _ℓ)
        low = min(Grad_Descent(_f, local_sol[1] + vec, _α, _β, _η, _κ)[2], Grad_Descent(_f, local_sol[1] - vec, _α, _β, _η, _κ)[2])

        while low >= local_sol[2] + _η && _ℓ > _η
            _ℓ *= _γ
            vec = Generate_ℓ_Vector(local_sol[1], k, _ℓ)
            low = min(Grad_Descent(_f, local_sol[1] + vec, _α, _β, _η, _κ)[2], Grad_Descent(_f, local_sol[1] - vec, _α, _β, _η, _κ)[2])


            push!(searchX, (local_sol[1] + vec)[1])
            push!(searchY, (local_sol[1] + vec)[2])
            push!(searchX, (local_sol[1] - vec)[1])
            push!(searchY, (local_sol[1] - vec)[2])
        end
    end



    if Grad_Descent(_f, local_sol[1] + vec, _α, _β, _η, _κ)[2] < local_sol[2] + _η
        return local_sol[1] + vec, _ℓ
    else
        return local_sol[1] - vec, _ℓ
    end
end


function Horizontal_Search(_f, _x, _α, _β, _η, _ϵ, _κ, _ℓ, _γ, width, maxIt)

    sol = Grad_Descent(_f, _x, _α, _β, _η, _κ)
    s = Search(sol, _f, _ℓ, _γ, _η, _κ, _α, _β)
    x_prev = s[1]
    #_ℓ = s[2]

    while abs(norm(_x) - norm(x_prev)) > _η

        s = Search(sol, _f, _ℓ, _γ, _η, _κ, _α, _β)
        x_prev = s[1]
        #_ℓ = s[2]

        sol = Grad_Descent(_f, x_prev, _α, _β, _η, _κ)
        s = Search(sol, _f, _ℓ, _γ, _η, _κ, _α, _β)
        _x = s[1]
        #_ℓ = s[2]

        println(sol)
    end

#    sol = Newton(_f, _x, _ϵ)
    println(sol)
    push!(finalSolX, sol[1][1])
    push!(finalSolY, sol[1][2])

    return sol
end



x0 = [4,0]
ϵ = 1e-8
η = 1e-2
α = 0.5
β = 0.8
κ = 1
ℓ = 3
γ = 0.5

searchWidth = 10

xPlot = []
yPlot = []
solPlotX = []
solPlotY = []
searchX = []
searchY = []
finalSolX = []
finalSolY = []

var = x0
maxIterations = 150

f(x) = x[1]^2 + x[2]^2 + 7*sin(x[1] + x[2]) + 10sin(5x[1])

minimum = Horizontal_Search(f, x0, α, β, η, ϵ, κ, ℓ, γ, searchWidth, maxIterations)

plotf(x,y) = f([x, y])
_x = -5.0:0.03:7.0
_y = -5.0:0.03:5.0
X = repeat(reshape(_x, 1, :), length(_y), 1)
Y = repeat(_y, 1, length(_x))
Z = map(plotf, X, Y)
p1 = Plots.contour(_x,_y, plotf, fill = true)
plot(p1, legend = false)
scatter!(searchX, searchY, markersize = 2.5, color = "pink")
plot!(xPlot, yPlot, color = "white")
scatter!(xPlot, yPlot, markersize = 2, color = "orange")
scatter!(solPlotX, solPlotY, color = "green")
scatter!(finalSolX, finalSolY, color = "red")
