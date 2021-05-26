using Plots
using Calculus
using LinearAlgebra

function Compute_Gradient(_f, _x)
    return Calculus.gradient(g -> _f(g),_x)
end


function Compute_Hessian(_f, _x)
    return hessian(h -> _f(h),_x)
end


function m(_p, _fx, _∇f, _∇2f)
    return _fx + transpose(_∇f)*_p + transpose(_p)*_∇2f*_p/2.0
end


function Rho(_fx, _f, _x, _p, _m, _∇f, _∇2f)
    return (_fx - _f(_x + _p))/(_m(zeros(length(_x)), _fx, _∇f, _∇2f) - _m(_p, _fx, _∇f, _∇2f))
end


function Subproblem_Cauchy_Point(_Δk, _∇f, _∇2f)

    τ = 1.0
    nrm_Δf = norm(_∇f)
    val = transpose(_∇f) * _∇2f * _∇f

    if val > 0
        τ = min(nrm_Δf^3 / (_Δk * val), 1.0)
    end

    return -τ * _Δk/nrm_Δf * _∇f
end


function Conj_Grad_Steihaug(_Δk, _∇f, _∇2f, _ϵ, _itt)
    z = 0
    r = _∇f
    d = -r

    if norm(r) < _ϵ
        return 0
    end

    for j = 0:_itt
        if transpose(d)*_∇2f*d <= 0

        end
    end

end



function Trust_Region(_f, _x, _Δk, _Δm, _η1, _η2, _η3, _t1, _t2, _ϵ, _δ, itt)

    println("\n\n")

    for k = 1:itt
        fx = _f(_x)
        ∇f = Compute_Gradient(_f, _x)
        ∇2f = Compute_Hessian(_f, _x)


        p = Subproblem_Cauchy_Point(_Δk, ∇f, ∇2f) #Solve trust region subproblem

        ρ = Rho(fx, _f, _x, p, m, ∇f, ∇2f)

        if ρ < _η2
            _Δk *= _t1
        elseif ρ > _η3 && abs(norm(p) - _Δk) <= _δ
            _Δk = min(_t2*_Δk, _Δm)
        else
            #do nothing i.e. _Δk remains the same
        end

        if ρ > _η1
            _x += p
        else
            #_x stays the same i.e. model is poor and we need to solve another subproblem
        end

        println("Iteration: ", k)
        println("x = ", _x)
        println("f(x) = ", _f(_x))

        push!(xPlot, _x[1])
        push!(yPlot, _x[2])

        if norm(∇f) <= _ϵ
            break
        end
    end

    return _x
end

flush(stdout)

f(x) = x[1]^2 + x[2]^2 + 7*sin(x[1] + x[2]) + 10*sin(5x[1])
#f(x) = (x[2] - 0.129*x[1]^2 + 1.6*x[1] - 6)^2 + 6.07*cos(x[1]) + 10

η1 = 0.2
η2 = 0.25
η3 = 0.75
t1 = 0.25
t2 = 2.0
x0 = [10.0, 10.0]
Δk = 2.0
Δm = 4.0
ϵ = 1e-2
δ = 1e-3
maxIterations = 1e3


xPlot = []
yPlot = []


Trust_Region(f, x0, Δk, Δm, η1, η2, η3, t1, t2, ϵ, δ, maxIterations)


plotf(x,y) = f([x, y])
# _x = -5.0:0.03:10.0
# _y = 0.0:0.03:15.0
_x = -10.0:0.03:10.0
_y = -10.0:0.03:10.0
X = repeat(reshape(_x, 1, :), length(_y), 1)
Y = repeat(_y, 1, length(_x))
Z = map(plotf, X, Y)
p1 = Plots.contour(_x,_y, plotf, fill = true)
plot(p1, legend = false, title = "Global Minimization With Trust Regions")
plot!(xPlot, yPlot, color = "white")
scatter!(xPlot, yPlot, color = "red", markersize = 2)
