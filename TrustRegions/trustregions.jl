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
    return _fx + transpose(_∇f)*_p + (1/2)*transpose(_p)*_∇2f*_p
end


function Rho(_fx, _f, _x, _p, _m, _∇f, _∇2f)
    return (_fx - _f(_x + _p))/(_m(0, _fx, _∇f, _∇2f) - _m(_p, _fx, _∇f, _∇2f))
end


function Subproblem_Cauchy_Point(_Δk, _∇f, _∇2f)

    τ = 1.0
    val = transpose(_∇f) * _∇2f * _∇f

    if val > 0
        τ = min(norm(_∇f)^3 / (_Δk * val), 1.0)
    end

    return -τ * _Δk/norm(_∇f) * _∇f
end


function Trust_Region(_f, _x, _Δk, _Δm, _η1, _η2, _η3, _t1, _t2, itt)

    for k = 1:itt
        fx = _f(_x)
        ∇f = Compute_Gradient(_f, _x)
        ∇2f = Compute_Hessian(_f, _x)

        ρ = Rho(fx, _f, _x, _p, m, ∇f, ∇2f)

    end

end
