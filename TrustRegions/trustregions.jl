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


function Subproblem_Cauchy_Point(_Δk, _∇f, _∇2f, _nrm_Δk)

    τ = 1.0
    val = transpose(_∇f) * _∇2f * _∇f

    if val > 0
        τ = min(norm(_∇f)^3 / (_Δk * val), 1.0)
    end

    return -τ * _Δk/norm(_∇f) * _∇f
end


function Trust_Region(_f, _x, _Δk, _Δm, _η1, _η2, _η3, _t1, _t2, _ϵ, itt)

    for k = 1:itt
        fx = _f(_x)
        ∇f = Compute_Gradient(_f, _x)
        ∇2f = Compute_Hessian(_f, _x)
        nrm_Δk = norm(_Δk)

        p = Subproblem_Cauchy_Point(_Δk, _∇f, _∇2f, nrm_Δk) #Solve trust region subproblem

        ρ = Rho(fx, _f, _x, _p, m, ∇f, ∇2f)

        if ρ < _η2
            _Δk *= _t1
        elseif ρ > _η3 && abs(p - nrm_Δk) <= _ϵ
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
        println("f(x) = ", _x)

    end

end
