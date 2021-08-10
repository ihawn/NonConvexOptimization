using Zygote
using LinearAlgebra
using Plots
using DataStructures
include("C:/Users/Isaac/Documents/Optimization/NonConvex/NonConvexOptimization/NonConvexOptimiztion/testobjectives.jl")

function SetB(_x, _fx, _∇f, _∇2f, _K, _v, _nv)
    B = [_fx + _∇2f*(_v[1] - _x) + _K*norm(_v[1] - _x)^2]
    for k = 2:_nv
        append!(B, _fx + _∇2f*(_v[k] - _x) + _K*norm(_v[k] - _x)^2)
    end
    return B
end

function SetAdj(_v, _nv)
    _adj = [v[1]]
    for k = 2:_nv
        append!(_adj, _v[k])
    end
    return _adj
end

f(x) = x[1]^2 + x[2]^2
x = [3, 3]
v = [[-10, -10], [10, -10], [10, 10], [-10, 10]]
nv = length(v)
K = 10


#Initialize
fx = f(x)
∇f = f'(x)
∇2f = Zygote.hessian(f, x)
B = SetB(x, fx, ∇f, ∇2f, K, v, nv)
adj = SetAdj(v, nv)
h = BinaryMinHeap(v)
maxIt = 1000

for i = 1:maxIt


    for k = 1:nv

    end
end
