using LinearAlgebra
include("C:/Users/Isaac/Documents/Optimization/NonConvex/NonConvexOptimization/NonConvexOptimiztion/testobjectives.jl")

function Σ(x1, x2)
    body
end

function Noisy_F(_f, _x, _noise)
    return _f(_x) + (rand() - 0.5)*_noise
end

bounds = [-1.0, 2.0]
f(x) = -sin(3x) - x^2 + 0.7x
noise = 0.2
x0 = [-0.9, 1.1] #first two initial points
y0 = Noisy_F.(f, x0, noise)


noiseX = range(bounds[1], bounds[2], step = 0.01)
noiseY = Noisy_F.(f, noiseX, noise)

plot(f, xlims = bounds, label = "f(x)")
scatter!(noiseX, noiseY, label = "Noise", markersize = 2)
scatter!(x0, y0, label = "Initial Sample")
