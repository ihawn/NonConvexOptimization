using Zygote
using LinearAlgebra
using CSV
using DataFrames
using Statistics
include("testobjectives.jl")

# L = 0.0
# f(x) = Himmelblau(x, 2)
# size = 10.0
# grads = []

# for i = 1:1000000
#     pt = rand(2)*2*size .- size
#     grad = norm(f'(pt))
#     append!(grads, grad)
#     if grad > L
#         global L = norm(grad)
#     end
# end

# println("L = ", L)


# funcs = [Griewank, Ackley, Himmelblau, Easom, Rastrigin, Holder_Table, Drop_Wave, Schaffer_N2, Styblinski_Tang]
# names = ["Griewank", "Ackley", "Himmelblau", "Easom", "Rastrigin", "Holder_Table", "Drop_Wave", "Shaffer N. 2", "Styblinski-Tang"]
# _σ = []
# ave = []
# maxS = []


# for i = 1:length(funcs)
#     global L = 0
#     f(x) = funcs[i](x, 2)
#     for i = 1:1000000
#         pt = rand(2)*2*size .- size
#         grad = norm(f'(pt))
#         append!(grads, grad)
#         if grad > L
#             global L = norm(grad)
#         end
#     end

#     append!(_σ, std(grads))
#     append!(ave, mean(grads))
#     append!(maxS, L)
# end

# df = DataFrame(Name = names, stdGrad = _σ, aveGrad = ave, maxGrad = maxS)
# CSV.write("csv/FunctionData.csv", df)

names1 = ["Himmelblau", "Griewank", "Easom", "Drop Wave", "Ackley", "Shaffer N. 2", "Rastrigin"]
Ls = [0.598377077, 1.003935648, 1.085059873, 6.006944355, 9.222083471, 19.6190065, 114.9895065]
totalRuntimes1 = [124.22, 215.89, 89.28, 573.15, 45.37, 3387.2, 562.33]
