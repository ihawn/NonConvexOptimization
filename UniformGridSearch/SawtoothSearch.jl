using Plots

SawTooth(x1, x2, z1, z2, L) = (z1 - z2)/2L + (x1 + x2)/2, (z1 + z2)/2 + L*(x1 - x2)/2

function Piyavskii(f, a, b, L, ϵ, m)
    fa = f(a); fb = f(b)
    x, y = SawTooth(a, b, fa, fb, L)

    x_list = [a, x, b]
    y_list = [fa, y, fb]
    X = [a, b]
    Y = [fa, fb]
    α_list = Float64[]
    β_list = Float64[]
    minpos = 0

    for i = 1:m
        minpos = argmin(y_list)
        y_list[minpos] = f(x_list[minpos])

        for j = 1:2
            pos1 = minpos - j%2
            pos2 = minpos + (j+1)%2
            x, y = SawTooth(x_list[pos1], x_list[pos2], y_list[pos1], y_list[pos2], L)
            append!(x_list, x); append!(y_list, y);
            append!(X, x); append!(Y, f(x));
        end

        p = sortperm(x_list); permute!(x_list, p); permute!(y_list, p)

        β = minimum(y_list)
        α, minpos = findmin(Y)
        append!(α_list, α); append!(β_list, β)

        if abs(α - β) <= ϵ
            break
        end
    end

    x = X[minpos]
    y = Y[minpos]
    @show x
    @show y

    if abs(α_list[end] - β_list[end]) > ϵ
        @warn("Failed to converge.")
    end

    return x, y, x_list, y_list, α_list, β_list
end


# f(x) = -sum(k*sin((k+1)*x + k) for k in 1:6)
# L = 109.321
# a, b = -8.0, 8.0

f(x) = sqrt(x^2 + 5) + sin(3x)
L = 4.0
a, b = -8.0, 8.0

# f(x) = x^2
# L = 4
# a, b= -2.0, 2.0

ϵ = 1e-8

x, y, x_list, y_list, α_list, β_list = Piyavskii(f, a, b, L, ϵ, 2)

p=plot(f, xlims = (a, b), dpi = 200, legend = false, title = "Piyavskii's Algorithm: i=1")
plot!(x_list, y_list)
scatter!([x], [y])
# show(p)
# savefig(p, "~/NonConvexOptimiztion/UniformGridSearch/piyavskii-i-1.png")

# p = plot(α_list, color=:black, dpi=200, title = "α, β Convergence", legend = false, xaxis=:log)
# plot!(β_list, color=:black)
# #hline!([0], color=:black)
# savefig(p, "~/NonConvexOptimiztion/UniformGridSearch/piyavskiiαβ.png")
