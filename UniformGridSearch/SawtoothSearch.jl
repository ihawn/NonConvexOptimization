using Plots

SawTooth(x1, x2, z1, z2, L) = (z1 - z2)/2L + (x1 + x2)/2, (z1 + z2)/2 + L*(x1 - x2)/2

function Piyavskii(f, a, b, L, ϵ, m)
    fa = f(a); fb = f(b)
    x, y = SawTooth(a, b, fa, fb, L)

    x_list = [a, x, b]
    y_list = [fa, y, fb]
    X = [a, b]
    Y = [fa, fb]
    α_list = []
    β_list = []

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

        p = sortperm(x_list); x_list = x_list[p]; y_list = y_list[p]


        β = (f(x_list[minpos - 1]) + f(x_list[minpos + 1]))/2 - L*(x_list[minpos - 1] - x_list[minpos - 1])/2
        α, minpos = findmin(Y)
        append!(α_list, α); append!(β_list, β)

        if abs(α - β) <= ϵ
            x = X[minpos]
            y = Y[minpos]
            @show x
            @show y
            return x, y, x_list, y_list, α_list, β_list
            break
        end
    end
end


f(x) = -sum(k*sin((k+1)*x + k) for k in 1:6)
L = 109.321
a, b = -8.0, 8.0
ϵ = 1e-8

x, y, x_list, y_list, α_list, β_list = Piyavskii(f, a, b, L, ϵ, 10000)

plot(f, xlims = (a, b), legend=false)
plot!(x_list, y_list)
scatter!([x], [y])

# plot(α_list, legend = false)
# plot!(β_list)
