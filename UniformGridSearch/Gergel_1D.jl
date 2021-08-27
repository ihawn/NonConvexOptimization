using Plots
using Symbolics

function SubdivideInterval(_interval, n)
    x_list = zeros(n + 1)
    x_list[1] = _interval[1]
    dif = _interval[2] - _interval[1]

    for i = 2:n
        x_list[i] = _interval[1] + dif*((i - 1)/n)
    end

    x_list[n + 1] = _interval[2]

    return x_list
end


function Compute_z(_f, x_list)
    n = length(x_list)
    z_list = zeros(n)
    for i = 1:n
        z_list[i] = _f(x_list[i])
    end

    return z_list
end


function Compute_ρ(x_list)
    n = length(x_list)
    ρ_list = zeros(n-1)
    for i = 2:n
        ρ_list[i-1] = (x_list[i] - x_list[i-1])
    end

    return ρ_list
end


function Compute_holder(x_list, z_list, _r, ρ_list)
    n = length(z_list)
    candidate_list = zeros(n - 1)

    for i = 2:n
        candidate_list[i-1] = abs(z_list[i] - z_list[i-1])/ρ_list[i-1]
    end

    M = maximum(candidate_list)

    return _r*M
end


function Compute_charList(z_list, ρ_list, m)
    n = length(z_list)
    R_list = zeros(n)

    for i = 1:n
        if i == 1
            R_list[i] = 2*ρ_list[i] - 4*z_list[i]/m
        elseif i < n
            R_list[i] = ρ_list[i] + (z_list[i] - z_list[i-1])^2 / (m^2*ρ_list[i]) - 2*(z_list[i] + z_list[i-1])/m
        else
            R_list[i] = 2*ρ_list[i-1] - 4*z_list[i-1]/m
        end
    end

    return R_list
end


function ArgMaxExcludeFirstAndLast(lst)
    _max = -Inf
    _loc = 2
    for i = 2:length(lst) - 1
        if lst[i] > _max
            _max= lst[i]
            _loc = i
        end
    end
    return _loc
end


function DetermineMaxChar(R_list, x_list)
    t = ArgMaxExcludeFirstAndLast(R_list)
    return x_list[t - 1], x_list[t], t
end


function Compute_x(x_list, z_list, t, _r, _m)
    n = length(x_list)
    N = length(x_list[1])
    _x = x_list[t]
    if t == 2 || t == n
        _x = (x_list[t] + x_list[t - 1])/2
    elseif t != 1
        _x -= sign(z_list[t] - z_list[t - 1])/(2*_r) * (r*abs(z_list[t] - z_list[t - 1])/_m)^N
    end

    return _x
end


function PerformTrials(_f, _𝚇, _z, _t, _r, _m, _k)
    trialPoints_x = zeros(_k)
    trialPoints_z = zeros(_k)

    for i = 1:_k
        trialPoints_x[i] = Compute_x(_𝚇, _z, _t, _r, _m)
        trialPoints_z[i] = f(trialPoints_x[i])
    end

    j = argmin(trialPoints_z)

    return trialPoints_x[j], trialPoints_z[j]
end


function Solve_F(_f, inter, _r, _ϵ, maxIt, start)
    minY = Inf
    minX = Inf
    solLoc = 0
    startInter = inter

    for k = start:maxIt
        append!(interval_list, inter[1])
        append!(interval_list, inter[2])

        𝚇 = SubdivideInterval(inter, k)
        z = Compute_z(_f, 𝚇)
        ρ = Compute_ρ(𝚇)
        m = Compute_holder(𝚇, z, _r, ρ)
        R = Compute_charList(z, ρ, m)
        𝑿 = DetermineMaxChar(R, 𝚇)
        t = 𝑿[3]
        x = Compute_x(𝚇, z, t, _r, m)
        diff = 0 #max(x - 𝑿[1], 𝑿[2] - x)
        inter = [𝑿[1] - diff, 𝑿[2] + diff]
        f = _f(x)
        solLoc = k

        if f < minY
            minY = f
            minX = x
        end
        println(t)
        println("Interval: ", inter)

        append!(x_plot, x)
        append!(y_plot, _f(x))

        if ρ[t-1] < _ϵ
            break
        end
    end

    println("\nSolution found at iteration ", solLoc - start)
    println("x = ", minX)
    println("f = ", minY)

    return minY
end

#f(x) = 0.0001x^2 + sin(20x)
#f(x) = -(x + sin(x))*exp(-x^2)
#f(x) = -(1.4 - 3x)*sin(18x)
#f(x) = -sum(k*cos((k + 1)*x + k) for k = 1:5)
#f(x) = -exp(-x)*sin(2*pi*x)
#f(x) = (x^2 - 5x + 6)/(x^2 + 1)
#f(x) = exp(-3x) - (sin(x))^3
f(x) = x^2

interval = [-1.0, 10.1]
#interval = [-10.0, 10.0]
#interval = [0.0, 1.2]
#interval = [-10.0, 1.0]
#interval = [0.0, 4.0]
#interval = [-5.0, 5.0]
#interval = [0.0, 20.0]
#interval = [-1.0, 1.0]

r = 1.1
ϵ = 1e-16
maxIterations = 20
startIteration = 15
x_plot = []
y_plot = []
interval_list = []

# fp(_t) = -abs(substitute(D, t => _t))

Solve_F(f, interval, r, ϵ, maxIterations, startIteration)


plot(f, interval[1],interval[2], label = "f", title = "Gergel's Minimization Algorithm in 1D")
scatter!(x_plot, y_plot, label = "f(x)")
plot!(interval_list, seriestype = :vline, label = "Intervals")
