using Plots

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
        ρ_list[i-1] = (x_list[i] - x_list[i-1])^(1/length(x_list[i]))
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
            R_list[i] = ρ_list[i] + (z_list[i] - z_list[i])^2 / (m^2*ρ_list[i]) - 2*(z_list[i] + z_list[i-1])/m
        else
            R_list[i] = 2*ρ_list[i-1] - 4*z_list[i-1]/m
        end
    end

    return R_list
end


function ArgMaxExcludeFirstAndLast(lst)
    _max = -Inf
    _loc = 0
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
    if t == 1
        return x_list[t], x_list[t + 1], t
    else
        return x_list[t - 1], x_list[t], t
    end
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


function Solve_F(_f, inter, _r, _ϵ, maxIt)
    minY = Inf
    minX = Inf
    for k = 2:maxIt
        𝚇 = SubdivideInterval(inter, k)
        z = Compute_z(_f, 𝚇)
        ρ = Compute_ρ(𝚇)
        m = Compute_holder(𝚇, z, _r, ρ)
        R = Compute_charList(z, ρ, m)
        𝑿 = DetermineMaxChar(R, 𝚇)
        t = 𝑿[3]
        x = Compute_x(𝚇, z, t, _r, m)
        f = _f(x)

        if f < minY
            minY = f
            minX = x
        end


        append!(x_plot, x)
        append!(y_plot, _f(x))

        if t > 1
            println()
            println("x: ", x)
            println("ρ: ", ρ[t - 1])
        end

        if t > 1 && ρ[t-1] < _ϵ
            break
        end
    end

    println("\nSolution: ")
    println("x = ", minX)
    println("f = ", minY)
end

f(x) = (x - 0.5)^2 + 0.03*sin(20x)
interval = [0, 1]
r = 10000
ϵ = 1e-8
maxIterations = 100
x_plot = []
y_plot = []

Solve_F(f, interval, r, ϵ, maxIterations)

plot(f, 0,1, label = "f")
scatter!(x_plot, y_plot, label = "f(x)")


𝚇 = SubdivideInterval(interval, 2)
z = Compute_z(f, 𝚇)
ρ = Compute_ρ(𝚇)
m = Compute_holder(𝚇, z, r, ρ)
R = Compute_charList(z, ρ, m)
𝑿 = DetermineMaxChar(R, 𝚇)
t = 𝑿[3]
x = Compute_x(𝚇, z, t, r, m)
println()
println(𝑿[1])
println(𝑿[2])
