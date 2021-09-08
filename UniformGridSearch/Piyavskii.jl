using Plots


function SortPlot(x, y, _f)
    n = length(x)
    t = 1
    for i = t:t:n
        j=0
        while i-j-t > 0 && (x[i-j][1] < x[i-j-t][1] || maximum(y[i-j]) < maximum(y[i-j-t]))
            x[i-j], x[i-j-t] = x[i-j-t], x[i-j]
            y[i-j], y[i-j-t] = y[i-j-t], y[i-j]
            j+=t
        end
    end
    return x, y
end

function Dist(x,y)
    return sqrt(sum((x[i]-y[i])^2 for i in 1:length(x)))
end

function RemoveSupersets(x, y)
    n = length(x)

    # for i = 1:2:n-1
    #     for j = 1:2:n-1
    #         if i != j && x[i][1] <= x[j][1] && x[i+1][2] >= x[j+1][2]
    #             x[i] = [0,0]
    #             y[i] = [0,0]
    #             x[i+1] = [0,0]
    #             y[i+1] = [0,0]
    #         end
    #     end
    # end
    for i = 1:n
        for j = 1:n
            t = j%2 + 1
            if i != j && x[i][t] == x[j][t] && y[i][t] < y[j][t]
                x[j] = [0,0]
                y[j] = [0,0]
            end
        end
    end
    return x, y
end

function RemoveSharedDuplicates(x, y)

    i = 1
    while i < length(x)
        j=0
        while i-j < length(x) && sum(x .== x[i-j]) >= 2 && sum(y .== y[i-j]) >= 2
            deleteat!(x, i-j)
            deleteat!(y, i-j)
            j-=1
        end
        i+=1
    end
    return x, y
end

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

function IntPoint(_f, x1, x2, _L)
    intX = (f(x2) - f(x1))/(2L) + (x1 + x2)/2
    return intX, f(intX), f(x1) + L*(intX - x1)
end

function FillZ(_f, _x, _L)
    n = length(_x)-1
    z = []
    for i = 1:n
        append!(z, (_f(_x[i]) + _f(_x[i+1]))/2 - _L*(_x[i+1] - _x[i])/2)
    end
    return z
end

function PiyavskiiSolve(_f, _L, _a, _b, _ϵ, maxIt, κ, plot)
    X = SubdivideInterval([_a, _b], κ)

    F = []
    for i in 1:length(X)
        append!(F, _f(X[i]))
    end

    Z = []
    Q = []
    P = []

    α = 0.0
    β = 0.0

    plotX = []
    plotY = []
    linePlotX = []
    linePlotY = []
    plotα = []
    plotβ = []
    maxPos = Int32(floor(length(X)/2)) + 1

    append!(plotX, X)
    append!(plotY, F)


    for n = 1:maxIt

        if n > 1
            maxPos = argmax(Q)#FillZ(_f, X, _L))
        end
        maxX = X[maxPos]

        sawPoint1 = IntPoint(_f, X[maxPos - 1], X[maxPos], _L)
        sawPoint2 = IntPoint(_f, X[maxPos], X[maxPos+1], _L)

        append!(Q, sawPoint1[3])
        append!(Q, sawPoint2[3])
        append!(P, sawPoint1[1])
        append!(P, sawPoint2[1])


        append!(X, sawPoint1[1])
        append!(X, sawPoint2[1])
        append!(F, sawPoint1[2])
        append!(F, sawPoint2[2])

        p = sortperm(X); X .= X[p]; F .= F[p]
        p = sortperm(P); P .= P[p]; Q .= Q[p]
        α = minimum(Q)
        β = (_f(X[maxPos - 1]) + _f(X[maxPos + 1]))/2 - _L*(X[maxPos + 1] - X[maxPos-1])/2

        if plot
            append!(plotX, sawPoint1[1])
            append!(plotY, sawPoint1[2])
            append!(plotX, sawPoint2[1])
            append!(plotY, sawPoint2[2])

            append!(linePlotX, [X[maxPos-1], sawPoint1[1], X[maxPos+1], X[maxPos+3], sawPoint2[1]])
            append!(linePlotY, [F[maxPos-1], sawPoint1[3], F[maxPos+1], F[maxPos+3], sawPoint2[3]])

            linePlotX, linePlotY = RemoveSharedDuplicates(linePlotX, linePlotY)
            p = sortperm(linePlotX); linePlotX .= linePlotX[p]; linePlotY .= linePlotY[p]
            # deleteat!(linePlotX, maxPos+1)
            # deleteat!(linePlotY, maxPos+1)
            @show maxPos



            # linePlotX, linePlotY = SortPlot(linePlotX, linePlotY)
            # linePlotX, linePlotY = RemoveSupersets(linePlotX, linePlotY)
            # linePlotX[maxPos] = [0,0]
            # linePlotY[maxPos] = [0,0]
            # if maxPos > 1
            #     linePlotX[maxPos-1] = [0,0]
            #     linePlotY[maxPos-1] = [0,0]
            # end



            append!(plotα, α)
            append!(plotβ, β)


        end


        if abs(α - β) <= _ϵ
            #@show n
            break
        end
    end

    z = max(α, β)
    #@show z

    return [z, plotX, plotY, linePlotX, linePlotY, plotα, plotβ]
end


f(x) = -sum(k*sin((k+1)*x + k) for k in 1:6)#-sqrt(x^2 + 5) - sin(3x)
L = 400#1.12677256487
ϵ = 1e-8
a = -8.0
b = 8


p = PiyavskiiSolve(f, L, a, b, ϵ, 500, 2, true)



plot(f, a, b, ylims = (-25,25), xlims = (a, b), legend = false)
plot!(p[4], p[5], color=:green)
#scatter!(p[4], p[5], color=:green)
scatter!(p[2], p[3], markersize=3, color =:orange)
hline!([p[1]], color=:purple, linestyle=:dash)



# plot(p[6], label = "α")
# plot!(p[7], label = "β")
