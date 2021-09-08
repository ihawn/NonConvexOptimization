using Plots


function IntPoint(_f, x1, x2, _L)
    intX = (f(x2) - f(x1))/(2L) + (x1 + x2)/2
    return intX, f(intX), f(x1) + L*(intX - x1)
end

function PiyavskiiSolve(_f, _L, _a, _b, _ϵ, maxIt, plot)
    X = []
    append!(X, _a)
    append!(X, (_a+_b)/2)
    append!(X, _b)

    F = []
    append!(F, _f(X[1]))
    append!(F, _f(X[2]))
    append!(F, _f(X[3]))

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

    append!(plotX, X)
    append!(plotY, F)

    for n = 1:maxIt

        if n == 1
            maxPos = 2
        else
            maxPos = argmax(Q)#maxPos = argmax(F)
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
        α = maximum(F)
        β = (_f(X[maxPos - 1]) + _f(X[maxPos + 1]))/2 - _L*(X[maxPos + 1] - X[maxPos-1])/2

        if plot
            append!(plotX, sawPoint1[1])
            append!(plotY, sawPoint1[2])
            append!(plotX, sawPoint2[1])
            append!(plotY, sawPoint2[2])
            push!(linePlotX, [X[maxPos-1], sawPoint1[1]], [sawPoint1[1], X[maxPos+1]])
            push!(linePlotX, [X[maxPos+3], sawPoint2[1]], [sawPoint2[1], X[maxPos+1]])
            push!(linePlotY, [F[maxPos-1], sawPoint1[3]], [sawPoint1[3], F[maxPos+1]])
            push!(linePlotY, [F[maxPos+3], sawPoint2[3]], [sawPoint2[3], F[maxPos+1]])
            append!(plotα, α)
            append!(plotβ, β)
        end

        if abs(α - β) <= _ϵ
            break
        end
    end

    z = (α + β)/2
    @show z

    return [z, plotX, plotY, linePlotX, linePlotY, plotα, plotβ]
end


f(x) = -sqrt(x^2 + 5) - sin(3x)#-(x + 5sin(x))*3^(-x^2)##
L = 4#1.12677256487
ϵ = 1e-8
a = -8.0
b = 7.0


p = PiyavskiiSolve(f, L, a, b, ϵ, 200, true)



plot(f, a, b, ylims = (-15,15), xlims = (a, b), legend = false)
plot!(p[4][:], p[5][:], color=:green)
scatter!(p[2], p[3])

plot(p[6], label = "α")
plot!(p[7], label = "β")
