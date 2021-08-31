using Plots

function IntPoint(f, x1, x2, L)
    intX = (f(x1) - f(x2))/(2L) + (x1 + x2)/2
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

    α = 0.0
    β = 0.0

    plotX = []
    plotY = []

    for n = 1:maxIt
        maxPos = argmax(F)
        maxX = X[maxPos]

        if maxPos == 1
            maxPos += 1
        elseif maxPos == length(X)
            maxPos -= 1
        end

        sawPoint = IntPoint(_f, X[maxPos - 1], X[maxPos+1], _L)

        append!(X, sawPoint[1])
        append!(F, sawPoint[2])

        p = sortperm(X); X .= X[p]; F .= F[p]
        α = sawPoint[2]
        β = (_f(X[maxPos - 1]) + _f(X[maxPos + 1]))/2 - _L*(X[maxPos + 1] - X[maxPos-1])/2

        if plot
            append!(plotX, sawPoint[1])
            append!(plotY, sawPoint[2])
        end

        if abs(α - β) <= _ϵ
            break
        end
    end

    return [(α + β)/2, plotX, plotY]
end


f(x) = -(x + sin(x))*3^(-x^2)#-sqrt(x^2 + 5) - sin(3x)
L = 200
ϵ = 1e-8
a = -8.0
b = 8.0


p = PiyavskiiSolve(f, L, a, b, ϵ, 10000, true)

plot(f)
scatter!(p[2], p[3])
