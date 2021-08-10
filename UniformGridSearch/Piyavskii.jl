using Plots

function IntPoint(f, x1, x2, L)
    intX = (f(x1) - f(x2))/(2L) + (x1 + x2)/2
    return intX, f(intX), f(x1) + L*(intX - x1)
end

f(x) = -sqrt(x^2 + 5) - sin(3x)
L = 3.97

# f(x) = 12x^3*exp(-x^4) + 20x*exp(-2x^2)
# L = 20

ϵ = 1e-8
a = -8.0
b = 8.0

X = []
append!(X, a)
append!(X, (a+b)/2)
append!(X, b)

F = []
append!(F, f(X[1]))
append!(F, f(X[2]))
append!(F, f(X[3]))

plotX = []
plotY = []

for n = 1:100
    maxPos = argmax(F)
    maxX = X[maxPos]

    sawPoint = IntPoint(f, X[maxPos - 1], X[maxPos+1], L)

    append!(X, sawPoint[1])
    append!(F, sawPoint[2])

    p = sortperm(X); X .= X[p]; F .= F[p]

    if sawPoint[3] < 10
        append!(plotX, sawPoint[1])
        append!(plotY, sawPoint[3])
    end

    α = sawPoint[2]
    β = (f(X[maxPos - 1]) + f(X[maxPos + 1]))/2 - L*(X[maxPos + 1] - X[maxPos-1])/2

    println("\nIteration ", n)
    println("x = ", sawPoint[1])
    println("f(x) = ", sawPoint[2])


    if abs(α - β) <= ϵ
        break
    end
end

append!(plotX, X)
append!(plotY, F)
p = sortperm(plotX); plotX .= plotX[p]; plotY .= plotY[p]

plot(f, a, b, dpi=:250, title = "Piyavskii's Algorithm", color=:greys, label = "f(x)")
plot!(plotX, plotY, label = "Intersection lines",color=:greys)
scatter!(X, F, label = "Iterative Points",color=:greys)
savefig("C:/Users/Isaac/Documents/Optimization/NonConvex/NonConvexOptimization/NonConvexOptimiztion/Plots/Piyavskii.png")
