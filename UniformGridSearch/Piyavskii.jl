using Plots

maxIt = 100

f(x) = x^2 + sin(3x^2)
L = 100000

ϵ = 1e-8

#step 0
k = 2
a = -5.0
b = 4.0
x = []
z = []
P = []
zLow = []
append!(x, a)
append!(x, b)
append!(z, f(xi) for xi in x)

α = minimum(z)
β = (z[1] + z[2])/2 - L*(x[2] - x[1])/2

for k = 2:100
    #step 1
    xLow = []
    zLow = []
    P = []

    append!(zLow, (z[i] + z[i+1])/2 - L*(x[i+1] - x[i])/2 for i = 1:k-1)

    #step 2
    append!(P, -z[i] for i = 1:k-1)

    #step 3
    t = argmax(P)

    #step 4
    xLow = (x[t] + x[t+1])/2 + (z[t] - z[t+1])/(2L)
    append!(x, xLow)
    append!(z, f(xLow))

    p = sortperm(x); x .= x[p]; z .= z[p]

    β = zLow[t]
    α = min(f(xLow), α)
    println("x = ", α)

    if abs(α - β) <= ϵ
        break
    end
end


plot(f, a, b)
scatter!(x,z)
plot!(x,z)
