using Plots

f(x) = sqrt(x^2 + 5) + sin(3x)
L = 4
κ = 1024
ϵ = 1e-5
a = -2.0
b = 2.0

αPlot = Float64[]
βPlot = Float64[]
gap = Float64[]

x = Float64[]
z = Float64[]
zlow = Float64[]

axis = Float64[]

k = 4
while k < κ

    append!(axis, k)
    k *= 2

    x = zeros(k+1)
    z = zeros(k+1)

    for i = 1:k+1
        # _x = a + (i-1)/k * (b - a)
        # append!(x, _x)
        # append!(z,f(_x))
        x[i] = a + (i-1)/k * (b - a)
        z[i] = f(x[i])
    end

    zlow = zeros(k)
    for i = 1:k
        #append!(zlow, (z[i] + z[i+1])/2 - L*(x[i+1] - x[i])/2)
        zlow[i] = (z[i] + z[i+1])/2 - L*(x[i+1] - x[i])/2
    end

    α = minimum(z)
    β = minimum(zlow)

    append!(αPlot, α)
    append!(βPlot, β)
    append!(gap, α - β)

    println("\nIteration ", k)
    println("f* = ", α)

    if abs(α - β) <= ϵ
        break
    end
end

# plot(f, a, b)
# scatter!(x, z)
# scatter!(x[1:κ], zlow)
plot(axis,αPlot, label = "α")
plot!(axis,βPlot, label = "β")
# plot(gap, yaxis=:log)
