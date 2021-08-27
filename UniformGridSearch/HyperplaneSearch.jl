using Plots
using LinearAlgebra

function GenerateHyperplane(_x, _fx, _L, _θ, sign, switch)
    a2 = _L/sqrt(1 + sin(_θ)^2) * sign
    a1 = -a2*sin(_θ)
    k = _fx - a1*_x[1] - a2*_x[2] #gives a1*x1 + a2*x2 + k = x3

    if switch
        return [a2 a1 -1 -k]
    else
        return [a1 a2 -1 -k]
    end
end

function Intersect(h)
    n = length(h[1])
    A = zeros(Float64, n-1, n-1)
    b = zeros(n-1)
    for i = 1:n-1
        for j = 1:n-1
            A[i,j] = h[i][j]
        end
    end
    for i = 1:n-1
        b[i] = h[i][4]
    end

    return A\b
end

function ShrinkBounds(Ω, _f, _L)


    for i = 1:length(Ω)
        append!(boxplot_x, Ω[i][1])
        append!(boxplot_y, Ω[i][2])
    end

    append!(urx, Ω[1][1])
    append!(ury, Ω[1][2])
    append!(lrx, Ω[2][1])
    append!(lry, Ω[2][2])
    append!(llx, Ω[3][1])
    append!(lly, Ω[3][2])
    append!(ulx, Ω[4][1])
    append!(uly, Ω[4][2])


    α = sum(Ω)/length(Ω)


    h = [] #boundary hyperplanes
    p = [] #midpoint hyperplanes
    q = [] #intersection points


    push!(h, GenerateHyperplane(Ω[1], _f(Ω[1]), _L, 2*atan((α[1] - Ω[1][1])/(α[2] - Ω[1][2])), -1, false))
    push!(h, GenerateHyperplane(Ω[2], _f(Ω[2]), _L, 2*atan((α[1] - Ω[2][1])/(α[2] - Ω[2][2])), -1, false))
    push!(h, GenerateHyperplane(Ω[3], _f(Ω[3]), _L, 2*atan((α[1] - Ω[3][1])/(α[2] - Ω[3][2])), 1, false))
    push!(h, GenerateHyperplane(Ω[4], _f(Ω[4]), _L, 2*atan((α[1] - Ω[4][1])/(α[2] - Ω[4][2])), 1, false))

    push!(p, GenerateHyperplane(α, _f(α), _L, 0, 1, false))
    push!(p, GenerateHyperplane(α, _f(α), _L, pi, -1, true))
    push!(p, GenerateHyperplane(α, _f(α), _L, -2*pi, -1, false))
    push!(p, GenerateHyperplane(α, _f(α), _L, -3*pi, 1, true))

    # println()
    # for i = 1:4
    #     println(h[i])
    # end
    # println()
    # for i = 1:4
    #     println(p[i])
    # end
    #

    for i = 1:length(h)
        #push!(p, GenerateHyperplane(α, _f(α), _L, (i-1)*pi/2, 1))

        offset = 1
        if i == length(h) #Wrap around on the last iteration to find intersect of h4, h1, f4
            offset = -length(h) + 1
        end
        push!(q, Intersect([h[i], h[i + offset], p[i]]))
    end


    vals = zeros(Float64, length(q))
    for i = 1:length(q)
        vals[i] = q[i][length(q[1])] #_f(q[i][1:2])
    end
    new_corner = argmax(vals) #find maximum intersection point


    #Ensure rectangular bounds by implementing wrap around
    right = 1
    left = -1
    if new_corner == length(Ω)
        right = -length(Ω) + 1
    end
    if new_corner == 1
        left = length(Ω) - 1
    end

    Ω[new_corner] = q[new_corner][1:2]
    Ω[new_corner + left][1] = q[new_corner][1]
    Ω[new_corner + right][2] = q[new_corner][2]

    append!(aveplot_x, α[1])
    append!(aveplot_y, α[2])
    println(norm(Ω[1] - Ω[3]))
    return Ω
end

f(x) = -0.03*(x[1]^2 + x[2]^2) + 4*sin(0.1*x[1]*x[2]) + cos(x[1] - 2)
L = 10
bounds = [[9.5, 8], [-9, 8], [-9, -9], [9.5, -9]]
boxplot_x = []
boxplot_y = []
aveplot_x = []
aveplot_y = []

urx = []
ury = []
lrx = []
lry = []
llx = []
lly = []
ulx = []
uly = []

h1 = GenerateHyperplane(bounds[1], f(bounds[1]), L, -pi/2, -1)
h2 = GenerateHyperplane(bounds[2], f(bounds[2]), L, pi/2, -1)
h3 = GenerateHyperplane([0, 0], f([0, 0]), L, 0, 1)


for i = 1:4
    bounds = ShrinkBounds(bounds, f, L)
end



zoomMult = 1
plotf(x,y) = f([x, y])
_x = -10:0.1:10
_y = -10:0.1:10
X = repeat(reshape(_x, 1, :), length(_y), 1)
Y = repeat(_y, 1, length(_x))
Z = map(plotf, X, Y)
p1 = Plots.contour(_x,_y, plotf, fill = true, aspect_ratio=:equal)
plot(p1, xrange = (-10, 10), yrange = (-10, 10))
scatter!(boxplot_x, boxplot_y, legend = false, color =:green)
plot!(aveplot_x, aveplot_y, legend = false, color =:black)

plot!(urx, ury, legend = false, color =:white)
plot!(lrx, lry, legend = false, color =:white)
plot!(llx, lly, legend = false, color =:white)
plot!(ulx, uly, legend = false, color =:white)
