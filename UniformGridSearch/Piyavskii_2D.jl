include("Piyavskii.jl")

f(_x) = -0.03*(_x[1]^2 + _x[2]^2) + 4*sin(0.1*_x[1]*_x[2]) + cos(x[1] - 3)
L = 3.97
ϵ = 1e-8
a = [-8.0, 0.0]
b = [8.0, 0.0]
x = sum([a, b])/2
plot_x = []
plot_y = []

_f(_x) = _x[1] + _x[2]

for i = 1:10
    #fix x
    p = 0
    if i%2 == 1
        _f(_x) = f([x[1], _x])
        p = PiyavskiiSolve(_f, L, a[2], b[2], ϵ, 100)
        println(p)
        x[2] = p[1]
    end

    #fix y
    if i%2 == 0
        _f(_x) = f([_x, x[2]])
        p = PiyavskiiSolve(_f, L, a[1], b[1], ϵ, 100)
        x[1] = p[1]
    end

    append!(plot_x, x[1])
    append!(plot_y, x[2])
end


plotf(x,y) = f([x, y])
_x = -10:0.1:10
_y = -10:0.1:10
X = repeat(reshape(_x, 1, :), length(_y), 1)
Y = repeat(_y, 1, length(_x))
Z = map(plotf, X, Y)
p1 = Plots.contour(_x,_y, plotf, fill = true, aspect_ratio=:equal)
plot(p1, xrange = (-10, 10), yrange = (-10, 10))
scatter!(plot_x, plot_y, color = :green)
