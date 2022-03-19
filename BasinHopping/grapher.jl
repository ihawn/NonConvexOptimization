include("testobjectives.jl")
using Plots
using CSV
using DataFrames
using StatsPlots
using LinearAlgebra

function GraphObjectiveConvergence(f, xx, yy, size, finalSolX, finalSolY, filename, method, name)
    plotf(x,y) = f([x, y])
    _x = -size:0.005*abs(size):size
    _y = -size:0.005*abs(size):size
    X = repeat(reshape(_x, 1, :), length(_y), 1)
    Y = repeat(_y, 1, length(_x))
    Z = map(plotf, X, Y)
   
    p1 = Plots.contour(_x,_y, plotf, fill = true, contour_labels = false, color = cgrad([:black,:gray,:white]), aspect_ratio=:equal, dx=10, dy=10)  
    plt = plot(p1, xrange = (-size, size), yrange = (-size, size), legendfontsize = 4, dpi = 400, title = name*" Function - "*method)
    plot!(xx, yy, color = "white", linewidth = 1, label = false)
    scatter!(xx, yy, color = "white", markersize = 2.5, label = "Best Solutions")
    scatter!([xx[length(xx)]], [yy[length(yy)]], color = "red", markersize = 3, label = "Final Solution")
    savefig(plt, filename);
end

function PlotRuntime(nn, tt, filename, ttl)
    plt = plot(nn, tt, title = ttl, label = "Dimensional Runtime", legend=:topleft, xlabel = "n", ylabel = "Runtime (seconds)", color = "black", dpi = 150)
    savefig(plt, filename)
end

function GetDataFromCSV(filename)
    xx = []
    yy = []
    data = CSV.File(filename)

    for row in data
        append!(xx, row.x)
        append!(yy, row.y)
    end

    return xx, yy
end

mx = 8.006562530948429
my = 7.853981614754402

hx = -2.805118
hy = 3.131312

f(x) = Rastrigin(x, 2)

name = "Rastrigin"
method = "Basin Hopping"
xx, yy = GetDataFromCSV("csv/basinHopping-"*name*".csv")
GraphObjectiveConvergence(f, xx, yy, 10, 0, 0, "png/"*method*name*".png", method, name)


# nn1, tt1 = GetDataFromCSV("csv/runtime-Ackley.csv")
# nn2, tt2 = GetDataFromCSV("csv/runtime-Rast.csv")
# nn3, tt3 = GetDataFromCSV("csv/runtime-Rosenbrock.csv")

# plt = plot(nn1, title = "Basin Hopping Dimensional Runtime", label = "Ackley Runtime", legend=:topleft, xlabel = "n", ylabel = "Runtime (seconds)", color = "black", dpi = 150)
# plot!(nn3, tt3, label = "Rosenbrock Runtime", legend=:topleft, xlabel = "n", ylabel = "Runtime (seconds)", color = "black", dpi = 150, linestyle=:dot)
# plot!(nn2, tt2, label = "Rastrigin Runtime", legend=:topleft, xlabel = "n", ylabel = "Runtime (seconds)", color = "black", dpi = 150, linestyle=:dash)

# savefig(plt, "png/BasinHoppingDimensionalRuntimn")

# xNorm = []
# diff = []
# data = CSV.File("csv/ConeSearchEasom.csv")

# for row in data
#     append!(xNorm, norm([row.x, row.y] - ones(2)*pi))
#     append!(diff, row.ab_diff)
# end

# plt = plot(xNorm, title = "Cone Search Convergenge - Easom Function", label = "||x - x*||", color = "black", xlabel = "Algorithm Step", yaxis=:log)
# plot!(diff, label = "β - α", color = "black", linestyle=:dash)
# savefig(plt, "png/EasomConvergence.png")

# data = CSV.File("csv/ConeSearchVariableLRuntimes.csv")
# xx = []
# yy = []
# for row in data
#     append!(yy, row.Runtime)
#     append!(xx, row.L)
# end
# plt = plot(xx, yy, xlabel = "L", ylabel = "Runtime (Seconds)", title = "Cone Search Runtime v. L", legend = false, dpi = 150, color = "Black")
# savefig(plt, "png/VariableConesearchRuntimes.png")

# function GetTimeData1(data, sums)
#     sum = 0
#     for row in data
#         sum += row.time
#     end
#     append!(sums, sum)
#     return sums
# end

# function GetTimeData2(data, sums)
#     sum = 0
#     for row in data
#         sum += row.step_time
#     end
#     append!(sums, sum)
#     return sums
# end

# sums_b = []
# sums_c = []

# gData_b = CSV.File("csv/basinHopping-Griewank.csv")
# sums_b = GetTimeData1(gData_b, sums_b)
# gData_c = CSV.File("csv/ConeSearchGriewank.csv")
# sums_c = GetTimeData2(gData_c, sums_c)

# aData_b = CSV.File("csv/basinHopping-Ackley.csv")
# sums_b = GetTimeData1(aData_b, sums_b)
# aData_c = CSV.File("csv/ConeSearchAckley.csv")
# sums_c = GetTimeData2(aData_c, sums_c)

# eData_b = CSV.File("csv/basinHopping-Easom.csv")
# sums_b = GetTimeData1(eData_b, sums_b)
# eData_c = CSV.File("csv/ConeSearchEasom.csv")
# sums_c = GetTimeData2(eData_c, sums_c)

# rData_b = CSV.File("csv/basinHopping-Rastrigin.csv")
# sums_b = GetTimeData1(rData_b, sums_b)
# rData_c = CSV.File("csv/ConeSearchRastrigin.csv")
# sums_c = GetTimeData2(rData_c, sums_c)

# dData_b = CSV.File("csv/basinHopping-Drop_Wave.csv")
# sums_b = GetTimeData1(dData_b, sums_b)
# dData_c = CSV.File("csv/ConeSearchDropWave.csv")
# sums_c = GetTimeData2(dData_c, sums_c)

# @show sums_c
# @show sums_b

# ticklabel = ["Griewank", "Ackley", "Easom", "Rastrigin", "Drop Wave"]
# plt = groupedbar([sums_c sums_b],legend=:outertop, xlabel = "Function", ylabel = "Time (Seconds)", color = ["Black" "Gray"], yaxis=:log, ylims=(1e-2, 1e3), xticks = (1:5, ticklabel), label = ["Cone Search" "Basin Hopping"], dpi = 150, title = "Runtime: Cone Search v. Basin Hopping")

# savefig(plt, "png/RuntimeBarComparision.png")