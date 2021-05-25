#see https://en.wikipedia.org/wiki/Test_functions_for_optimization


function Rastrigin(x,n)
    sum = 0

    for i = 1:n
        sum += x[i]^2 - 10cos(2pi*x[i])
    end

    return 10n + sum
end


function Ackley(x)
    return -20exp(-0.2*sqrt(0.5*(x[1]^2 + x[2]^2))) - exp(0.5*(cos(2*pi*x[1]) + cos(2*pi*x[2]))) + exp(1) + 20
end
