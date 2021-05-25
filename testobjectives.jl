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


function Rosenbrock(x, n)
    sum = 0

    for i = 1:n-1
        sum += 100*(x[i+1] - x[i]^2)^2 + (1 - x[1])^2
    end

    return sum
end


function Beale(x)
    return (1.5 - x[1] + x[1]*x[2])^2 + (2.25 - x[1] + x[1]*x[2]^2)^2 + (2.625 - x[1] + x[1]*x[2]^3)^2
end


#This one is really hard to solve
function Bukin(x)
    return 100*sqrt(abs(x[2] - 0.01x[1]^2)) + 0.01*abs(x[1] + 10.0)
end


#Not differentiable
function Holder_Table(x)
    return -abs( sin(x[1])*cos(x[2]) * exp(abs(1 - sqrt(x[1]^2 + x[2]^2)/pi)))
end


function Schaffer_N2(x)
    0.5 + ((sin(x[1]^2 + x[2]^2))^2 - 0.5)/(1 + 0.001*(x[1]^2 + x[2]^2))^2
end


function Styblinski_Tang(x, n)
    sum = 0

    for i = 1:n
        sum += x[i]^4 - 16x[i]^2 + 5x[i]
    end

    return sum/2.0
end
