using JuMP
using Ipopt

m4 = Model(Ipopt.Optimizer)

@variable(m4, -100 <= z1 <= 100)
@variable(m4, -100 <= z2 <= 100)
@variable(m4, -100 <= z3 <= 100)
@variable(m4, -100 <= z4 <= 100)
@variable(m4, -100 <= z5 <= 100)


@NLobjective(m4, Min, z2^2 - z3^2 + z4^2 + ℯ^(z1*z5^2*z1^2) + sin(z1*z2*z3*z4*z5 + 1)^2)

JuMP.optimize!(m4)

print(value(z1)," | ", value(z2), " | ",value(z3)," | ", value(z4), " | ",value(z5))
