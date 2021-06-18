using Plots
using Calculus

function Compute_Gradient(_f, _x)
    return Calculus.gradient(g -> _f(g),_x)
end

function One_Vec(n, pos)
   vec = zeros(n)
   vec[pos] = 1
   return vec
end

function Complex_Differentiation(_f, _x, _h)
   temp = zeros(length(_x))

   for i = 1:length(_x)
      temp[i] = imag(f(x + im*_h*One_Vec(length(_x), i)))/_h
   end

   return temp
end

f(x) = exp(x[1])/(cos(x[2])^3 + sin(x[3])^3)
x = [pi/4, 0, 10]
h = 1
ϵ = 1e-8

while h > ϵ
   println("Complex: ", Complex_Differentiation(f, x, h))
   h*=0.5
end


println("Finite diff: ", Compute_Gradient(f, x))
