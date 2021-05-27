include("basin_hopping.jl")

max_dim = 150
it_per_dim = 5
acc = zeros(max_dim)


for i = 2:max_dim

    correct_count = 0.0
    f(x) = Rastrigin(x, i)

    for j = 1:it_per_dim
        x0 = rand(i) * 10.0
        minimum = Basin_Hopping(f, x0, α, β, η, ϵ, κ, ℓ, ℓ_range, γ, ϕ, T,
                                    target_acc_rate, maxIterations, static_threshold)
        if norm(minimum[1] - zeros(i)) <= 1e-6 && abs(minimum[2]) <= 1e-6
            correct_count += 1
        end
    end

    acc[i] = 100*correct_count/it_per_dim
end

plot(acc)
