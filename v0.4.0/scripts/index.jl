using Plots
θ = 2
k = 10
x = range(0, 4; length=100)
y1 = (x .> θ)
y2 = x .^ k ./ (x .^ k .+ θ^k)
y3 = 1 .- 1 ./ (1 .+ exp.(k*(x .- θ)))
plot(x, [y1, y2, y3]; label=["s⁺" "Hill" "Exponential"], xlabel="x")

k = 5
u = range(0, 2; length=100)
y1 = abs.(u .- 1)
y2 = (u .- 1) .* (u .^ k .- 1) ./ (u .^ k .+ 1)
y3 = (u .- 1) .* (1 .- 2 ./ (1 .+ exp.(k .* (u .- 1))))
plot(u, [y1, y2, y3]; label=["|u-1|" "Hill" "Exponential"], xlabel="u")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
