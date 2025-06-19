using Plots
using Plots.PlotMeasures
using OptimalControl
using NLPModelsIpopt
nothing # hide

# Regularization of the PWL dynamics
function s⁺(x, θ, regMethod)
    if regMethod == 1 # Hill
        out = x^k/(x^k + θ^k)
    elseif regMethod == 2 # Exponential
        out = 1 - 1/(1 + exp(k*(x-θ)))
    end
    return out
end

# Regularization of |u(t) - 1|
function abs_m1(u, regMethod)
    if regMethod == 1 # Hill
        out = (u^k - 1)/(u^k + 1)
    elseif regMethod == 2 # Exponential
        out = 1 - 2/(1 + exp(k*(u-1)))
    end
    return out*(u - 1)
end
nothing # hide

# Constant definition
k₁    = 2;    k₂    = 3     # Production rates
γ₁    = 0.2;  γ₂    = 0.3   # Degradation rates
θ₁    = 4;    θ₂    = 3     # Transcriptional thresholds
uₘᵢₙ  = 0.6;  uₘₐₓ  = 1.4   # Control bounds
x₁ᶜ   = 4.7                 # Cycle point (initial and final)
λ     = 0.5                 # Trade-off fuel/time

# Initial guest for the NLP
tf    = 1
u(t)  = 1
sol = (control=u, variable=tf)

# Optimal control problem definition
ocp = @def begin

    tf ∈ R,                variable
    t ∈ [ 0, tf ],         time
    x = ( x₁, x₂ ) ∈ R²,   state
    u ∈ R,                 control

    x₁(0) == x₁ᶜ
    x₂(0) == θ₂
    x₁(tf) == x₁ᶜ
    x₂(tf) == θ₂

    uₘᵢₙ ≤ u(t) ≤ uₘₐₓ
    tf ≥ 1 # Force the state out of the confort zone

    ẋ(t) == [ - γ₁*x₁(t) + k₁*u(t)*(1 - s⁺(x₂(t),θ₂,regMethod))  ,
              - γ₂*x₂(t) + k₂*s⁺(x₁(t),θ₁,regMethod) ]

    ∫(λ*abs_m1(u(t),regMethod) + 1-λ) → min

end
nothing # hide

regMethod = 1       # Hill regularization
ki = 10             # Value of k for the first iteration
N = 400
maxki = 30          # Value of k for the last iteration
while ki < maxki
    global ki += 10  # Iteration step
    local print_level = (ki == maxki) # Only print the output on the last iteration
    global k = ki
    global sol = solve(ocp; grid_size=N, init=sol, print_level=4*print_level)
end
nothing # hide

plt1 = plot()
plt2 = plot()

tf    = variable(sol)
tspan = range(0, tf, N)   # time interval
x₁(t) = state(sol)(t)[1]
x₂(t) = state(sol)(t)[2]
u(t)  = control(sol)(t)

xticks = ([0, θ₁, x₁ᶜ], ["0", "θ₁", "x₁ᶜ"])
yticks = ([0, θ₂], ["0", "θ₂"])

plot!(plt1, x₁.(tspan), x₂.(tspan), label="optimal trajectory", xlabel="x₁", ylabel="x₂", xlimits=(θ₁/1.5, 1.1*x₁ᶜ), ylimits=(θ₂/2, 1.75*θ₂))
xticks!(xticks)
yticks!(yticks)
plot!(plt2, tspan, u, label="optimal control", xlabel="t")
plot(plt1, plt2; layout=(1,2), size=(800,300))

regMethod = 2       # Exponential regularization
ki = 100             # Value of k for the first iteration
N = 400
maxki = 400          # Value of k for the last iteration
while ki < maxki
    global ki += 100  # Iteration step
    local print_level = (ki == maxki) # Only print the output on the last iteration
    global k = ki
    global sol = solve(ocp; grid_size=N, init=sol, print_level=4*print_level)
end
nothing # hide

plt1 = plot()
plt2 = plot()

tf    = variable(sol)
tspan = range(0, tf, N)   # time interval
x₁(t) = state(sol)(t)[1]
x₂(t) = state(sol)(t)[2]
u(t)  = control(sol)(t)

xticks = ([0, θ₁, x₁ᶜ], ["0", "θ₁", "x₁ᶜ"])
yticks = ([0, θ₂], ["0", "θ₂"])

plot!(plt1, x₁.(tspan), x₂.(tspan), label="optimal trajectory", xlabel="x₁", ylabel="x₂", xlimits=(θ₁/1.5, 1.1*x₁ᶜ), ylimits=(θ₂/2, 1.75*θ₂))
xticks!(xticks)
yticks!(yticks)
plot!(plt2, tspan, u, label="optimal control", xlabel="t")
plot(plt1, plt2; layout=(1,2), size=(800,300))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
