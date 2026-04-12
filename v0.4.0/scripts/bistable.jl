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
k₁ = 1;
k₂ = 1     # Production rates
γ₁ = 1.4;
γ₂ = 1.6   # Degradation rates
θ₁ = 0.6;
θ₂ = 0.4   # Transcriptional thresholds
uₘᵢₙ = 0.6;
uₘₐₓ = 1.4   # Control bounds
x₀ = [0.65, 0.2]         # Initial point
x₂ᶠ = 0.55                # Final point
λ = 0.25                # Trade-off fuel/time

# Initial guest for the NLP
tf = 1.5
u(t) = 0
sol = (control=u, variable=tf)

# Optimal control problem definition
ocp = @def begin
    tf ∈ R, variable
    t ∈ [0, tf], time
    x = (x₁, x₂) ∈ R², state
    u ∈ R, control

    x(0) == x₀
    x₁(tf) ≤ θ₁
    x₂(tf) == x₂ᶠ

    uₘᵢₙ ≤ u(t) ≤ uₘₐₓ
    tf ≥ 0

    ẋ(t) == [
        - γ₁*x₁(t) + k₁*u(t)*(1 - s⁺(x₂(t), θ₂, regMethod)),
        - γ₂*x₂(t) + k₂*u(t)*(1 - s⁺(x₁(t), θ₁, regMethod)),
    ]

    ∫(λ*abs_m1(u(t), regMethod) + 1-λ) → min
end
nothing # hide

regMethod = 1       # Hill regularization
ki = 50             # Value of k for the first iteration
N = 400
maxki = 200          # Value of k for the last iteration
while ki < maxki
    global ki += 50  # Iteration step
    local print_level = (ki == maxki) # Only print the output on the last iteration
    global k = ki
    global sol = solve(ocp; grid_size=N, init=sol, print_level=4*print_level)
end
nothing # hide

plt1 = plot()
plt2 = plot()

tf = variable(sol)
tspan = range(0, tf, N)   # time interval
x₁(t) = state(sol)(t)[1]
x₂(t) = state(sol)(t)[2]
u(t) = control(sol)(t)

xticks = ([0, θ₁], ["0", "θ₁"])
yticks = ([0, θ₂, x₂ᶠ], ["0", "θ₂", "x₂ᶠ"])

plot!(
    plt1,
    x₁.(tspan),
    x₂.(tspan);
    label="optimal trajectory",
    xlabel="x₁",
    ylabel="x₂",
    xlimits=(θ₁/3, k₁/γ₁),
    ylimits=(0, k₂/γ₂),
)
scatter!(plt1, [x₀[1]], [x₀[2]]; label="x₀", color=:deepskyblue)
xticks!(xticks)
yticks!(yticks)
plot!(plt2, tspan, u; label="optimal control", xlabel="t")
plot(plt1, plt2; layout=(1, 2), size=(800, 300))

regMethod = 2       # Exponential regularization
ki = 50             # Value of k for the first iteration
N = 400
maxki = 300          # Value of k for the last iteration
while ki < maxki
    global ki += 50  # Iteration step
    local print_level = (ki == maxki) # Only print the output on the last iteration
    global k = ki
    global sol = solve(ocp; grid_size=N, init=sol, print_level=4*print_level)
end
nothing # hide

plt1 = plot()
plt2 = plot()

tf = variable(sol)
tspan = range(0, tf, N)   # time interval
x₁(t) = state(sol)(t)[1]
x₂(t) = state(sol)(t)[2]
u(t) = control(sol)(t)

xticks = ([0, θ₁], ["0", "θ₁"])
yticks = ([0, θ₂, x₂ᶠ], ["0", "θ₂", "x₂ᶠ"])

plot!(
    plt1,
    x₁.(tspan),
    x₂.(tspan);
    label="optimal trajectory",
    xlabel="x₁",
    ylabel="x₂",
    xlimits=(θ₁/3, k₁/γ₁),
    ylimits=(0, k₂/γ₂),
)
scatter!(plt1, [x₀[1]], [x₀[2]]; label="x₀", color=:deepskyblue)
xticks!(xticks)
yticks!(yticks)
plot!(plt2, tspan, u; label="optimal control", xlabel="t")
plot(plt1, plt2; layout=(1, 2), size=(800, 300))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
