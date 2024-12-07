var documenterSearchIndex = {"docs":
[{"location":"bistable.html#Bistable-toggle-switch","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"","category":"section"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"The simplest GRN exhibiting bistability can be modeled through two variables x_1 and x_2 that mutually repress each othe. We suppose that the system can be externally controlled by a chemical inducer that targets the synthesis rates of both genes. The model is defined as","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"left beginarrayl\ndotx_1 = -gamma_1 x_1 + u(t) k_1 s^-(x_2theta_2)  \ndotx_2 = -gamma_2 x_2 + u(t) k_2  s^-(x_1theta_1)\nendarrayright","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"where a detail of each term can be found in Introduction. The domain K=0k_1  gamma_1times 0k_2  gamma_2 is forward invariant by the dynamics, which divides the state space into four regions:","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"beginarrayl\nB_00=left(x_1x_2)in R^2 mid 0x_1theta_1  0x_2theta_2right\nB_01=left(x_1x_2)in R^2 mid 0x_1theta_1  theta_2x_2frack_2gamma_2right\nB_10=left(x_1x_2)in R^2 mid theta_1x_1frack_1gamma_1  0x_2theta_2right\nB_11=left(x_1x_2)in R^2 mid theta_1x_1frack_1gamma_1  theta_2x_2frack_2gamma_2right\nendarray","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"and two locally asymptotically stable steady states","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"beginalign\n    nonumber phi_10 = left(frack_1gamma_10right)in barB_10 \n    nonumber phi_01 = left(0frack_2gamma_2right)in barB_01\nendalign","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"as shown in the following figure:","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"(Image: Alt Text)","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"The control objective is to induce a transition between an initial point in B_10 and a final value of x_2 in B_01, which can be written through the initial and terminal constraints:","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"    x(0) = x_0 in B_10 qquad x_1(t_f)  theta_1 qquad x_2(t_f) = x_2^f ","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"for free final time t_f  0 and for x_2^f  theta_2.","category":"page"},{"location":"bistable.html#Problem-definition","page":"Bistable toggle switch","title":"Problem definition","text":"","category":"section"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"using Plots\nusing Plots.PlotMeasures\nusing OptimalControl\nusing NLPModelsIpopt\nnothing # hide","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"We define the regularization functions, where the method is decided through the argument regMethod.","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"# Regularization of the PWL dynamics\nfunction s⁺(x, θ, regMethod)\n    if regMethod == 1 # Hill\n        out = x^k/(x^k + θ^k)\n    elseif regMethod == 2 # Exponential\n        out = 1 - 1/(1 + exp(k*(x-θ)))\n    end\n    return out\nend\n\n# Regularization of |u(t) - 1|\nfunction abs_m1(u, regMethod)\n    if regMethod == 1 # Hill\n        out = (u^k - 1)/(u^k + 1)\n    elseif regMethod == 2 # Exponential\n        out = 1 - 2/(1 + exp(k*(u-1)))\n    end\n    return out*(u - 1)\nend\nnothing # hide","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"Definition of the OCP:","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"# Constant definition\nk₁    = 1;    k₂    = 1     # Production rates\nγ₁    = 1.4;  γ₂    = 1.6   # Degradation rates\nθ₁    = 0.6;  θ₂    = 0.4   # Transcriptional thresholds\nuₘᵢₙ  = 0.6;  uₘₐₓ  = 1.4   # Control bounds\nx₀    = [0.65, 0.2]         # Initial point\nx₂ᶠ   = 0.55                # Final point\nλ     = 0.25                # Trade-off fuel/time\n\n# Initial guest for the NLP\ntf    = 1.5\nu(t)  = 0\nsol = (control=u, variable=tf)\n\n# Optimal control problem definition\nocp = @def begin\n\n    tf ∈ R,                variable\n    t ∈ [ 0, tf ],         time\n    x = [ x₁, x₂ ] ∈ R²,   state\n    u ∈ R,                 control \n\n    x(0) == x₀\n    x₁(tf) ≤ θ₁\n    x₂(tf) == x₂ᶠ\n\n    uₘᵢₙ ≤ u(t) ≤ uₘₐₓ\n    tf ≥ 0\n\n    ẋ(t) == [ - γ₁*x₁(t) + k₁*u(t)*(1 - s⁺(x₂(t),θ₂,regMethod))  ,\n              - γ₂*x₂(t) + k₂*u(t)*(1 - s⁺(x₁(t),θ₁,regMethod)) ]\n\n    ∫(λ*abs_m1(u(t),regMethod) + 1-λ) → min      \n\nend\nnothing # hide","category":"page"},{"location":"bistable.html#Resolution-through-Hill-regularization","page":"Bistable toggle switch","title":"Resolution through Hill regularization","text":"","category":"section"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"In order to ensure convergence of the solver, we solve the OCP by iteratively increasing the parameter k while using the i-1-th solution as the initialization of the i-th iteration.","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"regMethod = 1       # Hill regularization\nki = 50             # Value of k for the first iteration\nN = 400\nmaxki = 200          # Value of k for the last iteration\nwhile ki < maxki\n    global ki += 50  # Iteration step\n    local print_level = (ki == maxki) # Only print the output on the last iteration\n    global k = ki\n    global sol = solve(ocp; grid_size=N, init=sol, print_level=4*print_level)\nend\nnothing # hide","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"Plotting of the results:","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"plt1 = plot()\nplt2 = plot()\n\ntf    = sol.variable\ntspan = range(0, tf, N)   # time interval\nx₁(t) = sol.state(t)[1]\nx₂(t) = sol.state(t)[2]\nu(t)  = sol.control(t)\n\nxticks = ([0, θ₁], [\"0\", \"θ₁\"])\nyticks = ([0, θ₂, x₂ᶠ], [\"0\", \"θ₂\", \"x₂ᶠ\"])\n\nplot!(plt1, x₁.(tspan), x₂.(tspan), label=\"optimal trajectory\", xlabel=\"x₁\", ylabel=\"x₂\", xlimits=(θ₁/3, k₁/γ₁), ylimits=(0, k₂/γ₂))\nscatter!(plt1, [x₀[1]], [x₀[2]], label=\"x₀\", color=:deepskyblue)\nxticks!(xticks)\nyticks!(yticks)\nplot!(plt2, tspan, u, label=\"optimal control\", xlabel=\"t\")\nplot(plt1, plt2; layout=(1,2), size=(800,300))","category":"page"},{"location":"bistable.html#Resolution-through-exponential-regularization","page":"Bistable toggle switch","title":"Resolution through exponential regularization","text":"","category":"section"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"The same procedure for iteratively increasing k is used.","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"regMethod = 2       # Exponential regularization\nki = 50             # Value of k for the first iteration\nN = 400\nmaxki = 300          # Value of k for the last iteration\nwhile ki < maxki\n    global ki += 50  # Iteration step\n    local print_level = (ki == maxki) # Only print the output on the last iteration\n    global k = ki\n    global sol = solve(ocp; grid_size=N, init=sol, print_level=4*print_level)\nend\nnothing # hide","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"Plotting of the results:","category":"page"},{"location":"bistable.html","page":"Bistable toggle switch","title":"Bistable toggle switch","text":"plt1 = plot()\nplt2 = plot()\n\ntf    = sol.variable\ntspan = range(0, tf, N)   # time interval\nx₁(t) = sol.state(t)[1]\nx₂(t) = sol.state(t)[2]\nu(t)  = sol.control(t)\n\nxticks = ([0, θ₁], [\"0\", \"θ₁\"])\nyticks = ([0, θ₂, x₂ᶠ], [\"0\", \"θ₂\", \"x₂ᶠ\"])\n\nplot!(plt1, x₁.(tspan), x₂.(tspan), label=\"optimal trajectory\", xlabel=\"x₁\", ylabel=\"x₂\", xlimits=(θ₁/3, k₁/γ₁), ylimits=(0, k₂/γ₂))\nscatter!(plt1, [x₀[1]], [x₀[2]], label=\"x₀\", color=:deepskyblue)\nxticks!(xticks)\nyticks!(yticks)\nplot!(plt2, tspan, u, label=\"optimal control\", xlabel=\"t\")\nplot(plt1, plt2; layout=(1,2), size=(800,300))","category":"page"},{"location":"oscillator.html#Damped-genetic-oscillator","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"","category":"section"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"The simplest GRN exhibiting oscillatory behaviors can be modeled through two variables x_1 and x_2 with opposite mutual effects: x_1 catalyzes the production of x_2, that in turn inhibits the production of x_1. We suppose that the system can be externally controlled by a chemical inducer that targets only one of the genes. The model is defined as","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"left beginarrayl\ndotx_1 = -gamma_1 x_1 + u(t) k_1 s^-(x_2theta_2)  \ndotx_2 = -gamma_2 x_2 + k_2  s^+(x_1theta_1)\nendarrayright","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"where a detail of each term can be found in Introduction.","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"For any initial condition, the openloop system (i.e, u equiv 1) converges to the equilibrium point (theta_1 theta_2) when t rightarrow infty, producing a damped oscillatory behavior[1]:","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"(Image: Alt Text)","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"The control objective is to induce a sustained oscillation. Thus, we can state the problem of producing a single cycle, which can be written through the initial and terminal constraints:","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"    x(0) = x(t_f) = (x_1^c theta_2 ) ","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"for free final time t_f  0 and for x_1^c  theta_1.","category":"page"},{"location":"oscillator.html#Problem-definition","page":"Damped genetic oscillator","title":"Problem definition","text":"","category":"section"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"using Plots\nusing Plots.PlotMeasures\nusing OptimalControl\nusing NLPModelsIpopt\nnothing # hide","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"We define the regularization functions, where the method is decided through the argument regMethod.","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"# Regularization of the PWL dynamics\nfunction s⁺(x, θ, regMethod)\n    if regMethod == 1 # Hill\n        out = x^k/(x^k + θ^k)\n    elseif regMethod == 2 # Exponential\n        out = 1 - 1/(1 + exp(k*(x-θ)))\n    end\n    return out\nend\n\n# Regularization of |u(t) - 1|\nfunction abs_m1(u, regMethod)\n    if regMethod == 1 # Hill\n        out = (u^k - 1)/(u^k + 1)\n    elseif regMethod == 2 # Exponential\n        out = 1 - 2/(1 + exp(k*(u-1)))\n    end\n    return out*(u - 1)\nend\nnothing # hide","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"Definition of the OCP:","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"# Constant definition\nk₁    = 2;    k₂    = 3     # Production rates\nγ₁    = 0.2;  γ₂    = 0.3   # Degradation rates\nθ₁    = 4;    θ₂    = 3     # Transcriptional thresholds\nuₘᵢₙ  = 0.6;  uₘₐₓ  = 1.4   # Control bounds\nx₁ᶜ   = 4.7                 # Cycle point (initial and final)\nλ     = 0.5                 # Trade-off fuel/time\n\n# Initial guest for the NLP\ntf    = 1\nu(t)  = 2\nsol = (control=u, variable=tf)\n\n# Optimal control problem definition\nocp = @def begin\n\n    tf ∈ R,                variable\n    t ∈ [ 0, tf ],         time\n    x = [ x₁, x₂ ] ∈ R²,   state\n    u ∈ R,                 control \n\n    x₁(0) == x₁ᶜ\n    x₂(0) == θ₂\n    x₁(tf) == x₁ᶜ\n    x₂(tf) == θ₂\n\n    uₘᵢₙ ≤ u(t) ≤ uₘₐₓ\n    tf ≥ 1 # Force the state out of the confort zone\n\n    ẋ(t) == [ - γ₁*x₁(t) + k₁*u(t)*s⁺(x₂(t),θ₂,regMethod)  ,\n              - γ₂*x₂(t) + k₂*(1 - s⁺(x₁(t),θ₁,regMethod)) ]\n\n    ∫(λ*abs_m1(u(t),regMethod) + 1-λ) → min      \n\nend\nnothing # hide","category":"page"},{"location":"oscillator.html#Resolution-through-Hill-regularization","page":"Damped genetic oscillator","title":"Resolution through Hill regularization","text":"","category":"section"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"In order to ensure convergence of the solver, we solve the OCP by iteratively increasing the parameter k while using the i-1-th solution as the initialization of the i-th iteration.","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"regMethod = 1       # Hill regularization\nki = 10             # Value of k for the first iteration\nN = 400\nmaxki = 30          # Value of k for the last iteration\nwhile ki < maxki\n    global ki += 10  # Iteration step\n    local print_level = (ki == maxki) # Only print the output on the last iteration\n    global k = ki\n    global sol = solve(ocp; grid_size=N, init=sol, print_level=4*print_level)\nend\nnothing # hide","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"Plotting of the results:","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"plt1 = plot()\nplt2 = plot()\n\ntf    = sol.variable\ntspan = range(0, tf, N)   # time interval\nx₁(t) = sol.state(t)[1]\nx₂(t) = sol.state(t)[2]\nu(t)  = sol.control(t)\n\nxticks = ([0, θ₁, x₁ᶜ], [\"0\", \"θ₁\", \"x₁ᶜ\"])\nyticks = ([0, θ₂], [\"0\", \"θ₂\"])\n\nplot!(plt1, x₁.(tspan), x₂.(tspan), label=\"optimal trajectory\", xlabel=\"x₁\", ylabel=\"x₂\", xlimits=(θ₁/1.25, 1.1*x₁ᶜ), ylimits=(θ₂/2, 1.5*θ₂))\nxticks!(xticks)\nyticks!(yticks)\nplot!(plt2, tspan, u, label=\"optimal control\", xlabel=\"t\")\nplot(plt1, plt2; layout=(1,2), size=(800,300))","category":"page"},{"location":"oscillator.html#Resolution-through-exponential-regularization","page":"Damped genetic oscillator","title":"Resolution through exponential regularization","text":"","category":"section"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"The same procedure for iteratively increasing k is used.","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"regMethod = 2       # Exponential regularization\nki = 20             # Value of k for the first iteration\nN = 400\nmaxki = 400          # Value of k for the last iteration\nwhile ki < maxki\n    global ki += 20  # Iteration step\n    local print_level = (ki == maxki) # Only print the output on the last iteration\n    global k = ki\n    global sol = solve(ocp; grid_size=N, init=sol, print_level=4*print_level)\nend\nnothing # hide","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"Plotting of the results:","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"plt1 = plot()\nplt2 = plot()\n\ntf    = sol.variable\ntspan = range(0, tf, N)   # time interval\nx₁(t) = sol.state(t)[1]\nx₂(t) = sol.state(t)[2]\nu(t)  = sol.control(t)\n\nxticks = ([0, θ₁, x₁ᶜ], [\"0\", \"θ₁\", \"x₁ᶜ\"])\nyticks = ([0, θ₂], [\"0\", \"θ₂\"])\n\nplot!(plt1, x₁.(tspan), x₂.(tspan), label=\"optimal trajectory\", xlabel=\"x₁\", ylabel=\"x₂\", xlimits=(θ₁/1.25, 1.1*x₁ᶜ), ylimits=(θ₂/2, 1.5*θ₂))\nxticks!(xticks)\nyticks!(yticks)\nplot!(plt2, tspan, u, label=\"optimal control\", xlabel=\"t\")\nplot(plt1, plt2; layout=(1,2), size=(800,300))","category":"page"},{"location":"oscillator.html","page":"Damped genetic oscillator","title":"Damped genetic oscillator","text":"[1]: E. Farcot, J.-L. Gouzé, Periodic solutions of piecewise affine gene network models with non uniform decay rates: the case of a negative feedback loop, Acta biotheoretica 57 (4) (2009) 429–455.","category":"page"},{"location":"index.html#PWL-models-of-gene-regulatory-networks","page":"Introduction","title":"PWL models of gene regulatory networks","text":"","category":"section"},{"location":"index.html#Context","page":"Introduction","title":"Context","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"This application case studies the resolution of OCPs through a direct numerical method in piecewise linear models of gene regulatory networks[1]. More specifically, we study state transitions between two points of the state space minimizing the nonsmooth cost function","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"J_lambda(u) = lambda int_0^t_f u(t)-1  dt + (1-lambda) t_f","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"where u(t) in u_min u_max is the external control, lambda in (01) a fixed parameter, and t_f is the free final time. In the general case, the concentration of the i-th gene is described by a dynamical equation of the form","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"dotx_i = -gamma_i x_i + k_i s^pm(x_jtheta_j)","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"where the positive constants gamma_i, k_i correspond, respectively, to the degradation and the production rates of x_i, and the gene expression rate s^pm is a piecewise constant function defined as:","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"        s^+(x theta) = left beginarrayll\n        0 quad textitif  x  theta \n        1 quad textitif  x  theta\n        endarray right quad\n        s^-(x theta) = 1 - s^+(x theta) = left beginarrayll\n        1 quad textitif  x  theta \n        0 quad textitif  x  theta\n        endarray right","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"where s^- models an inhibiting effect, s^+ a catalyzing effect, and theta represents in both cases a threshold for transcriptional repression or activation, respectively.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"In order to solve this numerical OCP, there are two difficulties:","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"The hybrid nature of the dynamics.\nThe non-smoothness of the L^1 Lagrangian cost.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"The latter can be tackled through a regularization scheme to be detailed in next section.","category":"page"},{"location":"index.html#Regularization-strategies","page":"Introduction","title":"Regularization strategies","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Two regularization strategies are compared: through Hill functions, and through exponential functions, both depending on a parameter k in mathbbN. Each case is detailed in the following table.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Function Hill Exponential\ns^+(x theta) fracx^kx^k + theta^k 1 - frac11 + e^k(x-theta)\nmid u-1 mid (u-1) fracu^k - 1u^k + 1 (u-1) left 1 - frac21 + e^k(u-1) right","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"A comparison for low values of k, for s^+(x theta):","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"using Plots\nθ = 2\nk = 10\nx = range(0, 4, length=100)\ny1 = (x .> θ)\ny2 = x.^k ./ (x.^k .+ θ^k)\ny3 = 1 .- 1 ./ (1 .+ exp.(k*(x.-θ)))\nplot(x, [y1, y2, y3], label=[\"s⁺\" \"Hill\" \"Exponential\"], xlabel=\"x\")","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"and for u-1:","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"using Plots\nk = 5\nu = range(0, 2, length=100)\ny1 = abs.(u .- 1)\ny2 = (u .- 1).*(u.^k .- 1)./(u.^k .+ 1)\ny3 = (u .- 1).*(1 .- 2 ./ (1 .+ exp.(k .* (u .- 1))))\nplot(u, [y1, y2, y3], label=[\"|u-1|\" \"Hill\" \"Exponential\"], xlabel=\"u\")","category":"page"},{"location":"index.html#Dependencies","page":"Introduction","title":"Dependencies","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"All the numerical simulations to generate this documentation are performed with the following packages.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"using Pkg\nPkg.status()","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"[1]: Agustín G. Yabo, Nicolas Augier. On L¹ and time-optimal state transitions in piecewise linear models of gene-regulatory networks. Preprint. 2024. https://hal.science/hal-04820387.","category":"page"}]
}