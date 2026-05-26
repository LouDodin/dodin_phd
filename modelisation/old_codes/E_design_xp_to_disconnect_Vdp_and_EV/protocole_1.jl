using DifferentialEquations
using Plots
using Measures
using LaTeXStrings

include(joinpath(@__DIR__, "models_list.jl"))

# =====================
# Gradient colors
# =====================

virus_palette = cgrad([:darkred, :orange])
col_Vi   = virus_palette[0.2]

col_partic = :purple

# =====================
# Models
# =====================

model_ctrl   = MODELS["SIViVdpEv"]
model_ev     = MODELS["SIViEv"]
model_vdp    = MODELS["SIViVdp"]
model_mix    = MODELS["SIViVdpEv"]

model_list = [model_ctrl, model_ev, model_vdp, model_mix]

# =====================
# Initial conditions (log-space)
# =====================

X0_ctrl = log.([1e-12, 1e-12, 1e-12, 1e10, 1e-12, 1e-12, 1e-12])
X0_ev   = log.([1e-12, 1e-12, 1e-12, 1e10, 1e-12, 1e-12, 1e7])
X0_vdp  = log.([1e-12, 1e-12, 1e-12, 1e10, 1e7, 1e-12, 1e-12])
X0_mix  = log.([1e-12, 1e-12, 1e-12, 1e10, 1e7, 1e-12, 1e7])

X0_list = [X0_ctrl, X0_ev, X0_vdp, X0_mix]

title_list = [
    "Control : Vi only (Model : SIViVdpEv)",
    "Particles : 100% EV (Model : SIViEv)",
    "Particles : 100 % Vdp (Model : SIViVdp)",
    "Particles : 50% Vdp + 50% EV (Model : SIViVdpEv)"
]

# =====================
# Time
# =====================

tspan = (0.0, 10.0)
t = range(0, 10, length=100)

# =====================
# Parameters
# =====================

const μ = 0.7
const k = 1e9
const φi = 2.2e-7
const β = 52
const δ = 1.5e-2
const η = 3
const εdp = 0.66
const σdp = 0.26
const μ_r = 0.7
const k_r = 1e9
const ν = 1

θ = (μ, k, φi, β, δ, η, εdp, σdp, μ_r, k_r, ν)

# =====================
# Layout 2 × 2
# =====================

plt = plot(layout=(2, 2),
           grid=true,
           yscale=:log10,
           ylim=(1e5, 1e11),
           yticks=([10.0^i for i in 5:1:10],
                   [L"10^{%$i}" for i in 5:1:10]),
           size=(3200, 2400),
           margins=20mm,
           legendfontsize=18,
           guidefontsize=18,
           tickfontsize=18,
           titlefontsize=24,
           plot_title="Dynamics after infection depending on particle type",
           plot_titlefontsize=28)

floor_val = 1e5
mask_floor(x) = ifelse.(x .< floor_val, NaN, x)

# =====================
# Simulation
# =====================

for i in 1:4

    model_spec = model_list[i]
    X0 = X0_list[i]

    prob = ODEProblem(model_spec.dynamics!, X0, tspan, θ)
    sol = solve(prob, Tsit5(), saveat=t)

    nvar = size(sol, 1)

    # Récupération des variables de base
    S  = nvar ≥ 1 ? exp.(sol[1, :]) : zeros(length(t))
    I  = nvar ≥ 2 ? exp.(sol[2, :]) : zeros(length(t))
    R  = nvar ≥ 3 ? exp.(sol[3, :]) : zeros(length(t))
    Vi = nvar ≥ 4 ? exp.(sol[4, :]) : zeros(length(t))

    # Variables supplémentaires selon le modèle
    Vdp  = nvar ≥ 5 ? exp.(sol[5, :]) : zeros(length(t))
    Vdip = nvar ≥ 6 ? exp.(sol[6, :]) : zeros(length(t))
    Ev   = nvar ≥ 7 ? exp.(sol[7, :]) : zeros(length(t))

    # Calcul des totaux
    partic = Vdp .+ Vdip .+ Ev

    # Plot
    plot!(plt[i], t, mask_floor(Vi),
          lw=4, color=col_Vi, label="Vi")

    if i != 1
        plot!(plt[i], t, mask_floor(partic),
              lw=4, color=col_partic, label="Particles")
    end

    title!(plt[i], title_list[i])
    xlabel!(plt[i], "Time (day)")
    ylabel!(plt[i], "Concentration (parts/mL)")
    plot!(plt[i], legend=:bottomright)
end

savefig(plt, joinpath(@__DIR__, "output_protocole_1/dynamics.png"))
display(plt)