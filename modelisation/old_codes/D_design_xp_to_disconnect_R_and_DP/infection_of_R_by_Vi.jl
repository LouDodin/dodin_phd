using DifferentialEquations
using Plots
using CSV
using DataFrames

MODELS_TO_PLOT = [
    "SIVi",
    # "SIRVi",
    # "SIViVdp",
]

include(joinpath(@__DIR__, "models_list.jl"))

best_df = CSV.read(joinpath(@__DIR__, "best.csv"), DataFrame)

# =========================================================
# --- Initial conditions (log-space, 7 variables now!)
# =========================================================
X0 = log.([1e-12, 1e-12, 1e7, 1e6, 1e-12, 1e-12, 1e-12])  
# S, I, R, Vi, Vdp, Vdip, Ev

tspan = (0.0, 240.0)
t = range(tspan[1], tspan[2], length = 1160)

# =========================================================
# --- Global constants (always present)
# =========================================================
const μ = 0.022416495821516925
const k = 9.78512085186654e7

println("Import des constantes fini")

# =========================================================
# --- Full ordered parameter list (IMPORTANT)
# =========================================================
const PARAM_ORDER = [
    :μ, :k, :φi, :β, :δ, :η,
    :εdp, :σdp,
    :μ_r, :k_r, :ν, :α
]

# =========================================================
# --- Plot setup
# =========================================================
n_models = length(MODELS_TO_PLOT)

plt = plot(
    layout = (2*n_models, 1),
    grid = true,
    yscale = :log10,
    size = (900, 450*n_models)
)

col_H = :royalblue3
col_S = :lightskyblue
col_I = :royalblue3
col_R = :navy
col_V = :darkorchid3
col_Vi = :plum2
col_Vdp = :darkorchid3
col_Vdip = :indigo

# =========================================================
# --- Loop over models
# =========================================================
for (idx, model_name) in enumerate(MODELS_TO_PLOT)

    println("Simulation du modèle : ", model_name)

    model_spec = MODELS[model_name]

    # =====================================================
    # --- Extract row from CSV
    # =====================================================
    row = best_df[best_df.model .== model_name, :]
    if nrow(row) == 0
        error("Modèle $model_name non trouvé dans best.csv")
    end

    # =====================================================
    # --- Build parameter dictionary (constants + fitted)
    # =====================================================
    θ_dict = Dict{Symbol, Float64}()

    # constantes
    θ_dict[:μ] = μ
    θ_dict[:k] = k

    for pname in propertynames(row)

        # 🔴 ignorer la colonne "model"
        pname == :model && continue

        val = row[1, pname]

        if ismissing(val)
            continue
        elseif val isa Number
            θ_dict[pname] = Float64(val)
        elseif val isa AbstractString
            if isempty(strip(val))
                continue
            end
            try
                θ_dict[pname] = parse(Float64, val)
            catch
                # ignore si non numérique
            end
        end
    end

    # =====================================================
    # --- Convert to ordered vector p
    # =====================================================
    p = Float64[]
    for pname in PARAM_ORDER
        push!(p, get(θ_dict, pname, 0.0))
    end

    # =====================================================
    # --- Solve ODE
    # =====================================================
    prob = ODEProblem(model_spec.dynamics!, X0, tspan, p)
    sol = solve(prob, Tsit5(), saveat = t)

    # =====================================================
    # --- Back to real space (IMPORTANT FIX)
    # =====================================================
    S    = exp.(sol[1, :])
    I    = exp.(sol[2, :])
    R    = exp.(sol[3, :])
    Vi   = exp.(sol[4, :])
    Vdp  = exp.(sol[5, :])
    Vdip = exp.(sol[6, :])

    # sécurité numérique
    S    = max.(S, 1e-12)
    I    = max.(I, 1e-12)
    R    = max.(R, 1e-12)
    Vi   = max.(Vi, 1e-12)
    Vdp  = max.(Vdp, 1e-12)
    Vdip = max.(Vdip, 1e-12)

    H = S .+ I .+ R
    V = Vi .+ Vdp .+ Vdip

    # =====================================================
    # --- Subplot index
    # =====================================================
    top_plot    = 2*idx - 1
    bottom_plot = 2*idx

    # =====================================================
    # --- Plot H & V
    # =====================================================
    plot!(plt[top_plot], t./24, H, lw=3, label="H ($model_name)", color=col_H)
    plot!(plt[top_plot], t./24, V, lw=3, label="V ($model_name)", color=col_V)

    ylabel!(plt[top_plot], "Abundance (log10)")
    xlabel!(plt[top_plot], "Time (day)")

    # =====================================================
    # --- Compartments
    # =====================================================
    plot!(plt[bottom_plot], t./24, S,    lw=2, label="S", color=col_S)
    plot!(plt[bottom_plot], t./24, I,    lw=2, label="I", color=col_I)
    plot!(plt[bottom_plot], t./24, R,    lw=2, label="R", color=col_R)
    plot!(plt[bottom_plot], t./24, Vi,   lw=2, label="Vi", color=col_Vi)
    plot!(plt[bottom_plot], t./24, Vdp,  lw=2, label="Vdp", color=col_Vdp)
    plot!(plt[bottom_plot], t./24, Vdip, lw=2, label="Vdip", color=col_Vdip)

    ylabel!(plt[bottom_plot], "Abundance (log10)")
    xlabel!(plt[bottom_plot], "Time (day)")
end

display(plt)