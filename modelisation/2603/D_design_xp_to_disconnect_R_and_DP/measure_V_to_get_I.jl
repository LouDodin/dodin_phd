using DifferentialEquations
using Plots
using CSV
using DataFrames
using Statistics

include(joinpath(@__DIR__, "models_list.jl"))

struct ModelSpec
    name::String
    dynamics!::Function
end

MODELS_TO_PLOT = ["SIVi"]

# S0 and I0 taken from t = 2.002002002002002j of coculture HV
X0 = log.([3.7166186443442135,11194.791057986662,1.000000000000001e-12,1.000000000000001e-12,1.000000000000001e-12,1.000000000000001e-12])

tspan = (0.0, 1/24)
t = range(0, 1/24, length=12)

const μ = 0.022416495821516925
const k = 9.78512085186654e7
const φi = 2.18e-7
const β = 51.7
const ω = 1.00e-8
const δ = 0.0150
const η = 3.01
const εdp = 0.664
const σdp = 0.257
const γ = 1.2e-5

# Colors
col_H    = :royalblue3
col_S    = :lightskyblue
col_I    = :royalblue3
col_R    = :navy
col_V    = :darkorchid3
col_Vi   = :plum2
col_Vdp  = :darkorchid3
col_Vdip = :indigo

# =====================
# Plot
# =====================
n_models = length(MODELS_TO_PLOT)
plt = plot(layout=(n_models,2), grid=true, yscale=:log10, size=(1500, 450*n_models))

for (idx, model_name) in enumerate(MODELS_TO_PLOT)
    println("Simulating model: $model_name")
    model_spec = MODELS[model_name]

    # Solve ODE
    prob = ODEProblem(model_spec.dynamics!, X0, tspan)
    sol = solve(prob, Tsit5(), saveat=t, abstol=1e-9, reltol=1e-6)

    # Recover real values
    S    = exp.(sol[1, :])
    I    = exp.(sol[2, :])
    R    = exp.(sol[3, :])
    Vi   = exp.(sol[4, :])
    Vdp  = exp.(sol[5, :])
    Vdip = exp.(sol[6, :])

    H = S .+ I .+ R
    V = Vi .+ Vdp .+ Vdip

    # Apply minimum floor
    floor = 1e-3
    H_plot   = max.(H, floor)
    V_plot   = max.(V, floor)
    S_plot   = max.(S, floor)
    I_plot   = max.(I, floor)
    R_plot   = max.(R, floor)
    Vi_plot  = max.(Vi, floor)
    Vdp_plot = max.(Vdp, floor)
    Vdip_plot= max.(Vdip, floor)

    left_plot    = 2*idx - 1
    right_plot = 2*idx

    # Total H and V
    plot!(plt[left_plot], t*24*60, H_plot, lw=3, label="H ($model_name)", color=col_H)
    plot!(plt[left_plot], t*24*60, V_plot, lw=3, label="V ($model_name)", color=col_V)
    ylabel!(plt[left_plot], "Abundance (log10)")
    xlabel!(plt[left_plot], "Time (min)")

    # Compartments
    plot!(plt[right_plot], t*24*60, S_plot,    lw=2, label="S", color=col_S)
    plot!(plt[right_plot], t*24*60, I_plot,    lw=2, label="I", color=col_I)
    plot!(plt[right_plot], t*24*60, R_plot,    lw=2, label="R", color=col_R)
    plot!(plt[right_plot], t*24*60, Vi_plot,   lw=2, label="Vi", color=col_Vi)
    plot!(plt[right_plot], t*24*60, Vdp_plot,  lw=2, label="Vdp", color=col_Vdp)
    plot!(plt[right_plot], t*24*60, Vdip_plot, lw=2, label="Vdip", color=col_Vdip)
    ylabel!(plt[right_plot], "Abundance (log10)")
    xlabel!(plt[right_plot], "Time (min)")


    # =====================
    # Save CSV per model
    # =====================
    df = DataFrame(
        time = t,
        H = H,
        S = S,
        I = I,
        R = R,
        V = V,
        Vi = Vi,
        Vdp = Vdp,
        Vdip = Vdip
    )

    csv_filename = joinpath(@__DIR__, "output_measure_V/$(model_name)_timeseries.csv")
    CSV.write(csv_filename, df)

    println("Saved CSV: $csv_filename")


    # Pente initiale

    t_small = t[1:2]
    Vi_small = Vi[1:2]
    println(t_small)

    slope = cov(t_small, Vi_small) / var(t_small)

    I0_est = slope / (β * η)
    println(slope)
    println(I0_est)
end

display(plt)



