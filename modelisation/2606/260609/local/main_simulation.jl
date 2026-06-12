## ===== Simulation SRVi / SViVd sur données expérimentales =====
#
# Simule cycle par cycle en propageant les proportions S/R et Vi/Vd
# comme dans le script d'optimisation d'origine.
# Les données expérimentales sont superposées aux simulations.
#
# Choisissez le modèle via MODEL_SELECT (~ligne 20)
# Usage : julia simulate.jl

using DifferentialEquations
using OrdinaryDiffEqRosenbrock
using CSV
using DataFrames
using Statistics
using Plots
using Measures
using LaTeXStrings
using SciMLBase


## ===== Choix du modèle =====
#   "SRVi"  — dS, dR, dVi
#   "SViVd" — dS, dVi, dVd

const MODEL_SELECT = "SRVi"


## ===== Paramètres =====

# --- SRVi ---
const SRVi_params = (
    r = 0.574619342477644,
    K = 6.675449070379925e7,
    β = 144.0,
    δ = 0.02,
    α = 7e-5,    # taux de mutation S→R (j⁻¹)
    φ = 6e-9,   # taux d'interaction (mL/(cell·j))
)

# --- SViVd ---
const SViVd_params = (
    r = 0.5592225270686286,
    K = 7.29695252684594e7,
    β = 144.0,
    δ = 0.02,
    φ = 1e-10,   # taux d'interaction (mL/(cell·j))
    ε = 0.5,     # fraction de virions infectieux
)


## ===== Choix des données =====
const replicates = ["A", "B", "C"]
const cycles_sim = 5


## ===== Chargement des données =====

raw_data = Dict{String, Vector{NamedTuple}}()

for rep in replicates
    entries    = NamedTuple[]
    t_offset   = nothing
    t_end_prev = nothing

    for cyc_idx in 1:cycles_sim
        path_H = "modelisation/input/xp_input_20/hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc_idx).csv"
        path_V = "modelisation/input/xp_input_20/virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc_idx).csv"
        (!isfile(path_H) || !isfile(path_V)) && continue

        df_H = CSV.read(path_H, DataFrame)
        df_V = CSV.read(path_V, DataFrame)

        tH = df_H[:, 1] ./ 24.0;  H = Vector{Float64}(df_H[:, 2])
        tV = df_V[:, 1] ./ 24.0;  V = Vector{Float64}(df_V[:, 2])

        if t_offset === nothing; t_offset = tH[1]; end
        tH .-= t_offset;  tV .-= t_offset

        if t_end_prev !== nothing
            gap = tH[1] - t_end_prev
            tH .-= gap;  tV .-= gap
        end
        t_end_prev = tH[end]

        push!(entries, (index=cyc_idx, tH=tH, H=H, tV=tV, V=V))
    end
    raw_data[rep] = entries
end

n_cycles_loaded = minimum(length(raw_data[rep]) for rep in replicates)
println("Cycles chargés : $n_cycles_loaded  (demandés : $cycles_sim)")


## ===== Définitions des modèles =====

isoutofdomain = (u, p, t) -> any(x -> x < 0, u)

function solve_cycle(ode!, p, u0, t0, t1)
    prob = ODEProblem(ode!, u0, (t0, t1), p)
    return solve(prob, Rodas5(), reltol=1e-6, abstol=1e-6,
                 isoutofdomain=isoutofdomain)
end

# ── SRVi ────────────────────────────────────────────────────────────────────
function SRVi!(dY, Y, p, t)
    r, K, β, δ, α, φ = p.r, p.K, p.β, p.δ, p.α, p.φ
    S, R, Vi = Y[1], Y[2], Y[3]
    H = S + R
    dY[1] = r*S*(1 - H/K) - φ*S*Vi - α*S
    dY[2] = r*R*(1 - H/K) + α*S
    dY[3] = β*φ*S*Vi - δ*Vi
end

host_sum_SRVi(u)       = u[1] + u[2]
virus_sum_SRVi(u)      = u[3]
s_val_SRVi(u)          = u[1]
vi_val_SRVi(u)         = u[3]
phi_equiv_SRVi(u, p)   = p.φ * u[1] / (u[1] + u[2])
mk_u0_SRVi(H0, V0, pS, pVi) = [pS * H0, (1 - pS) * H0, V0]

# ── SViVd ────────────────────────────────────────────────────────────────────
function SViVd!(dY, Y, p, t)
    r, K, β, δ, φ, ε = p.r, p.K, p.β, p.δ, p.φ, p.ε
    S, Vi, Vd = Y[1], Y[2], Y[3]
    dY[1] = r*S*(1 - S/K) - φ*S*Vi
    dY[2] = ε*β*φ*S*Vi - δ*Vi
    dY[3] = (1 - ε)*β*φ*S*Vi - δ*Vd
end

host_sum_SViVd(u)       = u[1]
virus_sum_SViVd(u)      = u[2] + u[3]
s_val_SViVd(u)          = u[1]
vi_val_SViVd(u)         = u[2]
phi_equiv_SViVd(u, p)   = p.φ * u[2] / (u[2] + u[3])
mk_u0_SViVd(H0, V0, pS, pVi) = [H0, pVi * V0, (1 - pVi) * V0]


## ===== Dispatch =====

if MODEL_SELECT == "SRVi"
    ode!      = SRVi!
    params    = SRVi_params
    mk_u0     = mk_u0_SRVi
    host_sum  = host_sum_SRVi
    virus_sum = virus_sum_SRVi
    s_val     = s_val_SRVi
    vi_val    = vi_val_SRVi
    phi_equiv = phi_equiv_SRVi

    # panel=1 → hôtes, panel=2 → virus
    state_labels = [
        (expr=u->u[1]+u[2], label="Model H",  color=:black, lw=4, ls=:solid, panel=1),
        (expr=u->u[1],       label="Model S",  color=:red,   lw=2, ls=:dash,  panel=1),
        (expr=u->u[2],       label="Model R",  color=:green, lw=2, ls=:dash,  panel=1),
        (expr=u->u[3],       label="Model V",  color=:black, lw=4, ls=:solid, panel=2),
        (expr=u->u[3],       label="Model Vi", color=:red,   lw=2, ls=:dash,  panel=2),
    ]

elseif MODEL_SELECT == "SViVd"
    ode!      = SViVd!
    params    = SViVd_params
    mk_u0     = mk_u0_SViVd
    host_sum  = host_sum_SViVd
    virus_sum = virus_sum_SViVd
    s_val     = s_val_SViVd
    vi_val    = vi_val_SViVd
    phi_equiv = phi_equiv_SViVd

    state_labels = [
        (expr=u->u[1],      label="Model H",  color=:black, lw=4, ls=:solid, panel=1),
        (expr=u->u[1],      label="Model S",  color=:red,   lw=2, ls=:dash,  panel=1),
        (expr=u->u[2]+u[3], label="Model V",  color=:black, lw=4, ls=:solid, panel=2),
        (expr=u->u[2],      label="Model Vi", color=:red,   lw=2, ls=:dash,  panel=2),
        (expr=u->u[3],      label="Model Vd", color=:green, lw=2, ls=:dash,  panel=2),
    ]

else
    error("MODEL_SELECT inconnu : \"$MODEL_SELECT\". Choisir \"SRVi\" ou \"SViVd\".")
end


## ===== Figure =====

replicate_colors = Dict(
    "A" => RGB(0.6, 0.8, 1.0),
    "B" => RGB(31/255, 119/255, 180/255),
    "C" => RGB(0.0, 0.3, 0.7)
)

ytick_vals1   = [10.0^i for i in -2:2:8]
ytick_labels1 = [L"10^{%$i}" for i in -2:2:8]
ytick_vals2   = [10.0^i for i in 2:2:10]
ytick_labels2 = [L"10^{%$i}" for i in 2:2:10]
ytick_vals3   = [10.0^i for i in -14:2:-6]
ytick_labels3 = [L"10^{%$i}" for i in -14:2:-6]

pl = plot(
    layout             = (3, 1),
    size               = (1800, 2000),
    left_margin        = 15mm,
    right_margin       = 10mm,
    top_margin         = 5mm,
    bottom_margin      = 10mm,
    grid               = true,
    yscale             = :log10,
    xlims              = (0, 67),
    ytickfontsize      = 22,
    legendfontsize     = 15,
    guidefontsize      = 20,
    xtickfontsize      = 20,
    titlefontsize      = 20,
    xlabel             = "Time (days)",
    legend             = :bottomright,
    plot_title         = "Model $MODEL_SELECT",
    plot_titlefontsize = 25,
)


## ===== Données expérimentales =====

for rep in replicates
    for cyc_idx in 1:n_cycles_loaded
        data = raw_data[rep][cyc_idx]
        lbl  = cyc_idx == 1 ? "Replicate $rep" : ""
        scatter!(pl[1], data.tH, data.H;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Host\n(cell/mL)",
            ylims=(1e-2, 3e8), yticks=(ytick_vals1, ytick_labels1),
            markershape=:circle, markersize=8)
        scatter!(pl[2], data.tV, data.V;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Virus\n(part/mL)",
            ylims=(1e2, 1e10), yticks=(ytick_vals2, ytick_labels2),
            markershape=:circle, markersize=8, legend=(0.13, 0.45))
    end
end


## ===== Simulation cycle par cycle =====

global prop_S  = 1.0
global prop_Vi = 1.0

for cycle in 1:n_cycles_loaded
    H0 = mean(raw_data[rep][cycle].H[1] for rep in replicates)
    V0 = mean(raw_data[rep][cycle].V[1] for rep in replicates)
    u0 = mk_u0(H0, V0, prop_S, prop_Vi)

    t0 = minimum(raw_data[rep][cycle].tH[1]   for rep in replicates)
    t1 = maximum(max(raw_data[rep][cycle].tH[end],
                     raw_data[rep][cycle].tV[end]) for rep in replicates)

    sol = solve_cycle(ode!, params, u0, t0, t1)

    sol.retcode != SciMLBase.ReturnCode.Success &&
        @warn "Cycle $cycle : solveur a échoué ($(sol.retcode))"

    # --- Propagation vers le cycle suivant ---
    u_end  = sol.u[end]
    H_end  = host_sum(u_end)
    V_end  = virus_sum(u_end)
    global prop_S  = H_end > 0 ? clamp(s_val(u_end)  / H_end, 0.0, 1.0) : 1.0
    global prop_Vi = V_end > 0 ? clamp(vi_val(u_end) / V_end, 0.0, 1.0) : 1.0

    lbl(s) = cycle == 1 ? s : ""

    # Panel 1 : hôtes
    for sl in filter(s -> s.panel == 1, state_labels)
        plot!(pl[1], sol.t, [sl.expr(u) for u in sol.u];
              label=lbl(sl.label), color=sl.color, lw=sl.lw, ls=sl.ls)
    end

    # Panel 2 : virus
    for sl in filter(s -> s.panel == 2, state_labels)
        plot!(pl[2], sol.t, [sl.expr(u) for u in sol.u];
              label=lbl(sl.label), color=sl.color, lw=sl.lw, ls=sl.ls)
    end

    # Panel 3 : φ_equiv
    phi_app = [phi_equiv(u, params) for u in sol.u]
    plot!(pl[3], sol.t, phi_app;
          label=lbl("φ_equiv"), color=:red, lw=3,
          ylabel="φ_equiv\n(mL/(cell·day))",
          ylims=(1e-14, 1e-6),
          title="φ_equiv",
          yticks=(ytick_vals3, ytick_labels3))

    println("Cycle $cycle : t=$t0→$t1  u0=$(round.(u0, sigdigits=3))  prop_S=$(round(prop_S, sigdigits=3))  prop_Vi=$(round(prop_Vi, sigdigits=3))")
end

display(pl)

fig_path = joinpath(@__DIR__, "output/simulation_only/$(MODEL_SELECT)/simulation.png")
savefig(pl, fig_path)
println("\nFigure sauvegardée : $fig_path")