## Import packages 
using Dates
using CSV
using DataFrames
using Optim
using DifferentialEquations
using LogExpFunctions
using Statistics
using Plots
using Measures
using BlackBoxOptim
using Base.Threads
using DataInterpolations
using Printf
using LaTeXStrings

## ===== Constants =====
r = 0.5592225270686286
K = 7.29695252684594e7
β = 144
δ = 0.02

replicates    = ["A", "B", "C"]
n_cycles      = 5
cycle_changes = [24.770833333333332, 34.4375, 43.104166666666664, 55.854166666666664]

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

## ===== Input =====
cycles = Dict{Tuple{String,Int}, NamedTuple}()

for rep in replicates
    t_H_prev_end = nothing
    for cycle in 1:n_cycles
        df_H = CSV.read(joinpath(@__DIR__, "../input/xp_input_20/hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cycle).csv"), DataFrame)
        df_V = CSV.read(joinpath(@__DIR__, "../input/xp_input_20/virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cycle).csv"), DataFrame)

        t_H = df_H[:, 1] ./ 24
        H   = df_H[:, 2]
        t_V = df_V[:, 1] ./ 24
        V   = df_V[:, 2]

        if t_H_prev_end !== nothing
            shift = t_H[1] - t_H_prev_end
            t_H = t_H .- shift
            t_V = t_V .- shift
        end
        t_H_prev_end = t_H[end]

        cycles[(rep, cycle)] = (tH=t_H, H=H, tV=t_V, V=V, u0=[H[1], V[1], 0])
    end
end

## ===== ODE MODELS =====
function model_SRVi(dY, Y, p, t)
    ϕ, gamma = p
    S  = Y[1]
    R  = Y[2]
    Vi = Y[3]
    dY[1] = r*S*(1 - (S+R)/K) - ϕ*S*Vi - gamma*S
    dY[2] = r*R*(1 - (S+R)/K) + gamma*S   # note: dépend de S/K pas (S+R)/K ?
    dY[3] = β*ϕ*S*Vi - δ*Vi
end

function model_SViVd(dY, Y, p, t)
    ϕ, epsilon = p
    S  = Y[1]
    Vi = Y[2]
    Vd = Y[3]
    dY[1] = r*S*(1 - S/K) - ϕ*S*Vi
    dY[2] = epsilon*β*ϕ*S*Vi - δ*Vi
    dY[3] = (1-epsilon)*β*ϕ*S*Vi - δ*Vd
end

## ===== INFOS DICT =====
infos = Dict(
    "SRVi" => (
        best_error  = 0.668845255701089,
        best_params = Dict(:phi => 5.953630401798023e-9, :gamma => 1.9631603447397458e-5),
        ode_fn      = model_SRVi,
        p_vec       = p -> [p[:phi], p[:gamma]],
        u0_fn       = (H0, V0, frac) -> [H0 * frac[1], H0 * frac[2], V0],
        frac_init   = [1.0, 0.0],   # [S_frac, R_frac]
        frac_update = sol -> begin
            S_tf = sol[1, end]; R_tf = sol[2, end]
            tot  = S_tf + R_tf
            tot > 0 ? [S_tf/tot, R_tf/tot] : [1.0, 0.0]
        end,
        phi_eff = (sol, t, p) -> begin
            S     = max(sol(t)[1], 1e-30)
            R     = max(sol(t)[2], 0.0)
            denom = S + R
            denom > 0 ? p[:phi] * S / denom : p[:phi]
        end
    ),
    "SViVd" => (
        best_error  = 24.449433952075378,
        best_params = Dict(:phi => 1.6269091475117133e-7, :epsilon => 1.4647654308910798e-7),
        ode_fn      = model_SViVd,
        p_vec       = p -> [p[:phi], p[:epsilon]],
        u0_fn       = (H0, V0, frac) -> [H0, V0 * frac[1], V0 * frac[2]],
        frac_init   = [1.0, 0.0],   # [Vi_frac, Vd_frac]
        frac_update = sol -> begin
            Vi_tf = sol[2, end]; Vd_tf = sol[3, end]
            tot   = Vi_tf + Vd_tf
            tot > 0 ? [Vi_tf/tot, Vd_tf/tot] : [1.0, 0.0]
        end,
        phi_eff = (sol, t, p) -> begin
            Vi    = max(sol(t)[2], 1e-30)
            Vd    = max(sol(t)[3], 0.0)
            denom = Vi + Vd
            denom > 0 ? p[:phi] * Vi / denom : p[:phi]
        end
    ),
)

## ===== PARSE POLYNOME =====
function parse_polynome(filepath::String)
    lines = readlines(filepath)

    function extract_float(s)
        m = match(r"([+-])?\s*([0-9]+\.[0-9]+(?:[eE][+-]?[0-9]+)?)", s)
        sign = (m.captures[1] == "-") ? "-" : ""
        return parse(Float64, sign * m.captures[2])
    end

    intervals = []
    i = 1
    while i <= length(lines)
        m = match(r"Interval\s+\[([0-9eE+\-.]+),\s*([0-9eE+\-.]+)\]", lines[i])
        if m !== nothing
            t0 = parse(Float64, m.captures[1])
            a  = extract_float(split(lines[i+1],"=")[2])
            b  = extract_float(split(lines[i+2],"*")[1])
            c  = extract_float(split(lines[i+3],"*")[1])
            d  = extract_float(split(lines[i+4],"*")[1])
            push!(intervals, (t0,a,b,c,d))
            i += 5
        else
            i += 1
        end
    end

    function φ(t)
        for i in 1:length(intervals)-1
            t0 = intervals[i][1]
            t1 = intervals[i+1][1]
            if t ≥ t0 && t < t1
                dt = t - t0
                a,b,c,d = intervals[i][2:end]
                return exp(a + b*dt + c*dt^2 + d*dt^3)
            end
        end
        a,b,c,d = intervals[end][2:end]
        dt = t - intervals[end][1]
        return exp(a + b*dt + c*dt^2 + d*dt^3)
    end

    return φ, intervals
end

n_interv  = 3
poly_file = joinpath(@__DIR__, "240426_output/knots_dilutions_$(n_interv)_polynome_3rep.txt")
φ_global, intervals = parse_polynome(poly_file)

## ===== PLOT =====
mkpath(joinpath(@__DIR__, "030526_output"))

color_phi_import = RGB(255/255, 127/255, 14/255)
model_colors = Dict(
    "SRVi"  => RGB(0.2, 0.6, 0.2),
    "SViVd" => RGB(214/255, 39/255,  40/255),
)

ytick_vals   = [10.0^i for i in -50:5:-8]
ytick_labels = [L"10^{%$i}" for i in -50:5:-8]

pl_phi = plot(
    size=(1800, 1800),
    left_margin=15mm, right_margin=10mm,
    top_margin=10mm,  bottom_margin=10mm,
    grid=true,
    xlabel="Time (days)",
    ylabel="ϕ and ϕ_effective\n(mL/(cell.day))",
    legend=:bottomright,
    yscale=:log10,
    xlims=(0, 67),
    ylims=(1e-50, 1e-8),
    ytickfontsize=26, legendfontsize=20,
    guidefontsize=22, xtickfontsize=22,
    yticks=(ytick_vals, ytick_labels)
)

# ---- phi importé ----
t_global     = collect(range(0.0, 67.0, length=2000))
phi_imported = [φ_global(t) for t in t_global]
plot!(pl_phi, t_global, phi_imported,
    label="ϕ(t)",
    color=color_phi_import,
    lw=6
)

t_knots = [it[1] for it in intervals]
scatter!(pl_phi, t_knots, φ_global.(t_knots), label="", color=color_phi_import, markersize=12, markerstrokewidth=0)


# ---- phi_eff par modèle ----
for (model_name, info) in infos
    p        = info.best_params
    phi_fn   = info.phi_eff
    ode_fn   = info.ode_fn
    p_vec    = info.p_vec(p)
    frac     = copy(info.frac_init)

    for cycle in 1:n_cycles
        H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
        V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
        u0_mean = info.u0_fn(H0_mean, V0_mean, frac)

        t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
        t1 = maximum(max(cycles[(rep, cycle)].tH[end], cycles[(rep, cycle)].tV[end]) for rep in replicates)

        sol_c = solve(
            ODEProblem(ode_fn, u0_mean, (t0, t1), p_vec),
            Rodas5(), reltol=1e-6, abstol=1e-6,
            saveat=collect(range(t0, t1, length=500)),
            isoutofdomain=isoutofdomain
        )

        frac = info.frac_update(sol_c)

        phi_eff_vals = [phi_fn(sol_c, t, p) for t in sol_c.t]

        lbl = cycle == 1 ? "ϕ_effective $model_name" : ""
        plot!(pl_phi, sol_c.t, phi_eff_vals,
            label=lbl,
            color=model_colors[model_name],
            lw=6, alpha=0.85
        )
    end
end

# ---- Lignes de dilution ----
for t_change in cycle_changes
    vline!(pl_phi, [t_change],
        color=:gray, linestyle=:dot, lw=5,
        label = t_change == cycle_changes[1] ? "Dilution" : nothing
    )
end

display(pl_phi)
savefig(pl_phi, joinpath(@__DIR__, "030526_output/phi_eff_comparison.png"))