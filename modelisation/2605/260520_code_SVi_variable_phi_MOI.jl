## ===== Import packages =====
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
using Sundials

println(Threads.nthreads())


## ===== Infos =====
n_cycles = 5

color_A = RGB(0.6, 0.8, 1.0)
color_B = RGB(31/255, 119/255, 180/255)
color_C = RGB(0.0, 0.3, 0.7) 
model_color = RGB(255/255, 127/255, 14/255)
data_color = RGB(31/255, 119/255, 180/255)
cycle_colors = [
    RGB(0.95, 0.45, 0.45),  # rose
    RGB(0.95, 0.70, 0.30),  # pêche
    RGB(0.40, 0.78, 0.40),  # vert
    RGB(0.35, 0.60, 0.95),  # bleu
    RGB(0.70, 0.45, 0.95),  # lavande
]
cycle_changes = [24.770833333333332, 34.4375, 43.104166666666664, 55.854166666666664]

replicate_colors = [color_A, color_B, color_C]
replicates = ["A", "B", "C"]

dir_output = "260520_output_MOI_SVi_variable_phi"

## ===== Input =====
cycles = Dict{Tuple{String,Int}, NamedTuple}()

t_H_all = Dict{String, Vector{Float64}}()
H_all   = Dict{String, Vector{Float64}}()
t_V_all = Dict{String, Vector{Float64}}()
V_all   = Dict{String, Vector{Float64}}()

for rep in replicates
    t_H_rep = Float64[]
    H_rep   = Float64[]
    t_V_rep = Float64[]
    V_rep   = Float64[]

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

        # Conditions initiales pour ce cycle
        u0 = [H[1], V[1]]

        cycles[(rep, cycle)] = (tH=t_H, H=H, tV=t_V, V=V, u0=u0)

        append!(t_H_rep, t_H)
        append!(H_rep,   H)
        append!(t_V_rep, t_V)
        append!(V_rep,   V)
    end

    t_H_all[rep] = t_H_rep
    H_all[rep]   = H_rep
    t_V_all[rep] = t_V_rep
    V_all[rep]   = V_rep
end


## ===== CHECK THE DATA =====
pl_data = plot(layout=(1,2), size=(900,300), margins=5mm, legend=:topright)

for (i, rep) in enumerate(replicates)
    scatter!(pl_data[1], t_H_all[rep], H_all[rep],
        color=replicate_colors[i], label="Rep $rep",
        xlabel="time (days)", ylabel="abundances (cell/ml)",
        yscale=:log10, ylims=(1e2, 1e8))

    scatter!(pl_data[2], t_V_all[rep], V_all[rep],
        color=replicate_colors[i], label="Rep $rep",
        xlabel="time (days)", ylabel="abundances (part/ml)",
        yscale=:log10, ylims=(1e3, 1e10))
end

#display(pl_data)




## ===== Constants =====
r = 0.5592225270686286
K = 7.29695252684594e7
β = 144
δ = 0.02


## ===== Import ϕ(t) =====
const poly_file = joinpath(@__DIR__, "polynome.txt")

function parse_polynome(filepath::String)
    isfile(filepath) || error("Polynomial file not found: $filepath")
    lines = readlines(filepath)

    metadata = Dict{String,Any}()
    for line in lines
        m = match(r"^(n_knots|n_intervalles|best_fitness|best_seed|n_runs)\s*=\s*(.+)$",
                  strip(line))
        m === nothing && continue
        key = m.captures[1]
        val = strip(m.captures[2])
        metadata[key] = key in ("n_knots", "n_intervalles", "best_seed", "n_runs") ?
                         parse(Int, val) : parse(Float64, val)
    end

    # Robust float extractor: keeps sign, digits, decimal point, exponent
    function extract_float(s::AbstractString)::Float64
        m = match(r"([+-]?\s*[0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)", strip(s))
        m === nothing && error("Cannot extract float from: \"$s\"")
        parse(Float64, replace(m.captures[1], " " => ""))
    end

    intervals = Vector{NTuple{6,Float64}}()
    i = 1
    while i <= length(lines)
        m_iv = match(
            r"Interval\s+\[([0-9eE+\-.]+),\s*([0-9eE+\-.]+)\]\s+days:",
            strip(lines[i]))
        if m_iv !== nothing
            # Need at least 4 more lines for coefficients
            if i + 4 > length(lines)
                @warn "Incomplete interval block at line $i — skipping"
                i += 1; continue
            end
            t0 = parse(Float64, m_iv.captures[1])
            t1 = parse(Float64, m_iv.captures[2])
            # Lines i+1..i+4: "a = ...", "b*...", "c*...", "d*..."
            a  = extract_float(split(lines[i+1], "=")[end])
            b  = extract_float(split(lines[i+2], "*")[1])
            c  = extract_float(split(lines[i+3], "*")[1])
            d  = extract_float(split(lines[i+4], "*")[1])
            push!(intervals, (t0, t1, a, b, c, d))
            i += 5; continue
        end
        i += 1
    end

    isempty(intervals) && error("No intervals found in $filepath")
    println("  Parsed $(length(intervals)) intervals from $(basename(filepath))")
    haskey(metadata, "best_fitness") &&
        println("  fitness = $(metadata["best_fitness"])")

    t_lo = intervals[1][1]
    t_hi = intervals[end][2]

    function phi_raw(t::Real)::Float64
        tc = clamp(Float64(t), t_lo, t_hi)
        # Find the interval containing tc
        idx = length(intervals)
        for k in eachindex(intervals)
            if tc <= intervals[k][2]
                idx = k; break
            end
        end
        t0, _, a, b, c, d = intervals[idx]
        dt = tc - t0
        return exp(a + b*dt + c*dt^2 + d*dt^3)
    end

    return phi_raw, metadata
end

phi_raw, phi_meta = parse_polynome(poly_file)

const N_GRID  = 10_000
const t_grid  = collect(range(t_H_all["A"][1], t_H_all["A"][end]; length=N_GRID))
const phi_grid = phi_raw.(t_grid)

# LinearInterpolation from DataInterpolations: LinearInterpolation(u, t)
const ϕ_interp = LinearInterpolation(phi_grid, t_grid)

t_check = range(t_H_all["A"][1], t_H_all["A"][end]; length=1_000)


## ===== Plot ϕ(t) =====
ytick_vals_phi   = [10.0^i for i in -14:-8]
ytick_labels_phi = [L"10^{%$i}" for i in -14:-8]

pl_interp = plot(
    size=(1800,800),
    left_margin=15mm, right_margin=10mm,
    top_margin=10mm,  bottom_margin=15mm,
    grid=true,
    yscale=:log10,
    xlims=(0, 67),
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="Time (days)",
    ylabel="ϕ (mL/(part.day))",
    legend=:bottomright,
    plot_title="Model SVi with variable phi",
    title="ϕ vs time"
)

plot!(pl_interp, collect(t_check), ϕ_interp.(collect(t_check)),
    color=model_color, lw=4, alpha=0.7, label="ϕ", ylims=(1e-14, 1e-8), yticks=(ytick_vals_phi, ytick_labels_phi))

for t_change in cycle_changes
    vline!(pl_interp, [t_change],
        color=data_color, linestyle=:dot, lw=2,
        label=(t_change == cycle_changes[1] ? "Dilution" : ""))
end

display(pl_interp)
savefig(pl_interp, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_phi_vs_t.png"))


## ===== Model =====
function model!(dY, Y, p, t)
    S, V = Y[1], Y[2]
    ϕt = ϕ_interp(clamp(t, t_grid[1], t_grid[end]))
    dY[1] = r*S*(1 - S/K) - ϕt*S*V
    dY[2] = β*ϕt*S*V - δ*V
end

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)


## ===== Plot S=f(t) and V=g(t) =====
ytick_vals1   = [10.0^i for i in 2:8]
ytick_labels1 = [L"10^{%$i}" for i in 2:8]

ytick_vals2   = [10.0^i for i in 5:10]
ytick_labels2 = [L"10^{%$i}" for i in 5:10]

pl_fit = plot(
    layout=(2,1),
    size=(1800,1100),
    left_margin=15mm, right_margin=10mm,
    top_margin=10mm,  bottom_margin=10mm,
    grid=true,
    yscale=:log10,
    xlims=(0, 67),
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="Time (days)",
    legend=:bottomright,
    plot_title="Model SVi with variable phi",
    title="S and V vs time"
)

# Scatter des données
for (i, rep) in enumerate(replicates)
    for cycle in 1:n_cycles
        cyc = cycles[(rep, cycle)]
        lbl = cycle == 1 ? "Replicate $rep" : ""
        scatter!(pl_fit[1], cyc.tH, cyc.H,
            label=lbl, color=replicate_colors[i], alpha=0.7,
            ylabel="Host abundance\n(cell/mL)",
            ylims=(1e2, 3e8), yticks=(ytick_vals1, ytick_labels1),
            markershape=:circle, markersize=8)
        scatter!(pl_fit[2], cyc.tV, cyc.V,
            label=lbl, color=replicate_colors[i], alpha=0.7,
            ylabel="Virus abundance\n(part/mL)",
            ylims=(1e5, 1e10), yticks=(ytick_vals2, ytick_labels2),
            markershape=:circle, markersize=8, legend=(0.13, 0.4))
    end
end

# Modèle : un segment par cycle, CI = moyenne des réplicats
for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1]   for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    lbl = cycle == 1 ? "Model" : ""
    plot!(pl_fit[1], sol_c.t, sol_c[1,:], label=lbl, color=model_color, lw=4, alpha=0.7)
    plot!(pl_fit[2], sol_c.t, sol_c[2,:], label=lbl, color=model_color, lw=4, alpha=0.7)
end

for t_change in cycle_changes
    for i in 1:2
        vline!(pl_fit[i], [t_change],
            color=data_color, linestyle=:dot, lw=2,
            label=(t_change == cycle_changes[1] ? "Dilution" : ""))
    end
end

display(pl_fit)
savefig(pl_fit, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_S_Vi_vs_t.png"))


## ===== Plot MOI=h(t) =====
ytick_vals3   = [10.0^i for i in -1:6]
ytick_labels3 = [L"10^{%$i}" for i in -1:6]

pl_MOI = plot(
    layout=(1,1),
    size=(1100,600),
    left_margin=15mm, right_margin=10mm,
    top_margin=10mm,  bottom_margin=10mm,
    grid=true,
    yscale=:log10,
    xlims=(0, 67),
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="Time (days)",
    ylabel="MOI (Vi/S)",
    legend=:topright,
    plot_title="Model SVi with variable phi",
    title="MOI vs time"
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    MOI = sol_c[2,:] ./ sol_c[1,:]
    lbl = cycle == 1 ? "MOI" : ""
    plot!(pl_MOI, sol_c.t, MOI, label=lbl, color=model_color, lw=4, alpha=0.7, ylims=(1e-1, 1e6), yticks=(ytick_vals3, ytick_labels3))
end

for t_change in cycle_changes
    vline!(pl_MOI, [t_change],
        color=data_color, linestyle=:dot, lw=2,
        label=(t_change == cycle_changes[1] ? "Dilution" : ""))
end

display(pl_MOI)
savefig(pl_MOI, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_MOI_vs_t.png"))


## ===== Plot MOI vs dV/dt / V =====
pl_MOI_rate = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="1/V dV/dt = βϕ(t)S−δ",
    ylabel="MOI (Vi/S)",
    yscale=:log10,
    legend=:topright,
    plot_title="Model SVi with variable phi",
    title="MOI vs virus growth rate",
    ylims=(1e-1, 1e6), yticks=(ytick_vals3, ytick_labels3)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec
    V_rate   = β .* ϕ_vec .* S_vec .- δ   # 1/V * dV/dt

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl_MOI_rate, V_rate, MOI,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl_MOI_rate, [V_rate[1]], [MOI[1]],
        label="",
        color=col,
        markersize=5,
        markershape=:circle,
        markerstrokewidth=0,
        alpha=1.0)
end

display(pl_MOI_rate)
savefig(pl_MOI_rate, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_MOI_vs_Vrate.png"))


## ===== Plot S vs MOI =====
xtick_vals   = [10.0^i for i in -1:6]
xtick_labels = [L"10^{%$i}" for i in -1:6]
ytick_vals4   = [10.0^i for i in 3:8]
ytick_labels4 = [L"10^{%$i}" for i in 3:8]

pl = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="MOI (Vi/S)",
    ylabel="S",
    xscale=:log10,
    yscale=:log10,
    legend=:topright,
    plot_title="Model SVi with variable phi",
    title="S vs MOI",
    xlims=(1e-1, 1e6), xticks=(xtick_vals, xtick_labels),
    ylims=(1e3, 1e8), yticks=(ytick_vals4, ytick_labels4)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl, MOI, S_vec,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl, [MOI[1]], [S_vec[1]],
    label="",
    color=col,
    markersize=5,
    markershape=:circle,
    markerstrokewidth=0,
    alpha=1.0)
end

display(pl)
savefig(pl, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_S_vs_MOI.png"))


## ===== Plot V vs MOI =====
xtick_vals   = [10.0^i for i in -1:6]
xtick_labels = [L"10^{%$i}" for i in -1:6]
ytick_vals4   = [10.0^i for i in 6:9]
ytick_labels4 = [L"10^{%$i}" for i in 6:9]

pl = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="MOI (Vi/S)",
    ylabel="V",
    xscale=:log10,
    yscale=:log10,
    legend=:bottomright,
    plot_title="Model SVi with variable phi",
    title="V vs MOI",
    xlims=(1e-1, 1e6), xticks=(xtick_vals, xtick_labels),
    ylims=(1e6, 3e9), yticks=(ytick_vals4, ytick_labels4)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl, MOI, V_vec,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl, [MOI[1]], [V_vec[1]],
    label="",
    color=col,
    markersize=5,
    markershape=:circle,
    markerstrokewidth=0,
    alpha=1.0)
end

display(pl)
savefig(pl, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_V_vs_MOI.png"))


## ===== Plot ϕ vs MOI =====
xtick_vals   = [10.0^i for i in -1:6]
xtick_labels = [L"10^{%$i}" for i in -1:6]
ytick_vals4   = [10.0^i for i in -14:-8]
ytick_labels4 = [L"10^{%$i}" for i in -14:-8]

pl = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="MOI (Vi/S)",
    ylabel="ϕ",
    xscale=:log10,
    yscale=:log10,
    legend=:bottomleft,
    plot_title="Model SVi with variable phi",
    title="ϕ vs MOI",
    xlims=(1e-1, 1e6), xticks=(xtick_vals, xtick_labels),
    ylims=(1e-14, 1e-8), yticks=(ytick_vals4, ytick_labels4)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl, MOI, ϕ_vec,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl, [MOI[1]], [ϕ_vec[1]],
    label="",
    color=col,
    markersize=5,
    markershape=:circle,
    markerstrokewidth=0,
    alpha=1.0)
end

display(pl)
savefig(pl, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_phi_vs_MOI.png"))


## ===== Plot ϕSV vs MOI =====
xtick_vals   = [10.0^i for i in -1:6]
xtick_labels = [L"10^{%$i}" for i in -1:6]
ytick_vals4   = [10.0^i for i in 1:7]
ytick_labels4 = [L"10^{%$i}" for i in 1:7]

pl = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="MOI (Vi/S)",
    ylabel="ϕSV",
    xscale=:log10,
    yscale=:log10,
    legend=:topright,
    plot_title="Model SVi with variable phi",
    title="ϕSV vs MOI",
    xlims=(1e-1, 1e6), xticks=(xtick_vals, xtick_labels),
    ylims=(1e1, 1e7), yticks=(ytick_vals4, ytick_labels4)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl, MOI, ϕ_vec.*S_vec.*V_vec,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl, [MOI[1]], [ϕ_vec[1]*S_vec[1]*V_vec[1]],
    label="",
    color=col,
    markersize=5,
    markershape=:circle,
    markerstrokewidth=0,
    alpha=1.0)
end

display(pl)
savefig(pl, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_phiSV_vs_MOI.png"))


## ===== Plot βϕSV vs MOI =====
xtick_vals   = [10.0^i for i in -1:6]
xtick_labels = [L"10^{%$i}" for i in -1:6]
ytick_vals4   = [10.0^i for i in 3:9]
ytick_labels4 = [L"10^{%$i}" for i in 3:9]

pl = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="MOI (Vi/S)",
    ylabel="βϕSV",
    xscale=:log10,
    yscale=:log10,
    legend=:topright,
    plot_title="Model SVi with variable phi",
    title="βϕSV vs MOI",
    xlims=(1e-1, 1e6), xticks=(xtick_vals, xtick_labels),
    ylims=(1e3, 1e9), yticks=(ytick_vals4, ytick_labels4)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl, MOI, β.*ϕ_vec.*S_vec.*V_vec,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl, [MOI[1]], [β.*ϕ_vec[1]*S_vec[1]*V_vec[1]],
    label="",
    color=col,
    markersize=5,
    markershape=:circle,
    markerstrokewidth=0,
    alpha=1.0)
end

display(pl)
savefig(pl, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_betaphiSV_vs_MOI.png"))


## ===== Plot dS/dt vs MOI =====
xtick_vals   = [10.0^i for i in -1:6]
xtick_labels = [L"10^{%$i}" for i in -1:6]
ytick_vals4   = [10.0^i for i in -7:7]
ytick_labels4 = [L"10^{%$i}" for i in -7:7]

pl = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="MOI (Vi/S)",
    ylabel="dS/dt",
    xscale=:log10,
    legend=:topright,
    plot_title="Model SVi with variable phi",
    title="dS/dt vs MOI",
    xlims=(1e-1, 1e6), xticks=(xtick_vals, xtick_labels),
    #ylims=(-2e7, 2e7), yticks=(ytick_vals4, ytick_labels4)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec

    dS = r.*S_vec.*(1 .-S_vec./K) .- ϕ_vec.*S_vec.*V_vec
    println(dS)

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl, MOI, dS,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl, [MOI[1]], [dS[1]],
    label="",
    color=col,
    markersize=5,
    markershape=:circle,
    markerstrokewidth=0,
    alpha=1.0)
end

display(pl)
savefig(pl, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_dS_vs_MOI.png"))


## ===== Plot 1/S dS/dt vs MOI =====
xtick_vals   = [10.0^i for i in -1:6]
xtick_labels = [L"10^{%$i}" for i in -1:6]
ytick_vals4   = [10.0^i for i in -7:7]
ytick_labels4 = [L"10^{%$i}" for i in -7:7]

pl = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="MOI (Vi/S)",
    ylabel="1/S dS/dt",
    xscale=:log10,
    legend=:bottomleft,
    plot_title="Model SVi with variable phi",
    title="1/S dS/dt vs MOI",
    xlims=(1e-1, 1e6), xticks=(xtick_vals, xtick_labels),
    #ylims=(-2e7, 2e7), yticks=(ytick_vals4, ytick_labels4)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec

    dS = r.*(1 .-S_vec./K) .- ϕ_vec.*V_vec
    println(dS)

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl, MOI, dS,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl, [MOI[1]], [dS[1]],
    label="",
    color=col,
    markersize=5,
    markershape=:circle,
    markerstrokewidth=0,
    alpha=1.0)
end

display(pl)
savefig(pl, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_dS_S_vs_MOI.png"))


## ===== Plot dV/dt vs MOI =====
xtick_vals   = [10.0^i for i in -1:6]
xtick_labels = [L"10^{%$i}" for i in -1:6]
ytick_vals4   = [10.0^i for i in -7:7]
ytick_labels4 = [L"10^{%$i}" for i in -7:7]

pl = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="MOI (Vi/S)",
    ylabel="dV/dt",
    xscale=:log10,
    legend=:topright,
    plot_title="Model SVi with variable phi",
    title="dV/dt vs MOI",
    xlims=(1e-1, 1e6), xticks=(xtick_vals, xtick_labels),
    #ylims=(-2e7, 2e7), yticks=(ytick_vals4, ytick_labels4)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec

    dV = β.*ϕ_vec.*S_vec.*V_vec - δ.*V_vec
    println(dV)

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl, MOI, dV,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl, [MOI[1]], [dV[1]],
    label="",
    color=col,
    markersize=5,
    markershape=:circle,
    markerstrokewidth=0,
    alpha=1.0)
end

display(pl)
savefig(pl, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_dV_vs_MOI.png"))


## ===== Plot 1/V dV/dt vs MOI =====
xtick_vals   = [10.0^i for i in -1:6]
xtick_labels = [L"10^{%$i}" for i in -1:6]
ytick_vals4   = [10.0^i for i in -7:7]
ytick_labels4 = [L"10^{%$i}" for i in -7:7]

pl = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="MOI (Vi/S)",
    ylabel="1/V dV/dt",
    xscale=:log10,
    legend=:topright,
    plot_title="Model SVi with variable phi",
    title="1/V dV/dt vs MOI",
    xlims=(1e-1, 1e6), xticks=(xtick_vals, xtick_labels),
    #ylims=(-2e7, 2e7), yticks=(ytick_vals4, ytick_labels4)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec

    dV = β.*ϕ_vec.*S_vec .- δ
    println(dV)

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl, MOI, dV,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl, [MOI[1]], [dV[1]],
    label="",
    color=col,
    markersize=5,
    markershape=:circle,
    markerstrokewidth=0,
    alpha=1.0)
end

display(pl)
savefig(pl, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_dV_V_vs_MOI.png"))


## ===== Plot rS(1-S/K) vs MOI =====
xtick_vals   = [10.0^i for i in -1:6]
xtick_labels = [L"10^{%$i}" for i in -1:6]
ytick_vals4   = [10.0^i for i in 3:8]
ytick_labels4 = [L"10^{%$i}" for i in 3:8]

pl2 = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,
    ytickfontsize=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="MOI (Vi/S)",
    ylabel="rS(1-S/K)",
    xscale=:log10,
    yscale=:log10,
    legend=:topright,
    plot_title="Model SVi with variable phi",
    title="rS(1-S/K) vs MOI",
    xlims=(1e-1, 1e6), xticks=(xtick_vals, xtick_labels),
    ylims=(1e3, 1e8), yticks=(ytick_vals4, ytick_labels4)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec

    val = r.*S_vec.*(1 .- S_vec./K)

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl2, MOI, val,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl2, [MOI[1]], [val[1]],
    label="",
    color=col,
    markersize=5,
    markershape=:circle,
    markerstrokewidth=0,
    alpha=1.0)
end

display(pl2)
savefig(pl2, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_rS1SK_vs_MOI.png"))


## ===== Plot δV vs MOI =====
xtick_vals   = [10.0^i for i in -1:6]
xtick_labels = [L"10^{%$i}" for i in -1:6]
ytick_vals4   = [10.0^i for i in 4:7]
ytick_labels4 = [L"10^{%$i}" for i in 4:7]

pl = plot(
    size=(900,700),
    left_margin=15mm, right_margin=10mm,
    top_margin=15mm,  bottom_margin=10mm,
    grid=true,çe=26, legendfontsize=17,
    guidefontsize=24, xtickfontsize=24,
    titlefontsize=18, plot_titlefontsize=24,
    xlabel="MOI (Vi/S)",
    ylabel="δV",
    xscale=:log10,
    yscale=:log10,
    legend=:bottomright,
    plot_title="Model SVi with variable phi",
    title="δV vs MOI",
    xlims=(1e-1, 1e6), xticks=(xtick_vals, xtick_labels),
    ylims=(1e4, 5e7), yticks=(ytick_vals4, ytick_labels4)
)

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end],
                     cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    S_vec   = sol_c[1,:]
    V_vec   = sol_c[2,:]
    t_vec   = sol_c.t

    ϕ_vec   = ϕ_interp.(clamp.(t_vec, t_grid[1], t_grid[end]))

    MOI      = V_vec ./ S_vec

    dV = δ.*V_vec
    println(dV)

    lbl = cycle == 1 ? "Cycle $cycle" : "Cycle $cycle"
    col = cycle_colors[cycle]

    plot!(pl, MOI, dV,
    label="Cycle $cycle",
    color=col,
    lw=4, alpha=0.7)

    scatter!(pl, [MOI[1]], [dV[1]],
    label="",
    color=col,
    markersize=5,
    markershape=:circle,
    markerstrokewidth=0,
    alpha=1.0)
end

display(pl)
savefig(pl, joinpath(@__DIR__, dir_output, "plot_SVi_variable_phi_deltaV_vs_MOI.png"))