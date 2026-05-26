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

println(Threads.nthreads())

n_runs = 10
n_knots_range = 0:2


## ===== Input =====
function load_replicate(rep::String)
    dfs = []
    for cyc in 1:5
        dfH = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv"), DataFrame)
        dfV = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv"), DataFrame)
        push!(dfs, (dfH=dfH, dfV=dfV))
    end
    return dfs
end

function build_cycles(dfs)
    cycles = []
    t_offset = 0.0

    for (i, df) in enumerate(dfs)
        tH = df.dfH[:,1]./24
        H  = df.dfH[:,2]
        tV = df.dfV[:,1]./24
        V  = df.dfV[:,2]

        if i == 1
            t_offset = tH[1]
        end

        tH = tH .- t_offset
        tV = tV .- t_offset

        if i > 1
            prev = cycles[end]
            gap = tH[1] - prev.tH[end]
            # colle les cycles bout à bout
            tH = tH .- gap
            tV = tV .- gap
            t_offset += gap
        end

        push!(cycles, (tH=tH, H=H, tV=tV, V=V, u0=[H[1], V[1]]))
    end
    return cycles
end

dfs_A = load_replicate("A")
dfs_B = load_replicate("B")
dfs_C = load_replicate("C")

cycles_A = build_cycles(dfs_A)
cycles_B = build_cycles(dfs_B)
cycles_C = build_cycles(dfs_C)

all_cycles = vcat(cycles_A, cycles_B, cycles_C)  # 15 cycles au total

general_cycles = []

    for i in 1:5

        cycA = cycles_A[i]
        cycB = cycles_B[i]
        cycC = cycles_C[i]

        ## temps (on prend A comme référence)
        tH = cycA.tH
        tV = cycA.tV

        ## moyenne des conditions initiales
        H0 = mean([cycA.u0[1], cycB.u0[1], cycC.u0[1]])
        V0 = mean([cycA.u0[2], cycB.u0[2], cycC.u0[2]])

        u0 = [H0, V0]

        push!(general_cycles, (tH=tH, tV=tV, u0=u0))
    end

# Vecteurs globaux pour les plots
t_H_A = vcat([c.tH for c in cycles_A]...); H_A = vcat([c.H for c in cycles_A]...)
t_V_A = vcat([c.tV for c in cycles_A]...); V_A = vcat([c.V for c in cycles_A]...)
t_H_B = vcat([c.tH for c in cycles_B]...); H_B = vcat([c.H for c in cycles_B]...)
t_V_B = vcat([c.tV for c in cycles_B]...); V_B = vcat([c.V for c in cycles_B]...)
t_H_C = vcat([c.tH for c in cycles_C]...); H_C = vcat([c.H for c in cycles_C]...)
t_V_C = vcat([c.tV for c in cycles_C]...); V_C = vcat([c.V for c in cycles_C]...)


## ===== CONSTANTS =====
r = 0.574619342477644
K = 6.675449070379925e7
β = 144
δ = 0.02


## ===== Bornes temporelles globales =====
t_start = minimum([c.tH[1]   for c in all_cycles])
t_end   = maximum([c.tH[end] for c in all_cycles])


## ===== MODÈLE =====
include("model_SV_2.jl")


## ===== COULEURS =====
color_A = RGB(0.6, 0.8, 1.0)
color_B = RGB(31/255, 119/255, 180/255)
color_C = RGB(0.0, 0.3, 0.7) 
model_color = RGB(255/255, 127/255,  14/255)
phi_color   = RGB(255/255, 127/255,  14/255)
cycle_color = RGBA(0.5, 0.5, 0.5, 0.4)

ytick_vals1   = [10.0^i for i in 3:1:9]
ytick_labels1 = [L"10^{%$i}" for i in 3:1:9]
ytick_vals2   = [10.0^i for i in 6:1:10]
ytick_labels2 = [L"10^{%$i}" for i in 6:1:10]
ytick_vals3   = [10.0^i for i in -15:-7]
ytick_labels3 = [L"10^{%$i}" for i in -15:-7]


## ===== ACCUMULATEURS =====
best_fitness_per_nknots = Float64[]
n_intervalles_vec       = Int[]


## ===== BOUCLE SUR n_knots =====
for n_knots in n_knots_range

    println("\n========== n_knots = $n_knots (n_intervalles = $(n_knots+1)) ==========")

    ## Nœuds de spline (basés sur les cycles du réplicat A comme référence temporelle)
    t_knots_cycles = [
        range(minimum(c.tH), maximum(c.tH), length=n_knots+2)[2:end]
        for c in cycles_A
    ]
    t_knots = vcat(t_start, vcat(t_knots_cycles...)[1:end-1], t_end)
    println(t_knots)
    total_knots = length(t_knots)
    println("  Nombre de nœuds de spline : ", total_knots)

    ## Fonction objectif — cumule l'erreur sur les 15 cycles (3 réplicats × 5 cycles)
    function objective_spline(θ)
        log_phi_spline = CubicSpline(θ, t_knots)
        t_lo = t_knots[1]
        t_hi = t_knots[end]
        phi_func(t) = exp(log_phi_spline(clamp(t, t_lo, t_hi)))

        total_err = 0.0
        for cyc in all_cycles
            t0 = cyc.tH[1]
            t1 = max(cyc.tH[end], cyc.tV[end])
            u0 = cyc.u0

            t_save = sort(unique(vcat(cyc.tH, cyc.tV)))

            prob = ODEProblem(model, u0, (t0, t1), phi_func)
            sol  = solve(prob, Rodas5(), reltol=1e-6, abstol=1e-6, saveat=t_save)

            if sol.retcode != :Success || any(u -> any(x -> !isfinite(x) || x < 0, u), sol.u)
                return 1e12
            end

            S_pred = [max(sol(t)[1], 1e-12) for t in cyc.tH]
            V_pred = [max(sol(t)[2], 1e-12) for t in cyc.tV]

            total_err += sum((log.(S_pred) .- log.(cyc.H)).^2)
            total_err += sum((log.(V_pred) .- log.(cyc.V)).^2)
        end
        return total_err
    end

    lower        = fill(log(1e-15), total_knots)
    upper        = fill(log(1e-6),  total_knots)
    search_range = [(lower[i], upper[i]) for i in eachindex(lower)]

    function run_DE(seed)
        res = bboptimize(
            objective_spline;
            SearchRange          = search_range,
            NumDimensions        = length(search_range),
            Method               = :xnes,
            PopulationSize       = 1000,
            MaxSteps             = 10000,
            DifferentialWeight   = 0.5,
            CrossoverProbability = 0.9,
            TraceMode            = :silent,
            RandomSeed           = seed
        )
        return (
            fitness = best_fitness(res),
            θ       = best_candidate(res),
            seed    = seed
        )
    end

    ## Optimisation parallèle
    results = Vector{NamedTuple{(:fitness, :θ, :seed), Tuple{Float64, Vector{Float64}, Int}}}(undef, n_runs)

    Threads.@threads for i in 1:n_runs
        println("  Thread ", threadid(), " démarre run ", i, " (seed=$(1000+i))")
        results[i] = run_DE(1000 + i)
    end

    best_idx    = argmin(r.fitness for r in results)
    best_result = results[best_idx]
    θbest       = best_result.θ
    best_seed   = best_result.seed
    best_fit    = best_result.fitness

    println("  Best fitness = $best_fit  (seed=$best_seed)")
    println("  Best phi at knots = ", exp.(θbest))

    push!(best_fitness_per_nknots, best_fit)
    push!(n_intervalles_vec,       n_knots + 1)


    ## ===== RECONSTRUCTION =====
    log_phi_spline_best = CubicSpline(θbest, t_knots)
    phi_func_best(t)    = exp(log_phi_spline_best(clamp(t, t_knots[1], t_knots[end])))


    sols = []

    for cyc in general_cycles
        t0 = cyc.tH[1]
        t1 = cyc.tH[end]

        prob = ODEProblem(model, cyc.u0, (t0, t1), phi_func_best)
        sol  = solve(prob, Rodas5(), reltol=1e-6, abstol=1e-6)

        push!(sols, sol)
    end

    ## ===== PLOT =====
    t_fine    = range(t_start, t_end, length=1000)
    phi_curve = [phi_func_best(t) for t in t_fine]

    n_int = n_knots + 1

    pl = plot(
        layout=(3,1),
        size=(2000,1500),
        left_margin=15mm, right_margin=10mm,
        top_margin=15mm,  bottom_margin=10mm,
        grid=true, yscale=:log10,
        xlims=(t_start, t_end+20),
        ytickfontsize=26, legendfontsize=20,
        guidefontsize=24, xtickfontsize=24,
        titlefontsize=24,
        xlabel="Time (days)",
        legend=:bottomright
    )

    plot!(pl, title="n_intervalles=$(n_int)  |  fitness=$(round(best_fit, sigdigits=5))  |  seed=$(best_seed)",
          subplot=1, titlefontsize=22)

    # Délimitations des cycles (basées sur réplicat A)
    cycle_starts = [c.tH[1] for c in cycles_A[2:end]]
    for sp in 1:3
        vline!(pl[sp], cycle_starts, lw=1, linestyle=:dash, color=cycle_color, label="Dilution")
    end

    # Données — Host
    scatter!(pl[1], t_H_A, H_A, label="Rep A", color=color_A,
             markershape=:circle, markersize=6, markerstrokewidth=0,
             ylabel="Host abundance\n(cell/mL)", ylims=(1e3,1e9),
             yticks=(ytick_vals1, ytick_labels1))
    scatter!(pl[1], t_H_B, H_B, label="Rep B", color=color_B,
             markershape=:circle, markersize=6, markerstrokewidth=0)
    scatter!(pl[1], t_H_C, H_C, label="Rep C", color=color_C,
             markershape=:circle, markersize=6, markerstrokewidth=0)

    # Données — Virus
    scatter!(pl[2], t_V_A, V_A, label="Rep A", color=color_A,
             markershape=:circle, markersize=6, markerstrokewidth=0,
             ylabel="Virus abundance\n(virion/mL)", ylims=(1e6,1e10),
             yticks=(ytick_vals2, ytick_labels2))
    scatter!(pl[2], t_V_B, V_B, label="Rep B", color=color_B,
             markershape=:circle, markersize=6, markerstrokewidth=0)
    scatter!(pl[2], t_V_C, V_C, label="Rep C", color=color_C,
             markershape=:circle, markersize=6, markerstrokewidth=0)

    # Courbes modèle
    for (i, sol) in enumerate(sols)
        plot!(pl[1], sol.t, sol[1,:], lw=4, color=model_color, label= i==1 ? "Model A" : nothing)
        plot!(pl[2], sol.t, sol[2,:], lw=4, color=model_color, label= i==1 ? "Model A" : nothing)
    end

    # φ(t)
    scatter!(pl[3], t_knots, exp.(θbest),
        label=" Knots", markershape=:circle, markersize=8,
        markerstrokewidth=0, color=model_color)
    plot!(pl[3], t_fine, phi_curve,
        label=" φ(t) spline", lw=4,
        ylabel="φ (ml/(cell.day))",
        color=phi_color,
        ylims=(1e-15, 1e-7),
        yticks=(ytick_vals3, ytick_labels3))

    plot_path = joinpath(@__DIR__, "240426_output/knots_dilutions_$(n_int)_plot_3rep.png")
    println("  Saving plot → $plot_path")
    savefig(pl, plot_path)


    ## ===== EXPORT POLYNÔME =====
    T = log_phi_spline_best.t
    A_sp = log_phi_spline_best.u
    z = log_phi_spline_best.z

    output_file = joinpath(@__DIR__, "240426_output/knots_dilutions_$(n_int)_polynome_3rep.txt")
    open(output_file, "w") do io
        @printf(io, "========== METADATA ==========\n")
        @printf(io, "n_knots       = %d\n",   n_knots)
        @printf(io, "n_intervalles = %d\n",   n_int)
        @printf(io, "best_fitness  = %.8e\n", best_fit)
        @printf(io, "best_seed     = %d\n",   best_seed)
        @printf(io, "n_runs        = %d\n\n", n_runs)
        @printf(io, "replicates    = A, B, C (joint fit)\n\n")

        @printf(io, "========== POLYNOMIAL DESCRIPTION OF φ(t) ==========\n")
        @printf(io, "φ(t) piecewise cubic on %d intervals\n\n", length(T)-1)

        for i in 1:(length(T)-1)
            ti  = T[i]; ti1 = T[i+1]; h = ti1 - ti
            a = A_sp[i]
            b = (A_sp[i+1] - A_sp[i])/h - h*(2*z[i] + z[i+1])/6
            c = z[i]/2
            d = (z[i+1] - z[i])/(6*h)

            @printf(io, "Interval [%.4f, %.4f] days:\n", ti, ti1)
            @printf(io, "  log(φ(t)) = %.6e\n", a)
            @printf(io, "       + %.6e * (t - %.4f)\n", b, ti)
            @printf(io, "       + %.6e * (t - %.4f)^2\n", c, ti)
            @printf(io, "       + %.6e * (t - %.4f)^3\n\n", d, ti)
        end
    end
    println("  Polynome exported → $output_file")

end  # fin boucle n_knots