using CSV
using DataFrames
using DifferentialEquations
using SciMLBase
using OrdinaryDiffEqRosenbrock
using DataInterpolations
using Statistics
using Plots
using Measures
using LaTeXStrings
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# plot_fit.jl
#
# Visualise, à partir d'un fichier polynomial.txt :
#   1. φ(t)       : la fonction de taux d'infection au cours du temps
#   2. S(t), V(t) : le fit du modèle ODE sur les données expérimentales
#
# Usage :
#   julia plot_fit.jl                              # cherche polynomial.txt dans le dossier courant
#   julia plot_fit.jl /chemin/polynomial.txt       # chemin explicite
#   julia plot_fit.jl /chemin/polynomial.txt /chemin/data_dir
# ─────────────────────────────────────────────────────────────────────────────

# ── Configuration ─────────────────────────────────────────────────────────────
combi = "3_2_2_3_2"
const POLY_PATH = length(ARGS) >= 1 ? ARGS[1] :
                  joinpath(@__DIR__, "../genotoul/output/nint_$(combi)/polynomial.txt") #"../genotoul/output/nint_$(combi)/polynomial.txt")
const DATA_DIR  = length(ARGS) >= 2 ? ARGS[2] :
                  joinpath(@__DIR__, "../../../input/xp_input_20")
const OUT_PATH  = joinpath(@__DIR__, "output/visualise_$(combi).png")

const REPLICATES = ["A", "B", "C"]
const N_CYCLES   = 5

# Paramètres biologiques (fixés)
const r = 0.574619342477644
const K = 6.675449070379925e7
const β = 144.0
const δ = 0.02

println("Polynomial : $POLY_PATH")
println("Data dir   : $DATA_DIR")
println("Output     : $OUT_PATH")


# Couleurs
const COLOR_REPS  = [RGB(0.6, 0.8, 1.0), RGB(31/255, 119/255, 180/255), RGB(0.0, 0.3, 0.7)]
const COLOR_MODEL = RGB(255/255, 127/255, 14/255)
const COLOR_PHI   = RGB(31/255, 119/255, 180/255)

yticks_log(lo, hi) = (
    [10.0^i for i in lo:hi],
    [L"10^{%$i}" for i in lo:hi]
)


# ── Parsing du polynomial.txt → (θ, t_knots) ─────────────────────────────────
# On relit les valeurs aux nœuds depuis la section OPTIMAL KNOTS,
# et on reconstruit un CubicSpline identique à celui de fit_phi.jl.
# C'est la seule façon d'avoir une phi_func autodiff-compatible.
 
function parse_knots(path::String)
    lines = readlines(path)
    t_knots = Float64[]
    θ       = Float64[]
 
    in_knots = false
    for line in lines
        if occursin("OPTIMAL KNOTS", line)
            in_knots = true
            continue
        end
        in_knots || continue
        isempty(strip(line)) && break          # bloc terminé
        occursin(r"^\s*i\s", line) && continue # ligne d'en-tête
 
        # Format : "i    t_i (days)    log(φ(t_i))    φ(t_i)"
        parts = split(strip(line))
        length(parts) >= 3 || continue
        try
            push!(t_knots, parse(Float64, parts[2]))
            push!(θ,       parse(Float64, parts[3]))
        catch
            continue
        end
    end
 
    isempty(t_knots) && error("Aucun nœud trouvé dans $path")
    println("$(length(t_knots)) nœuds lus, t ∈ [$(t_knots[1]), $(t_knots[end])]")
    return θ, t_knots
end
 
θ, t_knots = parse_knots(POLY_PATH)
 
# Reconstruction du spline — exactement comme dans fit_phi.jl
spline   = CubicSpline(θ, t_knots)
phi_func = t -> exp(spline(clamp(t, t_knots[1], t_knots[end])))
 
 
# ── Chargement des données — identique à fit_phi.jl ──────────────────────────
 
struct Cycle
    rep   :: String
    index :: Int
    tH    :: Vector{Float64}
    H     :: Vector{Float64}
    tV    :: Vector{Float64}
    V     :: Vector{Float64}
    u0    :: Vector{Float64}
end
 
function load_replicate(rep::String)::Vector{Cycle}
    cycles     = Cycle[]
    t_offset   = nothing
    t_end_prev = nothing
 
    for cyc in 1:N_CYCLES
        path_H = joinpath(DATA_DIR,
            "hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv")
        path_V = joinpath(DATA_DIR,
            "virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv")
 
        if !isfile(path_H) || !isfile(path_V)
            @warn "Fichier manquant pour rep=$rep cyc=$cyc — ignoré"
            continue
        end
 
        df_H = CSV.read(path_H, DataFrame)
        df_V = CSV.read(path_V, DataFrame)
 
        tH = df_H[:, 1] ./ 24.0
        H  = df_H[:, 2]
        tV = df_V[:, 1] ./ 24.0
        V  = df_V[:, 2]
 
        if t_offset === nothing
            t_offset = tH[1]
        end
        tH .-= t_offset
        tV .-= t_offset
 
        if t_end_prev !== nothing
            gap = tH[1] - t_end_prev
            tH .-= gap
            tV .-= gap
        end
        t_end_prev = tH[end]
 
        push!(cycles, Cycle(rep, cyc, tH, H, tV, V, [H[1], V[1]]))
    end
    return cycles
end
 
println("\nChargement des données...")
all_cycles     = vcat([load_replicate(rep) for rep in REPLICATES]...)
cycles_per_rep = Dict(rep => filter(c -> c.rep == rep, all_cycles) for rep in REPLICATES)
println("$(length(all_cycles)) cycles chargés")
 
 
# ── Modèle ODE — identique à fit_phi.jl ──────────────────────────────────────
 
function ode_model!(dY, Y, p, t)
    S, Vi = Y[1], Y[2]
    ϕ = p(t)
    dY[1] = r * S * (1 - S/K) - ϕ * S * Vi
    dY[2] = β * ϕ * S * Vi - δ * Vi
end
 
# Intégration identique à integrate_cycle de fit_phi.jl,
# avec u0/t0/t1 explicites (pas de struct Cycle fictif)
function integrate_cycle(u0::Vector{Float64}, t0::Float64, t1::Float64,
                         p; n_points=400)
    sv   = collect(range(t0, t1; length=n_points))
    prob = ODEProblem(ode_model!, u0, (t0, t1), p)
    try
        sol = solve(prob, Rodas5(),
                    reltol=1e-6, abstol=1e-6,
                    saveat=sv,
                    isoutofdomain=(u, p, t) -> any(x -> x < 0, u))
        return SciMLBase.successful_retcode(sol) ? sol : nothing
    catch err
        return nothing
    end
end
 
 
# ── φ(t) sur grille fine pour visualisation ───────────────────────────────────
 
t_global_min = minimum(c.tH[1]   for c in all_cycles)
t_global_max = maximum(c.tH[end] for c in all_cycles)
t_fine    = collect(range(t_global_min, t_global_max; length=2000))
phi_vals  = phi_func.(t_fine)
phi_finite = phi_vals[isfinite.(phi_vals) .& (phi_vals .> 0)]
phi_lo = floor(Int, log10(minimum(phi_finite)))
phi_hi = ceil(Int,  log10(maximum(phi_finite)))
 
 
# ── Figure ────────────────────────────────────────────────────────────────────
 
pl = plot(
    layout         = (3, 1),
    size           = (1200, 1200),
    left_margin    = 15mm,
    right_margin   = 10mm,
    top_margin     = 15mm,
    bottom_margin  = 10mm,
    guidefontsize  = 20,
    tickfontsize   = 18,
    titlefontsize  = 18,
    legendfontsize = 16,
    link           = :x,
)
 
# ── Panneau 1 : φ(t) ──────────────────────────────────────────────────────────
 
plot!(pl[1], t_fine, phi_vals;
    yscale    = :log10,
    color     = COLOR_PHI,
    lw        = 3,
    label     = "",
    ylabel    = L"\varphi(t)\ \mathrm{(mL\ cell^{-1}\ day^{-1})}",
    title     = "Taux d'infection φ(t)",
    yticks    = yticks_log(phi_lo, phi_hi),
    grid      = true,
    gridalpha = 0.3,
)
 
# Nœuds du spline
scatter!(pl[1], t_knots, exp.(θ);
    color=COLOR_MODEL, markersize=6, markerstrokewidth=0, label="nœuds")
 
# Frontières de cycles
for cyc in cycles_per_rep["A"]
    vline!(pl[1], [cyc.tH[1]];
        color=:gray, lw=1.2, ls=:dash, alpha=0.5, label="")
    annotate!(pl[1], cyc.tH[1] + 0.3,
        maximum(phi_finite) * 10^0.5,
        text("C$(cyc.index)", 10, :gray, :left))
end
 
# ── Panneaux 2 & 3 : données expérimentales ───────────────────────────────────
 
for (i, rep) in enumerate(REPLICATES)
    cycs = cycles_per_rep[rep]
    isempty(cycs) && continue
    col = COLOR_REPS[i]
    for cyc in cycs
        scatter!(pl[2], cyc.tH, cyc.H;
            color=col, alpha=0.7, markersize=8, markerstrokewidth=0,
            label=(cyc.index == 1 ? "Rep $rep" : ""),
        )
        scatter!(pl[3], cyc.tV, cyc.V;
            color=col, alpha=0.7, markersize=8, markerstrokewidth=0, label="",
        )
    end
end
 
# ── Une simulation par cycle, u0 = moyenne des réplicats ─────────────────────
 
for cyc_idx in 1:N_CYCLES
    cycs_this = filter(c -> c.index == cyc_idx, all_cycles)
    isempty(cycs_this) && continue
 
    u0_mean = [mean(c.u0[1] for c in cycs_this),
               mean(c.u0[2] for c in cycs_this)]
    t0 = minimum(c.tH[1]                   for c in cycs_this)
    t1 = maximum(max(c.tH[end], c.tV[end]) for c in cycs_this)
 
    sol = integrate_cycle(u0_mean, t0, t1, phi_func)
    if sol !== nothing
        plot!(pl[2], sol.t, max.(sol[1, :], 1e-12);
            color=COLOR_MODEL, lw=3,
            label=(cyc_idx == 1 ? "Modèle" : ""),
        )
        plot!(pl[3], sol.t, max.(sol[2, :], 1e-12);
            color=COLOR_MODEL, lw=3, label="",
        )
    else
        @warn "ODE non résolue pour cycle $cyc_idx"
    end
end
 
# ── Axes ──────────────────────────────────────────────────────────────────────
 
plot!(pl[2];
    yscale    = :log10,
    yticks    = yticks_log(3, 9),
    ylims     = (1e3, 1e9),
    ylabel    = L"S(t)\ \mathrm{(cells\ mL^{-1})}",
    title     = "Hôtes S(t)",
    legend    = :topright,
    grid      = true,
    gridalpha = 0.3,
)
plot!(pl[3];
    yscale    = :log10,
    yticks    = yticks_log(6, 10),
    ylims     = (1e6, 1e10),
    ylabel    = L"V(t)\ \mathrm{(virions\ mL^{-1})}",
    xlabel    = "Temps (jours)",
    title     = "Virus V(t)",
    legend    = :none,
    grid      = true,
    gridalpha = 0.3,
)
 
# ── Sauvegarde ────────────────────────────────────────────────────────────────
 
mkpath(dirname(OUT_PATH))
savefig(pl, OUT_PATH)
println("\nSaved → $OUT_PATH")