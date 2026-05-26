## ===== Package import =====
using Dates
using CSV
using DataFrames
using Optim
using DifferentialEquations
using LogExpFunctions
using Measures
using Plots
using Statistics

t_global_start = now()

# ┌─────────────────────────────────────────────────────────────────────┐
# │  Mettre RUN_FIT = true  → lance l'optimisation                      │
# │  Mettre RUN_FIT = false → utilise les valeurs directement           |
# └─────────────────────────────────────────────────────────────────────┘
const RUN_FIT_1 = false
const RUN_FIT_OTHERS = true


## ===== Imports =====

# --- Data import ---
dir_input = joinpath(@__DIR__, "input/xp_input_20°")

t_H = Vector{Vector{Float64}}()
H   = Vector{Vector{Float64}}()
t_V = Vector{Vector{Float64}}()
V   = Vector{Vector{Float64}}()

for data in ("host", "virus"), cycle in 1:5, rep in ("A", "B", "C")
    df = CSV.read(joinpath(dir_input, "$(data)Data_coevoCondition_Temperature20_Replicate$(rep)_cycle$(cycle).csv"), DataFrame)
    t = collect(skipmissing(df[:, 1])) ./ 24
    x = collect(skipmissing(df[:, 2]))
    if data == "host"
        push!(t_H, t); push!(H, x)
    else
        push!(t_V, t); push!(V, x)
    end
end

t_all_rep    = [sort(unique(vcat(t_H[i], t_V[i]))) for i in eachindex(t_H)]
t_end_cycle1 = t_all_rep[1][end]
Y0_global    = log.([mean([first(H[1]), first(H[2]), first(H[3])]), 1e-6, 1e-6, mean([first(V[1]), first(V[2]), first(V[3])]), 1e-6, 1e-6, 1e-6])

# --- Model import ---
include(joinpath(@__DIR__, "models.jl"))
model1 = MODELS["SIVi_mu"]
model2 = MODELS["SIVi_mu_2"]


## ===== Constants =====
const K = 9.784708604680645e7


## ===== KNOWN PARAMS =====
const KNOWN_PARAMS = [0.538001865245161, 1.0814785163571014e-8, 69.99999999999996, 0.00010000000000000026, 39.99999999999998]
const KNOWN_PARAMS_OTHERS = [1e-100]


## ===== Utils =====
function reconstruct(
    logtheta :: AbstractVector{<:Real},
    model    :: ModelSpec,
    fixed    :: Dict{Symbol,Float64} = Dict{Symbol,Float64}()
)
    fitted = Dict{Symbol,Float64}()
    idx = 1
    for p in model.fit_params
        if haskey(fixed, p)
            fitted[p] = fixed[p]
        else
            fitted[p] = exp(logtheta[idx])
            idx += 1
        end
    end
    return Float64[
        p == :k          ? K        :
        haskey(fixed, p) ? fixed[p] :
        get(fitted, p, 0.0)
        for p in model.full_params
    ]
end

function extract_all(sol, tvec)
    S    = [exp(sol(t)[1]) for t in tvec]
    I    = [exp(sol(t)[2]) for t in tvec]
    R    = [exp(sol(t)[3]) for t in tvec]
    Vi   = [exp(sol(t)[4]) for t in tvec]
    Vdp  = [exp(sol(t)[5]) for t in tvec]
    Vdip = [exp(sol(t)[6]) for t in tvec]
    Ev   = [exp(sol(t)[7]) for t in tvec]
    return S, I, R, Vi, Vdp, Vdip, Ev
end

function compute_error(sol, j, t_cut, t_H, H, t_V, V)
    err = 0.0
    idx_H = findall(t -> t <= t_cut, t_H[j])
    if !isempty(idx_H)
        H_pred = [logsumexp(sol(t)[1:3]) for t in t_H[j][idx_H]]
        err += sum((H_pred .- log.(H[j][idx_H])).^2)
    end
    idx_V = findall(t -> t <= t_cut, t_V[j])
    if !isempty(idx_V)
        V_pred = [logsumexp(sol(t)[4:6]) for t in t_V[j][idx_V]]
        err += sum((V_pred .- log.(V[j][idx_V])).^2)
    end
    return err
end

function compute_residue(sol, j, t1_loc, t2_loc, t_H, H, t_V, V)
    err = 0.0
    idx_H = findall(t -> t > t1_loc && t <= t2_loc, t_H[j])
    if !isempty(idx_H)
        H_pred = [logsumexp(sol(t)[1:3]) for t in t_H[j][idx_H]]
        err += sum((H_pred .- log.(H[j][idx_H])).^2)
    end
    idx_V = findall(t -> t > t1_loc && t <= t2_loc, t_V[j])
    if !isempty(idx_V)
        V_pred = [logsumexp(sol(t)[4:6]) for t in t_V[j][idx_V]]
        err += sum((V_pred .- log.(V[j][idx_V])).^2)
    end
    return err
end

function compute_derivative(t, y)
    n = length(y)
    dy = zeros(n-1)

    for i in 1:n-1
        dt = t[i+1] - t[i]
        dy[i] = (y[i+1] - y[i]) / dt
    end

    return dy
end

function detect_breakpoint(t, residue_local; factor=5.0, window=5)
    dresidue = compute_derivative(t, residue_local)
    # baseline robuste
    base = median(dresidue[1:min(window, length(dresidue))])
    threshold = base * factor
    idx = findfirst(x -> x > threshold, dresidue)
    if isnothing(idx)
        return t[end]
    else
        return t[idx + 1]
    end
end

# --- Colors ---
host_palette  = cgrad([:darkgreen, :chartreuse])
virus_palette = cgrad([:darkred, :orangered])
col_S   = host_palette[0.0];  col_I   = host_palette[0.5];  col_R    = host_palette[1.0]
col_Vi  = virus_palette[0.0]; col_Vdp = virus_palette[0.5]; col_Vdip = virus_palette[1.0]
col_Ev  = :blue


## ===== FIT 1 =====

theta0_fit1 = log.([0.5,   1.3e-7, 50.0, 1e-3, 10.0])
lb_fit1     = log.([0.4,   1e-7,   45.0, 1e-5,  1.0])
ub_fit1     = log.([0.7,   1e-6,   70.0, 1e-2, 40.0])

function obj_fit1(theta)
    p   = reconstruct(theta, model1)
    err = 0.0
    for i in 1:3
        t_start = t_all_rep[i][1]
        prob = ODEProblem(model1.dynamics!, Y0_global, (t_start, t_end_cycle1), p)
        sol  = solve(prob, Rodas5())
        err += compute_error(sol, i, t_end_cycle1, t_H, H, t_V, V)
    end
    return err
end

if RUN_FIT_1
    println("Fit 1 - optimisation en cours...")
    res_fit1   = optimize(obj_fit1, lb_fit1, ub_fit1, theta0_fit1, Fminbox(BFGS()))
    theta_opt1 = exp.(Optim.minimizer(res_fit1))
    println("Optimisation terminee.")
else
    println("Fit 1 - parametres charges manuellement (RUN_FIT_1 = false).")
    theta_opt1 = KNOWN_PARAMS
    res_fit1   = optimize(obj_fit1, lb_fit1, ub_fit1, theta0_fit1, Fminbox(BFGS()),
                          Optim.Options(iterations = 0))
end

println("mu=$(theta_opt1[1]), phi=$(theta_opt1[2]), beta=$(theta_opt1[3]), delta=$(theta_opt1[4]), eta=$(theta_opt1[5])")


## ===== Simulate fit 1 =====
t_start1     = t_all_rep[1][1]
t_all_sorted = t_all_rep[1]
t1           = range(t_start1, t_end_cycle1, length=500)

p_fit1    = reconstruct(log.(theta_opt1), model1)
prob_fit1 = ODEProblem(model1.dynamics!, Y0_global, (t_start1, t_end_cycle1), p_fit1)
sol_fit1  = solve(prob_fit1, Rodas5(); saveat=t1)

S1,I1,R1,Vi1,Vdp1,Vdip1,Ev1 = extract_all(sol_fit1, t1)
H1 = S1 .+ I1 .+ R1
V1 = Vi1 .+ Vdp1 .+ Vdip1


## ===== Residue for fit 1 =====
residue1 = Float64[]
for i in 2:length(t_all_sorted)
    local t_cut = t_all_sorted[i]
    err   = sum(compute_residue(sol_fit1, j, 0.0, t_cut, t_H, H, t_V, V) for j in 1:3)
    push!(residue1, err)
end

threshold1   = 20.0
i_threshold1 = findlast(<(threshold1), residue1)
t_threshold1 = t_all_sorted[i_threshold1 + 1]
println("Fit 1 threshold = $threshold1  =>  t_threshold = $t_threshold1")


## ===== Iterative fits =====

theta_fixed = Dict{Symbol,Float64}(:μ => theta_opt1[1], :β => theta_opt1[3], :δ => theta_opt1[4], :η => theta_opt1[5])

mask1    = t1 .<= t_threshold1
idx_end1 = findlast(mask1)

struct FitResult
    fit_id    :: Int                   # identifiant du fit (segment)
    t_start   :: Float64               # temps de début du segment
    t_end     :: Float64               # temps de fin du segment
    tvec      :: Vector{Float64}       # vecteur des temps simulés sur le segment

    S    :: Vector{Float64}      
    I    :: Vector{Float64}           
    R    :: Vector{Float64}         
    Vi   :: Vector{Float64}       
    Vdp  :: Vector{Float64}           
    Vdip :: Vector{Float64}            
    Ev   :: Vector{Float64}       
    H    :: Vector{Float64}      
    V    :: Vector{Float64}      

    t_residue :: Vector{Float64}       # temps associé aux résidus
    residue   :: Vector{Float64}       # erreur de fit sur chaque intervalle de données
    dresidue      :: Vector{Float64}       # dérivée de la résiduelle (variation du résidu)
    
    theta_opt :: Vector{Float64}       # paramètres optimisés du modèle sur ce segment
end

fit_results = FitResult[]

theta0_list = log.([9.5e-9, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20])
lb_list = log.([9e-9, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25])
ub_list = log.([1e-8, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18])
threshold_list   = [20.0, 4000.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0]

global Y0_cur = log.([S1[idx_end1], I1[idx_end1], R1[idx_end1], Vi1[idx_end1], Vdp1[idx_end1], Vdip1[idx_end1], Ev1[idx_end1]])
global t_start_cur = t_threshold1
global fit_id_cur  = 2
global param_idx   = 1          # index des paramètres à fitter
global data_offset = 0          # décalage global dans t_all_sorted

while t_start_cur < t_end_cycle1 && param_idx <= length(theta0_list)

    println("\n--- Fit $fit_id_cur  t_start = $t_start_cur ---")

    if RUN_FIT_OTHERS
        theta0_fi = [theta0_list[param_idx]]
        lb_fi     = [lb_list[param_idx]] 
        ub_fi     = [ub_list[param_idx]]

        Y0_loc      = copy(Y0_cur)
        t_start_loc = t_start_cur

        obj_fi = function(theta)
            p   = reconstruct(theta, model2, theta_fixed)
            prob = ODEProblem(model2.dynamics!, Y0_loc, (t_start_loc, t_end_cycle1), p)
            sol  = solve(prob, Rodas5())

            err = 0.0
            for j in 1:3
                err += compute_error(sol, j, t_end_cycle1, t_H, H, t_V, V)
            end

            println(err)
            return err
        end

        res_fi       = optimize(obj_fi, lb_fi, ub_fi, theta0_fi, Fminbox(BFGS()))
        theta_opt_fi = exp.(Optim.minimizer(res_fi))
    else
        theta_opt_fi = KNOWN_PARAMS_OTHERS
    end

    println("  phi_i_$(fit_id_cur) = $(theta_opt_fi[1])")

    # ===== Simulation =====
    t_fi_dense = range(t_start_cur, t_end_cycle1, length=500)

    p_fi    = reconstruct(log.(theta_opt_fi), model2, theta_fixed)
    prob_fi = ODEProblem(model2.dynamics!, Y0_cur, (t_start_cur, t_end_cycle1), p_fi)
    sol_fi  = solve(prob_fi, Rodas5(); saveat=t_fi_dense)

    S_fi,I_fi,R_fi,Vi_fi,Vdp_fi,Vdip_fi,Ev_fi = extract_all(sol_fi, t_fi_dense)
    H_fi = S_fi .+ I_fi .+ R_fi
    V_fi = Vi_fi .+ Vdp_fi .+ Vdip_fi

    # ===== Résidus locaux =====
    idx_start = findfirst(t -> t >= t_start_cur, t_all_sorted)

    if idx_start == 1
        t_fi_data = t_all_sorted
    else
        t_fi_data = t_all_sorted[idx_start-1:end]
    end

    residue_local = Float64[]
    for i in 2:length(t_fi_data)
        t_lo = t_fi_data[i-1]
        t_hi = t_fi_data[i]

        err = sum(compute_residue(sol_fi, j, t_lo, t_hi, t_H, H, t_V, V) for j in 1:3)
        push!(residue_local, err)
    end
    t_residue_local = t_fi_data[2:end]
    dresidue = compute_derivative(t_residue_local, residue_local)

    println(residue_local)

    # ===== Détection du seuil (corrigée) =====
    threshold = threshold_list[param_idx]

    idx_local = findlast(x -> x < threshold, residue_local)

    if isnothing(idx_local)
        t_end_cur = t_all_sorted[end]
    else
        t_end_cur = t_fi_data[idx_local + 1]
    end

    println("  t_end = $t_end_cur")

    # ===== Extraction segment =====
    idx_end_dense = findfirst(t -> t >= t_end_cur, t_fi_dense)
    if isnothing(idx_end_dense)
        idx_end_dense = length(t_fi_dense)
    end

    t_plot = collect(t_fi_dense[1:idx_end_dense])
    n = length(t_plot)

    push!(fit_results, FitResult(
        fit_id_cur,
        t_start_cur, t_end_cur,
        t_plot,
        S_fi[1:n], I_fi[1:n], R_fi[1:n],
        Vi_fi[1:n], Vdp_fi[1:n], Vdip_fi[1:n], Ev_fi[1:n],
        H_fi[1:n], V_fi[1:n],
        t_residue_local,
        residue_local,
        dresidue,
        theta_opt_fi
    ))

    # ===== Mise à jour =====
    global Y0_cur = log.([S_fi[idx_end_dense], I_fi[idx_end_dense], R_fi[idx_end_dense], Vi_fi[idx_end_dense], Vdp_fi[idx_end_dense], Vdip_fi[idx_end_dense], Ev_fi[idx_end_dense]])
    global t_start_cur = t_end_cur
    global fit_id_cur += 1
    global param_idx  += 1

    t_start_cur >= t_end_cycle1 && break
end

n_iter = length(fit_results)
println("\nTotal iterative fits: $n_iter  (fit 1 + $n_iter subsequent)")


## ===== Plotting =====

n_rows = 1 + n_iter + 1
fig_h  = 800 * n_rows
fig    = plot(layout=(n_rows, 3), grid=true, size=(2200, fig_h),
              xlabel="Time (days)", legend=:bottomright, margins=15mm)

# Global y limits
all_vals_global = vcat(
    H1, V1, S1, I1, R1, Vi1, Vdp1, Vdip1, Ev1,
    vcat(H...), vcat(V...),
    [vcat(fr.H, fr.V, fr.S, fr.I, fr.R, fr.Vi, fr.Vdp, fr.Vdip, fr.Ev)
     for fr in fit_results]...
)


y_min_g = max(1e-5, minimum(filter(x -> isfinite(x) && x > 0, all_vals_global)))
y_max_g = maximum(filter(x -> isfinite(x) && x > 0, all_vals_global))

t1_vec = collect(t1[mask1])

# ── Row 1: Fit 1 ─────────────────────────────────────────────────────
global row = 1
for i in 1:length(H)
    scatter!(fig[(row-1)*3+1], t_H[i], H[i], color=:green, marker=:circle, alpha=0.5,
             label=(i==1 ? "H data" : false))
    scatter!(fig[(row-1)*3+1], t_V[i], V[i], color=:red,   marker=:square, alpha=0.5,
             label=(i==1 ? "V data" : false))
end
plot!(fig[(row-1)*3+1], t1_vec, H1[mask1], lw=3, color=:green, label="H fit 1", yscale=:log10)
plot!(fig[(row-1)*3+1], t1_vec, V1[mask1], lw=3, color=:red,   label="V fit 1")
ylabel!(fig[(row-1)*3+1], "Concentration (parts/mL)")
title!(fig[(row-1)*3+1],
       "Fit 1 - H and V  [$(round(t_start1,digits=2)) to $(round(t_threshold1,digits=2)) d]")
ylims!(fig[(row-1)*3+1], (y_min_g, y_max_g))

plot!(fig[(row-1)*3+2], t1_vec, S1[mask1],    lw=3, label="S",    color=col_S,    yscale=:log10)
plot!(fig[(row-1)*3+2], t1_vec, I1[mask1],    lw=3, label="I",    color=col_I)
plot!(fig[(row-1)*3+2], t1_vec, R1[mask1],    lw=3, label="R",    color=col_R)
plot!(fig[(row-1)*3+2], t1_vec, Vi1[mask1],   lw=3, label="Vi",   color=col_Vi)
plot!(fig[(row-1)*3+2], t1_vec, Vdp1[mask1],  lw=3, label="Vdp",  color=col_Vdp)
plot!(fig[(row-1)*3+2], t1_vec, Vdip1[mask1], lw=3, label="Vdip", color=col_Vdip)
plot!(fig[(row-1)*3+2], t1_vec, Ev1[mask1],   lw=3, label="Ev",   color=col_Ev)
ylabel!(fig[(row-1)*3+2], "Concentration (parts/mL)")
title!(fig[(row-1)*3+2], "Fit 1 - S,I,R,Vi,Vdp,Vdip,Ev")
ylims!(fig[(row-1)*3+2], (y_min_g, y_max_g))

t_residue1 = [(t_all_sorted[i] + t_all_sorted[i-1]) / 2 for i in 2:length(t_all_sorted)]
plot!(fig[(row-1)*3+3], t_residue1, residue1, lw=3, label="Residue", color=:black)
hline!(fig[(row-1)*3+3], [threshold1],   lw=1, ls=:dash, color=:grey, label="threshold")
vline!(fig[(row-1)*3+3], [t_threshold1], lw=1, ls=:dot,  color=:blue, label="cut")
ylabel!(fig[(row-1)*3+3], "Residue")
title!(fig[(row-1)*3+3], "Fit 1 - Residue")

# ── Rows 2..N: iterative fits ─────────────────────────────────────────
for (k, fr) in enumerate(fit_results)
    global row = 1 + k
    lab = "Fit $(fr.fit_id)"

    for i in 1:length(H)
        scatter!(fig[(row-1)*3+1], t_H[i], H[i], color=:green, marker=:circle, alpha=0.5,
                 label=(i==1 ? "H data" : false))
        scatter!(fig[(row-1)*3+1], t_V[i], V[i], color=:red,   marker=:square, alpha=0.5,
                 label=(i==1 ? "V data" : false))
    end
    plot!(fig[(row-1)*3+1], fr.tvec, fr.H, lw=3, color=:green, label="H $lab", yscale=:log10)
    plot!(fig[(row-1)*3+1], fr.tvec, fr.V, lw=3, color=:red,   label="V $lab")
    ylabel!(fig[(row-1)*3+1], "Concentration (parts/mL)")
    title!(fig[(row-1)*3+1],
           "$lab - H and V  [$(round(fr.t_start,digits=2)) to $(round(fr.t_end,digits=2)) d]")
    ylims!(fig[(row-1)*3+1], (y_min_g, y_max_g))
    xlims!(fig[(row-1)*3+1], (0.0, t1[end]))

    plot!(fig[(row-1)*3+2], fr.tvec, fr.S,    lw=3, label="S",    color=col_S,    yscale=:log10)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.I,    lw=3, label="I",    color=col_I)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.R,    lw=3, label="R",    color=col_R)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.Vi,   lw=3, label="Vi",   color=col_Vi)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.Vdp,  lw=3, label="Vdp",  color=col_Vdp)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.Vdip, lw=3, label="Vdip", color=col_Vdip)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.Ev,   lw=3, label="Ev",   color=col_Ev)
    ylabel!(fig[(row-1)*3+2], "Concentration (parts/mL)")
    title!(fig[(row-1)*3+2],
           "$lab - S,I,R,...  phi_i=$(round(fr.theta_opt[1], sigdigits=3))")
    ylims!(fig[(row-1)*3+2], (y_min_g, y_max_g))
    xlims!(fig[(row-1)*3+2], (0.0, t1[end]))

    plot!(fig[(row-1)*3+3], fr.t_residue, fr.residue, lw=2, label="residue", color = :black)
    plot!(fig[(row-1)*3+3], fr.t_residue[1:end-1], fr.dresidue, lw=2, ls=:dash, label="d residue/dt", color = :black)
    #hline!(fig[(row-1)*3+3], [thresh_fi], lw=1, ls=:dash, color=:grey, label="threshold")
    vline!(fig[(row-1)*3+3], [fr.tvec[end]],  lw=1, ls=:dot,  color=:blue, label="cut")
    ylabel!(fig[(row-1)*3+3], "Residue")
    title!(fig[(row-1)*3+3], "$lab - Residue")
    xlims!(fig[(row-1)*3+3], (0.0, t1[end]))
end

# ── Last row: all fits combined ───────────────────────────────────────
row = n_rows

for i in 1:length(H)
    scatter!(fig[(row-1)*3+1], t_H[i], H[i], color=:green, marker=:circle, alpha=0.5,
             label=(i==1 ? "H data" : false))
    scatter!(fig[(row-1)*3+1], t_V[i], V[i], color=:red,   marker=:square, alpha=0.5,
             label=(i==1 ? "V data" : false))
end
plot!(fig[(row-1)*3+1], t1_vec, H1[mask1], lw=3, color=:green, label="H fit 1", yscale=:log10)
plot!(fig[(row-1)*3+1], t1_vec, V1[mask1], lw=3, color=:red,   label="V fit 1")
plot!(fig[(row-1)*3+2], t1_vec, S1[mask1],    lw=3, label="S fit1",    color=col_S,    yscale=:log10)
plot!(fig[(row-1)*3+2], t1_vec, I1[mask1],    lw=3, label="I fit1",    color=col_I)
plot!(fig[(row-1)*3+2], t1_vec, R1[mask1],    lw=3, label="R fit1",    color=col_R)
plot!(fig[(row-1)*3+2], t1_vec, Vi1[mask1],   lw=3, label="Vi fit1",   color=col_Vi)
plot!(fig[(row-1)*3+2], t1_vec, Vdp1[mask1],  lw=3, label="Vdp fit1",  color=col_Vdp)
plot!(fig[(row-1)*3+2], t1_vec, Vdip1[mask1], lw=3, label="Vdip fit1", color=col_Vdip)
plot!(fig[(row-1)*3+2], t1_vec, Ev1[mask1],   lw=3, label="Ev fit1",   color=col_Ev)

for fr in fit_results
    lab = "fit $(fr.fit_id)"
    plot!(fig[(row-1)*3+1], fr.tvec, fr.H, lw=3, color=:green, label="H $lab", ls=:dash)
    plot!(fig[(row-1)*3+1], fr.tvec, fr.V, lw=3, color=:red,   label="V $lab", ls=:dash)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.S,    lw=2, label=false, color=col_S,    ls=:dash)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.I,    lw=2, label=false, color=col_I,    ls=:dash)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.R,    lw=2, label=false, color=col_R,    ls=:dash)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.Vi,   lw=2, label=false, color=col_Vi,   ls=:dash)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.Vdp,  lw=2, label=false, color=col_Vdp,  ls=:dash)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.Vdip, lw=2, label=false, color=col_Vdip, ls=:dash)
    plot!(fig[(row-1)*3+2], fr.tvec, fr.Ev,   lw=2, label=false, color=col_Ev,   ls=:dash)
    vline!(fig[(row-1)*3+1], [fr.t_start], lw=1, ls=:dot, color=:black, label=false)
    vline!(fig[(row-1)*3+2], [fr.t_start], lw=1, ls=:dot, color=:black, label=false)
end

ylabel!(fig[(row-1)*3+1], "Concentration (parts/mL)")
title!(fig[(row-1)*3+1], "All fits combined - H and V")
ylims!(fig[(row-1)*3+1], (y_min_g, y_max_g))

ylabel!(fig[(row-1)*3+2], "Concentration (parts/mL)")
title!(fig[(row-1)*3+2], "All fits combined - S,I,R,Vi,...")
ylims!(fig[(row-1)*3+2], (y_min_g, y_max_g))

t_residue1 = t_residue1 = t_all_sorted[2:end]
plot!(fig[(row-1)*3+3], t_residue1[1:i_threshold1], residue1[1:i_threshold1], lw=2, label="Fit 1", color=:black)
vline!(fig[(row-1)*3+3], [t_threshold1], lw=1, ls=:dot, color=:black, label=false)
for (i, fr) in enumerate(fit_results)

    # limite droite = début du fit suivant
    if i < length(fit_results)
        t_stop = fit_results[i+1].t_start
    else
        t_stop = fr.t_end
    end

    mask = fr.t_residue .<= t_stop

    plot!(fig[(row-1)*3+3],
        fr.t_residue[mask],
        fr.residue[mask],
        lw=2,
        label="Fit $(fr.fit_id)"
    )

    vline!(fig[(row-1)*3+3], [t_stop],
        lw=1, ls=:dot, color=:black, label=false)
end
ylabel!(fig[(row-1)*3+3], "Residue")
title!(fig[(row-1)*3+3], "All fits - Residue")

display(fig)

run_id    = Dates.format(t_global_start, "yyyymmdd-HHMMSS")
plot_path = joinpath(@__DIR__, "100426_output/$(run_id)_iterative_plot.png")
savefig(fig, plot_path)
println("Plot saved: $plot_path")


## ===== Save log =====
log_path = joinpath(@__DIR__, "100426_output/$(run_id)_log.txt")

open(log_path, "w") do io
    sep = "=" ^ 60
    println(io, sep)
    println(io, "  RUN LOG - iterative fits")
    println(io, "  Date     : $(Dates.format(t_global_start, "yyyy-mm-dd HH:MM:SS"))")
    println(io, "  Mode fit1: $(RUN_FIT_1 ? "optimisation" : "parametres manuels")")
    println(io, sep)

    println(io, "\n--- FIT 1 ---")
    println(io, "  Segment  : $(t_start1) to $(t_threshold1) d")
    println(io, "  Threshold: $threshold1")
    for (i, p) in enumerate(model1.fit_params)
        println(io, "  $(rpad(string(p),8)) = $(round(theta_opt1[i], sigdigits=6))")
    end

    for fr in fit_results
        println(io, "\n--- FIT $(fr.fit_id) ---")
        println(io, "  Segment  : $(fr.t_start) to $(fr.t_end) d")
        println(io, "  phi_i    = $(round(fr.theta_opt[1], sigdigits=6))")
    end

    t_total = (now() - t_global_start).value / 1000
    println(io, "\n$sep")
    println(io, "  Total runtime : $(round(t_total, digits=1))s")
    println(io, "  END OF LOG")
    println(io, sep)
end

println("\n" * "="^60)
println("ALL DONE in $(round((now() - t_global_start).value / 1000, digits=1))s")
println("Log saved: $log_path")
println("Plot saved: $plot_path")
println("="^60 * "\n")