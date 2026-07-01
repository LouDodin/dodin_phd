## =====================================================================
## Burst size — méthode 2 : ajustement d'un modèle ODE (SIV), par cycle
##
## Principe : pour chaque (réplique, cycle), on ajuste un modèle SIV
## simple aux données observées (H = hôtes totaux, V = virus libres),
## en estimant le burst size β (et les autres paramètres) SPECIFIQUEMENT
## pour ce cycle. On obtient ainsi une série β(cycle) ≈ β(temps), que
## l'on compare à la méthode directe.
##
## Modèle SIV classique :
##   dS/dt = -phi*S*V
##   dI/dt =  phi*S*V - delta*I
##   dV/dt =  beta*delta*I - phi*S*V - mu*V
##
## avec :
##   S     = hôtes sains (susceptibles)
##   I     = hôtes infectés (pas encore lysés)
##   V     = virus libres
##   phi   = taux d'adsorption (mL/jour)
##   delta = 1/période de latence (jour^-1)
##   beta  = burst size (virions produits par cellule lysée)  <-- paramètre d'intérêt
##   mu    = taux de décroissance des virus libres (jour^-1)
##
## H_observé (cytométrie, ne distingue pas S et I) = S + I
##
## Adapter les bornes de paramètres à votre système biologique.
## =====================================================================

using DifferentialEquations
using OrdinaryDiffEqRosenbrock
using CSV
using DataFrames
using Statistics
using Plots
using Measures
using BlackBoxOptim
using Random



const replicates  = ["A", "B", "C"]
const cycles_sim  = 5

output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)


## ===== Input: S & V data =====

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

# Check plot
pl_data = plot(layout=(1,2), size=(900,350), margins=5mm, legend=:bottomright)

for rep in replicates
    for cyc in 1:cycles_sim
        data = raw_data[rep][cyc]
        scatter!(pl_data[1], data.tH, data.H;
            color=replicate_colors[rep], label=cyc==1 ? "Rep $rep" : "",
            xlabel="Time (days)", ylabel="Host abundance (cell/mL)",
            yscale=:log10, ylims=(1e2,1e8), title="H")
        scatter!(pl_data[2], data.tV, data.V;
            color=replicate_colors[rep], label=cyc==1 ? "Rep $rep" : "",
            xlabel="Time (days)", ylabel="Virus abundance (part/mL)",
            yscale=:log10, ylims=(1e3,1e10), title="Vi")
    end
end
#display(pl_data)




## ----- Définition du modèle -----
function siv_model!(du, u, p, t)
    S, I, V = u
    phi, delta, beta, mu = p

    infection = phi * S * V
    lysis     = delta * I

    du[1] = -infection
    du[2] =  infection - lysis
    du[3] =  beta * lysis - infection - mu * V
end

const isoutofdomain_siv = (u, p, t) -> any(x -> x < 0, u)

## ----- Bornes des paramètres (à ajuster selon votre système) -----
## ordre : [phi, delta, beta, mu]
const param_names = ["phi", "delta", "beta", "mu"]
const lower_b = log.([1e-10, 0.5,  1.0,  0.01])
const upper_b = log.([1e-6,  10.0, 2000.0, 5.0])

## ----- Fonction de simulation -----
function simulate_siv(p_log::Vector{Float64}, S0::Float64, I0::Float64, V0::Float64,
                       tspan::Tuple{Float64,Float64}, t_eval::Vector{Float64})
    p = exp.(p_log)
    u0 = [S0, I0, V0]
    prob = ODEProblem(siv_model!, u0, tspan, p)
    sol = solve(prob, Rodas5P(); isoutofdomain=isoutofdomain_siv,
                saveat=t_eval, reltol=1e-8, abstol=1e-10, maxiters=1e6)
    return sol
end

## ----- Fonction coût : écart log-résidus entre simulation et observations -----
## On compare H_sim = S+I à H_obs, et V_sim = V à V_obs, sur échelle log
## (cohérent avec le fait que vos données sont typiquement tracées en log10).
function cost_siv(p_log::Vector{Float64}, tH::Vector{Float64}, H::Vector{Float64},
                   tV::Vector{Float64}, V::Vector{Float64})

    # estimation des CI : S0 ~ H[1] (population quasi saine au début du cycle),
    # I0 ~ 0, V0 ~ V[1]
    S0 = H[1]
    I0 = 1e-3 * H[1]     # petite valeur non nulle pour amorcer l'infection
    V0 = V[1]

    t_all = sort(unique(vcat(tH, tV)))
    tspan = (t_all[1], t_all[end])

    sol = simulate_siv(p_log, S0, I0, V0, tspan, t_all)

    if sol.retcode != SciMLBase.ReturnCode.Success || any(isnan, sol.u[end])
        return 1e10
    end

    # extraire H_sim et V_sim aux temps observés
    idxH = [findfirst(==(t), t_all) for t in tH]
    idxV = [findfirst(==(t), t_all) for t in tV]

    H_sim = [sol.u[i][1] + sol.u[i][2] for i in idxH]
    V_sim = [sol.u[i][3] for i in idxV]

    eps = 1.0  # évite log(0)
    res_H = log10.(H .+ eps) .- log10.(H_sim .+ eps)
    res_V = log10.(V .+ eps) .- log10.(V_sim .+ eps)

    return sum(res_H.^2) + sum(res_V.^2)
end

## ----- Ajustement pour un (rep, cycle) donné -----
"""
    fit_burst_ode(tH, H, tV, V; n_runs=3, max_steps=20000)

Ajuste le modèle SIV à un seul cycle d'infection et retourne :
- p_fit  : vecteur des paramètres ajustés [phi, delta, beta, mu]
- fitness: valeur de la fonction coût au minimum trouvé
"""
function fit_burst_ode(tH::Vector{Float64}, H::Vector{Float64},
                        tV::Vector{Float64}, V::Vector{Float64};
                        n_runs::Int=3, max_steps::Int=20000)

    best_fit_val = Inf
    best_p = nothing

    for run in 1:n_runs
        Random.seed!(run)
        res = bboptimize(p -> cost_siv(p, tH, H, tV, V);
                          SearchRange = collect(zip(lower_b, upper_b)),
                          NumDimensions = 4,
                          Method = :adaptive_de_rand_1_bin_radiuslimited,
                          MaxSteps = max_steps,
                          TraceMode = :silent)

        fit_val = best_fitness(res)
        if fit_val < best_fit_val
            best_fit_val = fit_val
            best_p = exp.(best_candidate(res))
        end
    end

    return best_p, best_fit_val
end

## ----- Boucle sur toutes les répliques / cycles -----
## (réutilise raw_data, replicates, replicate_colors, output_dir)

ode_burst_results = DataFrame(replicate=String[], cycle=Int[], t_cycle_start=Float64[],
                               phi=Float64[], delta=Float64[], beta=Float64[],
                               mu=Float64[], fitness=Float64[])

for rep in replicates
    for data in raw_data[rep]
        println("Fitting rep=$rep cycle=$(data.index) ...")
        p_fit, fit_val = fit_burst_ode(data.tH, data.H, data.tV, data.V;
                                        n_runs=n_runs, max_steps=15000)
        if p_fit !== nothing
            push!(ode_burst_results, (replicate=rep, cycle=data.index,
                                       t_cycle_start=data.tH[1],
                                       phi=p_fit[1], delta=p_fit[2],
                                       beta=p_fit[3], mu=p_fit[4],
                                       fitness=fit_val))
        end
    end
end

CSV.write(joinpath(output_dir, "burst_ODE_per_cycle.csv"), ode_burst_results)

## ----- Visualisation : beta (burst size) estimé par cycle / temps -----
pl_beta = plot(size=(700,450), margins=5mm,
               xlabel="Cycle d'infection", ylabel="Burst size β (virions/cellule)",
               title="Burst size estimé par ajustement ODE (par cycle)",
               legend=:topright)

for rep in replicates
    sub = ode_burst_results[ode_burst_results.replicate .== rep, :]
    sort!(sub, :cycle)
    plot!(pl_beta, sub.cycle, sub.beta;
          color=replicate_colors[rep], label="Rep $rep",
          marker=:circle, markersize=5, lw=2)
end

savefig(pl_beta, joinpath(output_dir, "burst_size_ODE_vs_cycle.pdf"))
display(pl_beta)

## ----- Idem mais en fonction du temps réel (début de cycle) -----
pl_beta_t = plot(size=(700,450), margins=5mm,
                  xlabel="Temps (jours, début de cycle)", ylabel="Burst size β",
                  title="Burst size estimé (ODE) en fonction du temps",
                  legend=:topright)

for rep in replicates
    sub = ode_burst_results[ode_burst_results.replicate .== rep, :]
    sort!(sub, :t_cycle_start)
    plot!(pl_beta_t, sub.t_cycle_start, sub.beta;
          color=replicate_colors[rep], label="Rep $rep",
          marker=:circle, markersize=5, lw=2)
end

savefig(pl_beta_t, joinpath(output_dir, "burst_size_ODE_vs_time.pdf"))
display(pl_beta_t)

println("Burst size (méthode ODE, par cycle) ajusté et exporté dans : $output_dir")