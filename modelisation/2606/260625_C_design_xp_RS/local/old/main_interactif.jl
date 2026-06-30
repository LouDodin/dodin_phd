###############################################################
#
#      Estimation de γ — Design expérimental
#      Version Julia / GLMakie
#
###############################################################

using GLMakie
using Makie
using Optim
using Random
using Statistics

###############################################################
# Constantes
###############################################################

const TRUE_r = 0.574619342477644
const TRUE_K = 6.675449070379925e7
const TRUE_alpha = 2.185e-5
const TRUE_gamma = 9.47815e-3

const H0 = 1e5
const prop_R0 = 1.0 - 1e-15

###############################################################
# Structure contenant les résultats
###############################################################

struct MCResult

    curve_t::Vector{Float64}
    curve_true::Vector{Float64}
    curve_fit::Vector{Float64}

    scatter_t::Vector{Float64}
    scatter_pR::Vector{Float64}

    tp_t::Vector{Float64}
    tp_pR::Vector{Float64}

    hist_x::Vector{Float64}
    hist_y::Vector{Int}

    gamma_median::Float64
    err_median::Float64
    err25::Float64
    err75::Float64

    times::Vector{Float64}
    true_pR::Vector{Float64}

end

###############################################################
# Simulation ODE
###############################################################

function simulate(
    t_end::Float64,
    gamma::Float64,
    query_times::Vector{Float64};
    dt=0.5
)

    base_grid = collect(0:dt:t_end)

    t_grid = sort(unique(vcat(
        [0.0],
        base_grid,
        query_times
    )))

    S = (1 - prop_R0) * H0
    R = prop_R0 * H0

    pR_grid = zeros(length(t_grid))

    pR_grid[1] = R/(S+R)

    for i in 2:length(t_grid)

        Δt = t_grid[i]-t_grid[i-1]

        H = S+R

        Snew =
            max(
                0,
                S +
                Δt*(
                    TRUE_r*S*(1-H/TRUE_K)
                    -TRUE_alpha*S
                    +gamma*R
                )
            )

        Rnew =
            max(
                0,
                R +
                Δt*(
                    TRUE_r*R*(1-H/TRUE_K)
                    +TRUE_alpha*S
                    -gamma*R
                )
            )

        S = Snew
        R = Rnew

        pR_grid[i] = R/(S+R)

    end

    pR_query = Float64[]

    for qt in query_times

        idx = findfirst(x->isapprox(x,qt,atol=1e-10),t_grid)

        push!(pR_query,pR_grid[idx])

    end

    return t_grid,pR_grid,pR_query

end

###############################################################
# Construction des timepoints
###############################################################

function build_timepoints(
    n_groups::Int,
    Tmax::Float64
)

    collect(
        range(
            Tmax/n_groups,
            Tmax,
            length=n_groups
        )
    )

end

###############################################################
# Utilitaire
###############################################################

function relative_error(est,trueval)

    abs(est-trueval)/trueval*100

end

###############################################################
# Couleurs de l'interface
###############################################################

const BG = RGBf(0.06,0.09,0.16)
const PANEL = RGBf(0.10,0.15,0.27)

const BLUE = RGBf(0.23,0.73,0.98)
const ORANGE = RGBf(0.98,0.60,0.24)
const GREEN = RGBf(0.20,0.83,0.60)
const PINK = RGBf(0.95,0.45,0.70)
const YELLOW = RGBf(0.95,0.75,0.20)
const GREY = RGBf(0.55,0.60,0.65)

###############################################################
# Paramètres par défaut
###############################################################

const DEFAULT_GROUPS = 5
const DEFAULT_TMAX = 50.0
const DEFAULT_CELLS = 96
const DEFAULT_REP = 3
const DEFAULT_MC = 30


###############################################################
# Tirage binomial
###############################################################

"""
Renvoie une proportion binomiale.

Exemple :
96 cellules
p = 0.73

→ retourne par exemple 0.71875
"""
function sample_binomial(n::Int, p::Float64)

    c = 0

    @inbounds for i in 1:n

        if rand() < p
            c += 1
        end

    end

    return c / n

end


###############################################################
# Fonction objectif
###############################################################

function objective_gamma(
    γ,
    times,
    observations
)

    _, _, pred = simulate(
        maximum(times),
        γ,
        times
    )

    s = 0.0

    @inbounds for i in eachindex(pred)

        d = observations[i] - pred[i]

        s += d*d

    end

    return s

end


###############################################################
# Ajustement de γ
###############################################################

function fit_gamma(
    times,
    observations
)

    result = optimize(

        γ -> objective_gamma(
            γ,
            times,
            observations
        ),

        1e-4,
        0.5,

        Brent()

    )

    return Optim.minimizer(result)

end


###############################################################
# Histogramme
###############################################################

function build_histogram(values; bins=12)

    m = maximum(values)

    m = max(m,1.0)

    width = m/bins

    xs = Float64[]
    ys = Int[]

    for i in 1:bins

        left = (i-1)*width
        right = i*width

        push!(xs,(left+right)/2)

        push!(
            ys,
            count(
                x -> left <= x < right,
                values
            )
        )

    end

    return xs,ys

end


###############################################################
# Courbes pour affichage
###############################################################

function build_curves(
    t_grid,
    pR_true,
    t_fit,
    pR_fit
)

    n = length(t_grid)

    step = max(1,fld(n,120))

    curve_t = Float64[]
    curve_true = Float64[]
    curve_fit = Float64[]

    for i in 1:step:n

        push!(curve_t,t_grid[i])

        push!(curve_true,pR_true[i])

        idx = round(
            Int,
            (i-1)/(n-1)*(length(t_fit)-1)+1
        )

        idx = clamp(
            idx,
            1,
            length(pR_fit)
        )

        push!(
            curve_fit,
            pR_fit[idx]
        )

    end

    return curve_t,curve_true,curve_fit

end


###############################################################
# Monte Carlo complet
###############################################################

function runMonteCarlo(

    N_GROUPS::Int,
    T_MAX::Float64,
    N_CELLS::Int,
    N_REP::Int,
    N_MC::Int=30

)

    Random.seed!(42)

    ###########################################################
    # Timepoints
    ###########################################################

    times = build_timepoints(
        N_GROUPS,
        T_MAX
    )

    ###########################################################
    # Courbe vraie
    ###########################################################

    t_grid,
    pR_grid,
    true_pR =
        simulate(
            T_MAX,
            TRUE_gamma,
            times
        )

    ###########################################################
    # Résultats MC
    ###########################################################

    gamma_est = Float64[]

    gamma_error = Float64[]

    scatter_t = Float64[]
    scatter_pR = Float64[]

    ###########################################################
    # Boucle Monte Carlo
    ###########################################################

    for mc in 1:N_MC

        observations = Float64[]

        #######################################################

        for (k,p) in enumerate(true_pR)

            reps = Float64[]

            for r in 1:N_REP

                value = sample_binomial(
                    N_CELLS,
                    p
                )

                push!(reps,value)

                push!(scatter_t,times[k])
                push!(scatter_pR,value)

            end

            push!(
                observations,
                mean(reps)
            )

        end

        #######################################################

        γ = fit_gamma(
            times,
            observations
        )

        push!(
            gamma_est,
            γ
        )

        push!(
            gamma_error,
            relative_error(
                γ,
                TRUE_gamma
            )
        )

    end

    ###########################################################
    # Statistiques
    ###########################################################

    err = sort(gamma_error)

    γsort = sort(gamma_est)

    γmedian = γsort[cld(length(γsort),2)]

    errMedian = median(err)

    err25 = quantile(err,0.25)

    err75 = quantile(err,0.75)

    ###########################################################
    # Courbe ajustée
    ###########################################################

    t_fit,
    pR_fit,
    _ =
        simulate(
            T_MAX,
            γmedian,
            times
        )

    curve_t,
    curve_true,
    curve_fit =
        build_curves(
            t_grid,
            pR_grid,
            t_fit,
            pR_fit
        )

    ###########################################################
    # Histogramme
    ###########################################################

    hist_x,
    hist_y =
        build_histogram(err)

    ###########################################################
    # Résultat
    ###########################################################

    return MCResult(

        curve_t,
        curve_true,
        curve_fit,

        scatter_t,
        scatter_pR,

        times,
        true_pR,

        hist_x,
        hist_y,

        γmedian,
        errMedian,
        err25,
        err75,

        times,
        true_pR

    )

end


###############################################################
# INTERFACE SIMPLIFIÉE
###############################################################

nGroups = Observable(DEFAULT_GROUPS)
tMax    = Observable(DEFAULT_TMAX)
nCells  = Observable(DEFAULT_CELLS)
nRep    = Observable(DEFAULT_REP)

result = Observable(runMonteCarlo(
    DEFAULT_GROUPS,
    DEFAULT_TMAX,
    DEFAULT_CELLS,
    DEFAULT_REP,
    DEFAULT_MC
))

###############################################################
# FIGURE
###############################################################

fig = Figure(size = (1200, 3000))

controls = fig[1, 1] = GridLayout(width = 250)
plots    = fig[1, 2] = GridLayout()

###############################################################
# SLIDERS
###############################################################

Label(controls[1,1], "Groups")
sGroups = Slider(controls[2,1], range=2:10, startvalue=DEFAULT_GROUPS)

Label(controls[3,1], "T max")
sTmax = Slider(controls[4,1], range=10:5:120, startvalue=DEFAULT_TMAX)

Label(controls[5,1], "Cells")
sCells = Slider(controls[6,1], range=12:12:384, startvalue=DEFAULT_CELLS)

Label(controls[7,1], "Rep")
sRep = Slider(controls[8,1], range=1:6, startvalue=DEFAULT_REP)

###############################################################
# AXES (2 PLOTS)
###############################################################

ax1 = Axis(plots[1, 1],
    xlabel="Temps",
    ylabel="pR(t)"
)

ax2 = Axis(plots[2, 1],
    xlabel="Erreur (%)",
    ylabel="Fréquence"
)

###############################################################
# OBSERVABLES
###############################################################

trueLine = Observable(Point2f[])
fitLine  = Observable(Point2f[])
scatterP = Observable(Point2f[])
histBars = Observable(Point2f[])
medianLine = Observable(0.0)

###############################################################
# PLOT 1 + LEGENDES
###############################################################

lines!(ax1, trueLine, label="Simulation")
lines!(ax1, fitLine, label="Fit γ")
scatter!(ax1, scatterP, markersize=3, label="Observations")

axislegend(ax1)

###############################################################
# PLOT 2 + HIST + MÉDIANE
###############################################################

barplot!(ax2, histBars, label="Erreur MC")

vlines!(ax2, medianLine,
    label="Médiane erreur",
    linewidth=2
)

axislegend(ax2)

###############################################################
# UPDATE
###############################################################

function update!()

    res = runMonteCarlo(
        Int(sGroups.value[]),
        Float64(sTmax.value[]),
        Int(sCells.value[]),
        Int(sRep.value[]),
        DEFAULT_MC
    )

    result[] = res

    trueLine[] = Point2f.(res.curve_t, res.curve_true)
    fitLine[]  = Point2f.(res.curve_t, res.curve_fit)
    scatterP[] = Point2f.(res.scatter_t, res.scatter_pR)
    histBars[] = Point2f.(res.hist_x, res.hist_y)

    medianLine[] = res.err_median

    autolimits!(ax1)
    autolimits!(ax2)
end

###############################################################
# CALLBACKS
###############################################################

on(sGroups.value) do _ update!() end
on(sTmax.value)   do _ update!() end
on(sCells.value)  do _ update!() end
on(sRep.value)    do _ update!() end

###############################################################
# INIT
###############################################################

update!()

display(fig)