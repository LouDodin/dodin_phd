using Colors
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
using LaTeXStrings

## ===== COLORS =====
color_A = RGB(0.6, 0.8, 1.0)                 # light blue
color_B = RGB(31/255, 119/255, 180/255)     # original blue (KEEP THIS ONE)
color_C = RGB(0.0, 0.3, 0.7)                # dark blue

data_color = RGB(31/255, 119/255, 180/255)            # cycle change lines


## ===== FUNCTION: LOAD + SHIFT CYCLES =====
function load_rep(replicate, cycles)

    tH_all = Float64[]
    tV_all = Float64[]
    H_all  = Float64[]
    V_all  = Float64[]

    cycle_changes = Float64[] # Float64[210.5./24]   # <-- store boundaries here
    prev_end = 0.0

    for c in cycles
        df_H = CSV.read(joinpath(@__DIR__, "input/xp_input_all/hostData_coevoCondition_Temperature20_Replicate$(replicate)_Cycle$c.csv"), DataFrame)
        df_V = CSV.read(joinpath(@__DIR__, "input/xp_input_all/virusData_coevoCondition_Temperature20_Replicate$(replicate)_Cycle$c.csv"), DataFrame)

        tH = df_H[:,1] ./ 24
        tV = df_V[:,1] ./ 24

        H = df_H[:,2]
        V = df_V[:,2]

        # shift time so cycles connect
        shift = tH[1] - prev_end
        tH .-= shift

        #shift2 = tV[1] - prev_end
        tV .-= tV[1]-tH[1]

        append!(tH_all, tH)
        append!(tV_all, tV)
        append!(H_all, H)
        append!(V_all, V)

        prev_end = tH[end]

        # store cycle boundary (end of this cycle)
        push!(cycle_changes, prev_end)
    end
    cycle_changes = cycle_changes[1:end-1]

    return tH_all, H_all, tV_all, V_all, cycle_changes
end






## ===== LOAD EACH REPLICATE =====
t_HA, HA, t_VA, VA, cA = load_rep("A", 1:5)
t_HB, HB, t_VB, VB, cB = load_rep("B", 1:5)
t_HC, HC, t_VC, VC, cC = load_rep("C", 1:5)

cycle_changes = unique(sort(vcat(cA, cB, cC)))


## ===== PLOT =====
ytick_vals1   = [10.0^i for i in 2:8]
ytick_labels1 = [L"10^{%$i}" for i in 2:8]

ytick_vals2   = [10.0^i for i in 3:9]
ytick_labels2 = [L"10^{%$i}" for i in 3:9]


pl_sim = plot(
    layout=(2,1),
    size=(1000,1100),

    left_margin=15mm,
    right_margin=10mm,
    top_margin=10mm,
    bottom_margin=10mm,

    grid=true,
    yscale=:log10,

    xlims=(0, 15),

    ytickfontsize=26,
    legendfontsize=20,
    guidefontsize=24,
    xtickfontsize=24,
    titlefontsize=24,
    ylabel="Time (days)"
)


## ===== HOST =====
plot!(pl_sim[1], t_HA, HA,
    label="Replicate A", lw=4,
    color=color_A, markershape=:circle, markersize=5,
    alpha=0.7,
    ylabel="Host abundance\n(cell/mL)",
    legend=:bottomleft,
    ylims=(1e2, 1e8),
    yticks=(ytick_vals1, ytick_labels1)
)

plot!(pl_sim[1], t_HB, HB,
    label="Replicate B",
    lw=4,
    color=color_B, markershape=:circle, markersize=6,
    alpha=0.7
)

plot!(pl_sim[1], t_HC, HC,
    label="Replicate C",
    lw=4,
    color=color_C, markershape=:circle, markersize=5,
    alpha=0.7
)


## ===== VIRUS =====
plot!(pl_sim[2], t_VA, VA,
    label="Replicate A",
    lw=4,
    color=color_A, markershape=:square, markersize=5,
    alpha=0.7,
    ylabel="Virus abundance\n(part/mL)",
    legend=:bottomleft,
    ylims=(1e3, 1e9),
    yticks=(ytick_vals2, ytick_labels2)
)

plot!(pl_sim[2], t_VB, VB,
    label="Replicate B",
    lw=4,
    color=color_B, markershape=:square, markersize=6,
    alpha=0.7
)

plot!(pl_sim[2], t_VC, VC,
    label="Replicate C",
    lw=4,
    color=color_C, markershape=:square, markersize=5,
    alpha=0.7
)


## ===== CYCLE CHANGE VLINES =====
for t_change in cycle_changes
    for i in 1:2
        vline!(pl_sim[i], [t_change],
            color=data_color,
            linestyle=:dot,
            lw=2,
            label = t_change == cycle_changes[1] ? "Dilution" : nothing
        )
    end
end

display(pl_sim)